# 交叉波导的伴随法形状优化

用有限元 + 伴随法（adjoint method）自动优化光子交叉波导的交叉口形状，
把直通端透射率从 **|S₂₁|² ≈ 0.71** 提升到 **≈ 0.96**。

设计变量不是每个网格单元的密度，而是交叉口边界上三次 Bézier 曲线的控制点坐标；
目标函数对全部控制点的梯度由**一次**伴随求解得到，与设计变量个数无关。

> 代码写在 2022 年 2 月（见 `git log`），2022 年 3–4 月又用它做了微环阵列的
> 能带与边缘态计算。构建系统与文档是后来整理的。

---

## 目录结构

```
.
├── waveguide_optimization/     C++ 求解器（libMesh + PETSc/SLEPc）
│   ├── main.cpp                主流程：读网格 → 给节点挂上 Bézier 参数 → 装配求解
│   │                           → 提取 S 参数 → 解伴随方程 → 组装 dK/dp → 输出
│   ├── 2D_waveguide_model_params.{h,cpp}   端口参数与解析波导模式 E0、β
│   ├── Bezier.{h,cpp}                      Bézier 求值/插值/最近点/控制点灵敏度
│   └── TriDiagonal.h                       三对角求解（Thomas 算法）
├── python/                     Python 侧（扁平模块，数据契约集中在一个文件里）
│   ├── make_mesh.py            入口：生成 gmsh 网格 + design.h5
│   ├── run_optimize.py         入口：跑优化循环
│   ├── design.py               几何与参数：DOF↔Bézier 控制点、端口表、gmsh 建模
│   │                           （网格尺寸场、几何自检、节点参数 t 都在这里）
│   ├── physics.py              材料/频率参数与平板波导色散（求 β）
│   ├── solver.py               求解器控制：写 infile.pot、调 mpiexec、校验、读结果
│   ├── optimizer.py            主优化器：目标/梯度封装 + 梯度下降
│   ├── postproc.py             后处理：读场、按 physical group 给网格上色、pyvista 可视化
│   ├── dataio.py               数据契约：design.h5 / result.h5
│   ├── paths.py                统一的路径/命令入口
│   └── tests/                  pytest：数据契约、设计变量映射、有限差分（可选）
├── data/                       运行输入（infile.pot、cross_waveguide.msh、design.h5）
├── docs/                       Typst 文档（理论、构建与运行、结果）
├── Makefile
├── config.mk.example
└── requirements.txt
```

## 依赖

| 依赖 | 用途 |
|---|---|
| libMesh（≥1.6，需 `--enable-complex`） | 有限元装配、网格、FEMContext、伴随求解接口 |
| PETSc / SLEPc | 复数线性求解（直接法 LU）、本征值问题 |
| Eigen | 矩阵与向量类型 |
| HDF5 | 前后端数据交换（`design.h5` / `result.h5`，经 `libmesh-config` 已在链接项里） |
| GetPot | 解析 `infile.pot`（libMesh 自带） |
| MPI | libMesh 的并行基础设施，实际以单进程运行 |
| Python: gmsh / pygmsh / pyvista / numpy / scipy / matplotlib / pybind11 | 网格生成、优化循环、绘图 |

## 编译与运行

```bash
# 0) 告诉 Makefile 你的 libMesh 装在哪（默认 $HOME/local）
cp config.mk.example config.mk && $EDITOR config.mk

# 1) 编译 C++ 求解器
make                 # 等价于 make METHOD=opt；也可以 make METHOD=dbg

# 2) 跑一次正问题（在 data/ 里读写，输出文件已被 .gitignore 忽略）
make run

# 3) 编译 Python 用的 Bézier 扩展模块，然后生成算例 / 跑优化循环
make python
python3 python/make_mesh.py          # 生成 data/cross_waveguide.msh + design.h5
cd data && python3 ../python/run_optimize.py
```

`make run` 的典型输出：

```
Bezier 参数: 356 个边界节点, 最大残差 7.58911e-05
Id of Ez: 0
end point mode:0
```

第二行是网格自检：`design.h5` 里的 `/design/node_t` 必须和这张网格的节点顺序对得上。
现在**网格完全由 Python 侧（gmsh）生成**：优化边界是真 Bézier 曲线，节点天生落在曲线上；
加密判据（端口 `h=0.0037`、介质内部 `λ/16`）也由 Python 的尺寸场给定。C++ 不再加密、不再
吸附节点，只负责求解。

产物是 `data/result.h5`（HDF5，h5py 直接读），含

| 数据集 | 含义 |
|---|---|
| `/scalars/S_params` | 四个端口的 S 参数（复数，存成 `(n,2)` + 属性 `complex=1`） |
| `/scalars/ex_opt_sens` / `in_opt_sens` | 目标函数对两条优化边界控制点坐标的灵敏度 |
| `/fields/Ez` | 电场自由度（`infile.pot` 里 `save_field = 1` 时才写） |
| `/meta` | 属性：`format_version` / `n_process` / `freq` / `epsilon_*` 等 |

各 rank 的局部量（S 参数的单元积分、按节点拼出的场）在 **C++ 侧用 MPI 归约**
后再由 rank 0 写出，所以文件里一律是全局量。几何输入在同一个目录的
`design.h5`（设计变量 + Bézier 控制点 + 端口表 + 完整物理/几何参数 + 网格尺寸参数
+ 每个节点的曲线参数 `node_t`）——它是自描述的，`design.Case.load()` 只靠这一个文件
就能还原出整个算例。

`data/infile.pot` 是求解器的全部运行参数（频率、介电常数、端口、Bézier 端点模式等），
格式是 GetPot。

## 方法

**物理模型。** 二维标量亥姆霍兹方程 $nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z mu_(x y) E^z = 0$，
波导芯 $epsilon = 11.56$（$n = 3.4$）、包层 $epsilon = 1$，
取 $epsilon_0 = mu_0 = 1$ 使所有量纲化成长度。端口用一阶透明边界
$diff E^z slash diff n - i beta E^z = -2 i beta E_0$（输入）/ $= 0$（输出）。

**几何参数化。** 交叉口的每条内凹边界用一条三次 Bézier 曲线描述。控制点并不自由：
为保证相邻曲线 $C^2$ 连续、端口处切向给定、节点处曲率连续，控制点由节点切矢
经一个三对角线性系统确定。因此真正的设计变量是节点坐标与切矢幅值。

**伴随法。** 离散系统 $bold(K)(bold(p)) bold(a) = bold(b)$，目标函数
$W = |S_(21)|^2$。先解一次伴随方程 $bold(K)^T overline(bold(a)) = (diff W)/(diff bold(x)) - i (diff W)/(diff bold(y))$，
再对每个设计变量组装 $(diff bold(K))/(diff p)$ 并做一次矩阵-向量乘：

$ (diff W)/(diff p) = - Re( overline(bold(a))^T (diff bold(K))/(diff p) bold(a) ) $

代价对比：有限差分需要 $N + 1$ 次线性求解，伴随法只要 2 次（正问题 + 伴随）。

**几何与网格的两道闸。** 优化器会奖励"贴着对称轴走"的形状，这类形状有两个
后果：壁面与自己的镜像相交（二维里表现为退化的针状单元，三维 slab 里就是一条
贯穿的漏光通道），以及网格质量下降。所以 `optimizer` 里加了三件事：

1. **端点角惩罚**。目标改成 $J = |S_(21)|^2 - lambda (|cos theta_"in"| + |sin theta_"ex"|)$，
   把两条边界的轴端切向轻轻推向"垂直于对称轴"（in 的端点在 x 轴上、ex 在 y 轴上，
   所以内积分别是 $cos$ / $sin$）。初始设计是 30° 斜着撞上轴的，所以不能硬固定成
   垂直，只能用罚项。实测端点角的**物理**灵敏度是 $10^(-4)$ 量级，罚项不会跟物理
   打架；它的导数解析可得（`design.end_axis_alignment_grad()`）。
2. **越轴截断 + 重拟合**。一旦检测到壁面越过对称轴，就在**首次越轴处**截断，轴端
   切向取原曲线在该处的切向，再用**同样的控制点个数**做 L-BFGS-B 拟合（按弧长
   对应、最小化与保留弧段的点对点距离）。实测拟合残差 $3.6 times 10^(-5)$，
   约等于局部单元尺寸的 0.3%。
3. **网格质量门**。`design.QUALITY_MIN = 0.05`，标定依据：初始直线设计的网格
   $q_min = 0.235$，而越轴设计会出现 $q approx 10^(-4)$ 的针状单元。质量门管的是
   **有限元精度**；几何切触造成的退化单元不归它管（加密救不了切触）。

**验证。** `python/tests/test_fd_sensitivity.py` 用随机扰动设计变量（扰动后**重划网格**，
因为网格是设计的函数）、拿中心差分的 $f(bold(p)+bold(delta)) - f(bold(p)-bold(delta))$ 与
伴随的 $bold(g) dot bold(delta)$ 比较（默认跳过，`WAVEGUIDE_FD_TEST=1` 才跑）。
实测 7 个方向/步长的符号全部一致，量级比（伴随/差分）落在 $0.95 ~ 1.27$，均值 $approx 1.12$、
标准差 $0.06$。残差主要来自"有限差分重划网格、伴随假设网格随节点连续变形"这一不一致，
以及 $lambda\/16$ 的离散误差。

## 结果

| 量 | 值 |
|---|---|
| 未优化交叉口的 \|S₂₁\|²（`make run` 复现） | 0.712 |
| 优化后 \|S₂₁\|² | ≈ 0.96 |

图与推导见 `docs/`。

## 说明

- `waveguide_optimization/TriDiagonal.cpp` 已在整理时删除：那份文件从未被 CodeBlocks
  编译过（`.cbp` 里没列它，`obj/` 里也没有对应的 `.o`），而且用到的 `TridiagonalMatrix`
  类型在任何地方都没有定义。实际在用的是 `TriDiagonal.h` 里的模板版本。
- 仓库历史里仍保留着 `waveguide_optimization.exe`（19 MB）和 `obj/*.o` 等编译产物，
  它们从某次提交起不再被跟踪。这是为了不重写历史。
- `data/` 之外的 `remote_*.msh/.mat` 等大体积中间数据不纳入版本控制。

## 文档

`docs/` 下是 Typst 文档，按内容分三部分：

| 文件 | 内容 |
|---|---|
| `docs/README.typ` | 文档索引、编译方法、符号约定 |
| `docs/overview.typ` | 项目总览：物理问题、数值方法、代码结构、数据流 |
| `docs/build-and-run.typ` | 依赖、编译、运行、参数说明、如何扩展 |
| `docs/theory/` | 有限元基础、二维麦克斯韦方程、波导模式、Bézier 参数化、伴随法、耦合微环 |
| `docs/results/` | 交叉波导优化的结果，以及微环阵列的能带与边缘态 |

```bash
# 需要先装 typst（https://typst.app/，或 cargo install typst-cli）
make docs          # 或 ./docs/build.sh，产物在 docs/build/
```

文档分三层：`overview.typ` 与 `build-and-run.typ` 讲项目与用法，`theory/` 讲物理与数值方法，
`results/` 放结果。每篇技术文档里的 `代码对应` 块会把它和具体函数挂钩。

## 参考文献

- H. Igarashi and K. Watanabe, *IEEE Trans. Magn.* **46**, 2739 (2010) —— 复数伴随变量法。
- M. Hafezi *et al.*, *Nature Photonics* **7**, 1001 (2013) —— 硅光子拓扑边缘态（作为后续工作的参照）。
- C. Johnson, *Numerical Solution of Partial Differential Equations by the Finite Element Method*,
  Cambridge University Press, 1987 —— 有限元基础。
