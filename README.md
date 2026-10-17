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
│   ├── main.cpp                主流程：读网格 → 加密 → 吸附到 Bézier → 装配求解
│   │                           → 提取 S 参数 → 解伴随方程 → 组装 dK/dp → 输出
│   ├── 2D_waveguide_model_params.{h,cpp}   端口参数与解析波导模式 E0、β
│   ├── EM_mesh_flagging.{h,cpp}            按频率与折射率自适应加密网格
│   ├── Bezier.{h,cpp}                      Bézier 求值/插值/最近点/控制点灵敏度
│   └── TriDiagonal.h                       三对角求解（Thomas 算法）
├── python/                     Python 侧
│   ├── original_mesh.py        Gmsh 生成交叉波导网格，写出 .msh/.mat/infile.pot
│   ├── optimize.py             优化主循环：DOF ↔ 控制点互转、调求解器、梯度下降
│   ├── test_sensitivity.py     有限差分校核伴随灵敏度
│   ├── propagation_constant0.py 平板波导色散方程的牛顿迭代
│   ├── interp_DOF_to_points.py  把有限元自由度插值到指定点
│   └── paths.py                统一的路径/命令入口
├── data/                       运行输入（infile.pot、cross_waveguide.msh/.mat）
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
| matio | 读写 MATLAB `.mat`（C++ 侧） |
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

# 3) 编译 Python 用的 Bézier 扩展模块，然后跑优化循环
make python
cd data && python3 ../python/optimize.py
```

`make run` 的典型输出：

```
Curr refinement: 1 / 2 / 3
Id of Ez: 0
end point mode:0
```

产物是 `data/cross_waveguide_result.mat`，里面含

| 键 | 含义 |
|---|---|
| `S_params_0` | 四个端口的 S 参数（复数） |
| `ex_opt_sens` / `in_opt_sens` | 目标函数对两条优化边界控制点坐标的灵敏度 |
| `Ez_0` | 电场自由度（`infile.pot` 里 `save_field = 1` 时才写） |
| `N_process` | MPI 进程数 |

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

**验证。** `python/test_sensitivity.py` 用随机扰动设计变量、比较有限差分与伴随梯度。

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
make docs          # 或 ./docs/build.sh，产物在 docs/build/
```

## 参考文献

- H. Igarashi and K. Watanabe, *IEEE Trans. Magn.* **46**, 2739 (2010) —— 复数伴随变量法。
- M. Hafezi *et al.*, *Nature Photonics* **7**, 1001 (2013) —— 硅光子拓扑边缘态（作为后续工作的参照）。
- C. Johnson, *Numerical Solution of Partial Differential Equations by the Finite Element Method*,
  Cambridge University Press, 1987 —— 有限元基础。
