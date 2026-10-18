#import "template.typ": project-doc, fig, fig-row, note, fix, code-ref

#show: project-doc.with(
  title: "项目总览：交叉波导的伴随法形状优化",
  subtitle: "Adjoint shape optimization of a photonic cross-waveguide with libMesh",
  date: "2022 年 2 月",
)

== 要解决的问题

两条介质波导正交交叉时，交叉口是一个强散射体：即使波导本身单模、低损耗，直通端也会因为
模式失配和四角处的场集中损失几个百分点，同时产生串扰。教科书做法是“切角”或加一段锥形过渡，
但形状参数只能手工试。

本项目把交叉口边界写成三次 Bézier 曲线，把*控制点的坐标*当作设计变量，用*伴随法*
（adjoint method）求出目标函数（直通端透射率 $|S_(21)|^2$）对全部设计变量的梯度，
再由 Python 侧的优化器迭代更新控制点，自动把交叉口“捏”成高透射的形状。

和像素级拓扑优化相比，这种“边界参数化 + 伴随梯度”的做法有两个好处：

- 设计变量只有十几个到几十个（每条边界 12 个控制点），而不是每个网格单元一个自由度，
  因此不需要密度滤波、不需要处理灰度单元，优化结果是光滑可制造的曲线；
- 网格随控制点连续变形（mesh morphing），不需要在每次迭代里重新做拓扑改变，
  有限元矩阵的稀疏结构保持不变。

== 结果

#table(
  columns: (auto, auto, 1fr),
  align: (left, left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*量*], [*值*], [*来源*]),
  [未优化交叉口的 $|S_(21)|^2$], [0.712], [本仓库 `make run` 实测（见 `build-and-run.typ`）],
  [优化后 $|S_(21)|^2$], [$approx 0.96$], [2022 年 2 月的记录，见 `results/01-cross-waveguide.typ`],
)

2022 年 3–4 月，同一套求解器被用来给微环阵列算场与 S 参数，进而提取耦合系数、拟合
Haldane 有效模型、计算能带与边缘态。那部分内容在 `results/02-ring-array.typ`。

== 物理模型

只考虑二维问题（结构在传播方向上不变），取横向电场分量 $E^z$ 满足标量亥姆霍兹方程

$ nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z mu_(x y) E^z = 0 $

其中 $epsilon_z$、$mu_(x y)$ 在波导芯和包层中分别取不同常数。完整的推导（从三维麦克斯韦
方程出发，在 $partial_z = i beta_z$ 的分层假设下分离出 TE/TM 两种模式）见
`theory/02-maxwell-2d.typ`。

=== 单位制

全文取 $epsilon_0 = mu_0 = 1$，于是所有物理量都量纲化为长度，$c = 1$。实际取值：

#table(
  columns: (auto, auto, 1fr),
  align: (left, left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*量*], [*取值*], [*说明*]),
  [波导芯介电常数], [$epsilon_"wg" = 11.56$], [对应 $n = 3.4$，硅],
  [包层介电常数], [$epsilon_"amb" = 1.0$], [空气],
  [磁导率], [$mu_"wg" = mu_"amb" = 1.0$], [非磁性],
  [归一化频率], [$f = 1.55$], [对应波长 1.55 µm 的归一化写法],
  [波导宽度], [$w_"wg" = 0.08$], [约 $0.11 lambda$，保证单模],
  [固定直波导长度], [$L_"fixed" = 0.3$], [交叉口两侧的直波导段],
  [中心优化区长度], [$L_"opt" = 0.2 ~ 0.5$], [Bézier 边界覆盖的区域],
)

=== 几何参数化

交叉口的四条内凹边界各用一条三次 Bézier 曲线描述，$B(t)$ 由两端点 $P_0, P_1$ 与两个控制点
$C_0^+, C_1^-$ 决定。为了让相邻曲线 $C^2$ 连续、并在端口处保持给定切向、在连接节点处曲率
连续，控制点并不是自由的：它们由节点上的切矢 $V_i$ 通过一个三对角线性系统确定，
所以真正的设计变量是节点坐标 $P_i$ 与切矢幅值，实现上统一编码成一维向量 `DOF`。
细节见 `theory/04-bezier.typ`。

=== 边界条件

- *端口*：一阶透明边界（与解析模式匹配）。输入端口

  $ (diff E^z)/(diff n) - i beta E^z = -2 i beta E_0 $

  输出端口 $(diff E^z)/(diff n) - i beta E^z = 0$。$E_0$ 是波导解析基模，$beta$ 由色散方程
  数值求解（`theory/03-slab-waveguide.typ`）。
- *侧壁*：$E^z = 0$（金属壁，仅用于平板波导的验证算例）。

== 数值方法

=== 有限元

用一个二阶拉格朗日单元（libMesh 变量 `Ez`，类型 `SECOND`）离散 $E^z$，五阶 Gauss 积分。
线性系统是复对称（非 Hermite）的，用 PETSc 直接法求解（`-pc_type lu`）。

网格来自 Gmsh（`.msh`），读入后由 `EMMeshFlagging` 按“每波长单元数足够、且按折射率加权”
的准则自适应加密，直到没有单元需要再分（实测加密 3 轮）。加密完成后再把优化边界上的节点
“吸附”到 Bézier 曲线上，让网格跟着控制点变形。

#code-ref[
  `waveguide_optimization/EM_mesh_flagging.cpp` 的 `flag_elements()` 给出加密判据，
  `main.cpp` 的 `main()` 里循环调用 `MeshRefinement::refine_elements()`；
  吸附在 `mesh_snap_to_Bezier()` 里做，用 `Bezier::Nearest_t()` 把每个边界节点投影回曲线。
]

=== S 参数提取

S 参数用模式分解得到：在端口边界上做重叠积分

$ S_(i) = integral_Gamma 1/4 (bold(D)_"mode"^* dot bold(E) + bold(H)_"mode"^* dot bold(B)) dif Gamma $

这个归一化就是复坡印廷功率的 $1/4$ 系数形式，因此积分出来的直接是模式幅度即 S 参数。
$H_x$、$H_y$ 由 $E^z$ 的切向梯度给出：$H_x = -i/(omega mu) partial_y E^z$、
$H_y = i/(omega mu) partial_x E^z$。

#code-ref[
  `main.cpp` 的 `mode_decompose()` 就是上面这个被积函数，
  由 `integrate_along_side()` 在选定的端口边界上积分。
  解析的 $bold(D)_"mode"$、$bold(B)_"mode"$ 来自
  `2D_waveguide_model_params.cpp` 的 `analytic_E_field()`。
]

=== 目标函数与伴随梯度

取直通端（端口 1）的透射率 $W = |S_(21)|^2$。对 $bold(K)(bold(p)) bold(a)(bold(p)) = bold(b)$，
先解一次伴随方程

$ bold(K)^T overline(bold(a)) = (diff W)/(diff bold(x)) - i (diff W)/(diff bold(y)) $

再对每个设计变量组装 $(diff bold(K))/(diff p)$，梯度为

$ (diff W)/(diff p) = - Re( overline(bold(a))^T (diff bold(K))/(diff p) bold(a) ) $

关键点是：*伴随方程只解一次*，与设计变量个数无关。这正是能用几十个设计变量做梯度优化的
原因——用有限差分的话，每次迭代至少要几十次全场求解。完整推导（含右端项
$bold(g) = 2 S overline((diff S)/(diff A))$ 与代码的逐项对应）见 `theory/05-adjoint-sensitivity.typ`。

#code-ref[
  伴随右端在 `main.cpp` 的 `assemble_QOI_rhs()` 里装配；
  `LinearSolver::adjoint_solve()` 解转置系统；
  逐设计变量的 $(diff bold(K))/(diff p)$ 由 `assemble_sensitivity_K()` 组装，
  其中边界对控制点的导数来自 `Bezier.cpp` 的 `Bezier_V_sensitivity()`。
]

== 代码结构

项目分成“Python 侧”和“C++ 侧”两部分，通过 `.mat` 文件与进程调用交换数据。

=== C++ 侧：`waveguide_optimization/`

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*职责*]),
  [`main.cpp`], [`main()`：读网格 → 加密 → 吸附到 Bézier → 装配求解 → 提取 S 参数 → 解伴随方程 → 组装 $(diff K)/(diff p)$ → 输出。另外含 `assemble_waveguide()`、`mode_decompose()`、`integrate_along_side()`、`assemble_QOI_rhs()`、`assemble_sensitivity_K()`、`mesh_snap_to_Bezier()`],
  [`2D_waveguide_model_params.{h,cpp}`], [`WaveguideParams`：读 `infile.pot` 与 `.mat`；`analytic_E_field()` 给出端口解析模式的 $E_0, beta$ 及其 $H$、$B$ 分量],
  [`EM_mesh_flagging.{h,cpp}`], [`EMMeshFlagging::flag_elements()`：按频率与折射率给出加密判据],
  [`Bezier.{h,cpp}`], [`Bezier_tangent_vector()`（三对角系统求切矢）、`Bezier_interp()`、`Nearest_t()`、`Bezier_V_sensitivity()`；定义 `EXPORT_PYLIB` 时可编成 Python 模块],
  [`TriDiagonal.h`], [`Solve_tridiagonal()`：Thomas 算法的模板实现],
)

=== Python 侧：`python/`

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*职责*]),
  [`make_mesh.py`], [入口：生成 gmsh 网格与 `design.h5`],
  [`run_optimize.py`], [入口：跑优化循环],
  [`design.py`], [几何与参数：DOF ↔ Bézier 控制点、端口表、gmsh 建模（唯一产生网格的地方）],
  [`physics.py`], [材料/频率参数与平板波导色散（求 $beta$，含频率扫描时用上一个解当初值）],
  [`solver.py`], [求解器控制：写 `infile.pot`、调 `mpiexec`、检查退出码、读回结果],
  [`optimizer.py`], [主优化器：目标/梯度封装 + 梯度下降],
  [`postproc.py`], [后处理：读场、`pyvista` 可视化],
  [`dataio.py`], [数据契约：`design.h5` / `result.h5` 的唯一读写入口],
  [`paths.py`], [统一的路径与外部命令入口，所有脚本都在 `data/` 下工作],
  [`tests/`], [pytest：数据契约、设计变量映射、有限差分（可选）],
)

=== 数据流

```text
python/make_mesh.py                    data/cross_waveguide.msh
  ├─ gmsh 建几何 ─────────────────────►  data/design.h5   (Bézier 控制点 + 端口表)
  └─ 解色散方程求 β                                        │
                                                           ▼
                             waveguide_optimization/waveguide_optimization  (C++)
                                      │  solver.py 写 infile.pot 并调用
                                      ▼
                             data/result.h5   (/scalars, /fields, /meta)
                                      │  S 参数 + 灵敏度
                                      ▼
python/run_optimize.py  ──►  更新控制点，回到第一步
```

构建与运行方式见 `build-and-run.typ`。

== 依赖

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*依赖*], [*用途*]),
  [libMesh], [有限元装配、网格管理、`FEMContext`、伴随求解接口（需 `--enable-complex`）],
  [PETSc / SLEPc], [复数线性求解（直接法 LU）、本征值问题],
  [Eigen], [矩阵与向量类型],
  [HDF5], [前后端数据交换（`design.h5` / `result.h5`）],
  [GetPot], [解析 `infile.pot` 运行参数],
  [MPI], [libMesh 的并行基础设施，实际以单进程运行],
  [Python: gmsh / numpy / scipy / matplotlib / pybind11], [网格生成、优化循环、拟合与绘图],
)

== 文档地图

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文档*], [*内容*]),
  [`build-and-run.typ`], [依赖、配置、编译、运行、参数、想改点什么时改哪里],
  [`theory/01-fem-basics.typ`], [有限元方法基础：变分形式、Galerkin、误差估计、Poisson 方程],
  [`theory/02-maxwell-2d.typ`], [二维麦克斯韦方程、TE/TM 分解、三类边界条件、金属波导端口边界、平板波导解析解],
  [`theory/03-slab-waveguide.typ`], [平板波导传播常数：色散方程与牛顿迭代的数值稳定性],
  [`theory/04-bezier.typ`], [Bézier 曲线参数化：$C^2$ 拼接、插值、控制点灵敏度],
  [`theory/05-adjoint-sensitivity.typ`], [伴随法灵敏度分析：推导、右端项、代价、验证],
  [`theory/06-coupled-resonators.typ`], [耦合微环：input–output formalism 与绝热消除],
  [`theory/07-bands-and-edge-states.typ`], [微环阵列的能带与边缘态模型],
  [`notes/libmesh-notes.typ`], [libMesh 使用笔记],
  [`results/01-cross-waveguide.typ`], [交叉波导优化的结果],
  [`results/02-ring-array.typ`], [微环阵列的数值结果：耦合系数、能带、边缘态、赝自旋耦合],
)

== 说明

#note[
  代码写在 2022 年 2 月（见 `git log`，首次提交 `dca51a9` 在 2022-02-21），
  2022 年 3–4 月又用它做了微环阵列的能带与边缘态计算。构建系统与文档是后来整理的：
  原来用 Code::Blocks 的 `.cbp` 工程（依赖绝对路径和一个已经不存在的自定义编译器），
  现在改成 `Makefile` + libMesh 自带的 `libmesh-config`。
  整理过程没有改写 git 历史，原始提交保持原样。
]
