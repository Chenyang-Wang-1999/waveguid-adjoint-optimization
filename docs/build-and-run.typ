#import "template.typ": project-doc, note, fix, code-ref

#show: project-doc.with(
  title: "构建与运行",
  subtitle: "依赖、编译、运行、参数与扩展",
  date: "2026 年整理",
)

这份文档说明怎么把仓库跑起来，以及各个入口在代码里的位置。

== 依赖

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*依赖*], [*说明*]),
  [libMesh], [有限元装配、网格、`FEMContext`、伴随求解接口。安装时必须开 `--enable-complex`，否则 `Number` 是实数，整套复数伴随法无法编译],
  [PETSc / SLEPc], [复数线性求解（这里用直接法 LU）与本征值问题。由 libMesh 的 `libmesh-config` 一并带出],
  [Eigen], [`Eigen::MatrixXd` / `VectorXcd`，用于 Bézier 控制点与灵敏度矩阵],
  [HDF5], [C++ 侧读写 `design.h5` / `result.h5`（`waveguide_optimization/h5io.h`），连接 Python 与求解器；经 `libmesh-config` 已在链接项里],
  [GetPot], [解析 `data/infile.pot`，libMesh 自带],
  [MPI], [libMesh 的并行基础设施。本程序以单进程运行，但仍然链接 MPI],
  [Python], [`gmsh`/`pygmsh` 建网格，`numpy`/`scipy` 做优化与拟合，`pybind11` 把 `Bezier.cpp` 编成扩展模块],
)

Python 依赖见 `requirements.txt`。

== 配置 libMesh 的位置

仓库里*不存任何绝对路径*。`Makefile` 通过 libMesh 自带的 `libmesh-config` 推导全部
编译与链接参数：

```make
LIBMESH_DIR ?= $(HOME)/local
LIBMESH_CONFIG := $(LIBMESH_DIR)/bin/libmesh-config
CXX      := $(shell $(LIBMESH_CONFIG) --cxx)
CPPFLAGS := $(shell $(LIBMESH_CONFIG) --include)
LIBS     := $(shell $(LIBMESH_CONFIG) --libs) -lmatio
```

本机路径写进 `config.mk`（已被 `.gitignore` 忽略）：

```bash
cp config.mk.example config.mk
$EDITOR config.mk        # 改 LIBMESH_DIR / METHOD / PYTHON
```

`METHOD` 选择 libMesh 的编译模式（`opt` / `dbg` / `devel`，取决于 libMesh 编译时的
configure 选项），插值到 `libmesh-config --libs` 给出的 `-lmesh_$(METHOD)` 上。

#note[
  这里没有用 CMake，因为 libMesh 不安装 CMake 的 package config
  （没有 `libmesh-config.cmake`），`find_package(libMesh)` 走不通；而 `libmesh-config`
  已经把 PETSc、SLEPc、Eigen、HDF5 的 include/lib 一次性带出来了。
]

== 编译

```bash
make                # 编译 C++ 求解器，产物 waveguide_optimization/waveguide_optimization
make METHOD=dbg     # 用 libMesh 的 dbg 模式
make python         # 编译 python/Bezier.cpython-*.so（pybind11 扩展模块）
make clean
```

对象文件与依赖关系由 `-MMD -MP` 自动生成，改头文件会自动重编。

=== Python 扩展模块的一个坑

`Bezier.cpp` 在定义了 `EXPORT_PYLIB` 时会被 pybind11 编成 Python 模块，
`python/design.py` 用 `import Bezier` 直接调用，这样 Python 侧生成网格时用的
Bézier 求值/灵敏度代码与 C++ 侧*完全同一份*，不会出现两边公式不一致。

坑在于：本机 `/usr/include/python3.8` 落在编译器的默认头文件搜索路径里，如果只加
pybind11 的 include 而不加解释器自己的 include，编译器会拿 *3.8* 的 `Python.h` 去编，
虽然输出文件名后缀是 `cpython-39`，`import` 时仍然报版本不匹配。所以 `Makefile` 里
显式加了 `$(PYTHON)-config --includes`：

```make
PYTHON_INC ?= $(shell $(PYTHON)-config --includes 2>/dev/null)
```

== 运行

=== 只跑正问题

```bash
make run                       # 等价于 cd data && mpiexec -n 1 ../waveguide_optimization/waveguide_optimization
make run ARGS="-pc_type lu"    # 额外参数直接透传给 PETSc
```

求解器从*当前目录*读 `infile.pot`，所以 `make run` 会先 `cd data`。
典型输出：

```text
Bezier 参数: 356 个边界节点, 最大残差 7.58911e-05
Id of Ez: 0
end point mode:0
```

第一行是网格自检：`design.h5` 的 `/design/node_t` 必须和这张网格的节点顺序一致，
C++ 会把每个节点按其参数 $t$ 反算出的曲线点比对一遍（残差 $< 10^(-3)$ 才放行）。

跑完在 `data/` 下得到：

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*内容*]),
  [`result.h5`], [`/scalars/S_params`（四个端口的 S 参数）、`/scalars/ex_opt_sens` / `in_opt_sens`（灵敏度）、可选的 `/fields/Ez`、`/meta`],
  [`out_cross_waveguide.msh`], [求解器写回的网格（节点已经在曲线上，不再有任何吸附/变形），可以用 ParaView 看],
)

#code-ref[
  整个流程在 `waveguide_optimization/main.cpp` 的 `main()` 里，顺序是：
  读 `infile.pot`（`WaveguideParams curr_params`）→ `mesh.read()`（关闭节点重排）→
  `tag_boundary_param()` 从 `design.h5` 的 `/design/node_t` 给节点挂上曲线参数 →
  `eq_systems.init()` / `solve()` → `integrate_along_side(..., mode_decompose)` →
  `assemble_QOI_rhs()` + `adjoint_solve()` → 双重循环 `assemble_sensitivity_K()` →
  归约后写 `result.h5`。
]

网格（几何、尺寸场、节点参数 $t$）全部由 `python/design.py` 的 `build_mesh()` 产出。

=== 跑优化循环

```bash
make python
cd data && python3 ../python/run_optimize.py
```

`optimizer.run()` 做的是：

+ `design.make_case()` 按当前设计变量建网格、写 `design.h5`（第一次用初始直线设计，
  之后若 `DOF_list.csv` 已存在则从最后一轮续跑）；
+ `solver.run()` 跑一次正问题 + 伴随问题，`dataio.read_result()` 读回 `S_params` 与两组灵敏度；
+ `design.sensitivity_to_dofs()` 把灵敏度拼回设计变量顺序；
+ 按 `DOF += alpha * sens` 更新控制点，把每步的 S 参数、梯度、DOF 追加到三个 csv；
+ 反复 10 轮。

注意**网格是设计变量的函数**：每轮都会重划（`build_mesh()` 约 1 s），所以"设计变量 →
网格 → 物理量"是一条确定的链，不存在"只动边界不重划"的路径。

=== 校核灵敏度

```bash
cd data && WAVEGUIDE_FD_TEST=1 python3 -m pytest ../python/tests/test_fd_sensitivity.py -s
```

`test_fd_sensitivity.py` 沿随机方向扰动设计变量（扰动后重划网格）再跑正问题，比较

$ bold(g) dot bold(delta) quad "（伴随）" <=> quad f(bold(p) + bold(delta)) - f(bold(p)) quad "（有限差分）" , quad f = |S_(21)|^2 $

实测 7 个方向/步长的符号全部一致，量级比（伴随/差分）落在 $0.95 ~ 1.27$，
均值 $approx 1.12$、标准差 $0.06$。这个残差**不是**符号或数量级问题，更可能来自
"有限差分重划网格、伴随假设网格随节点连续变形"这一不一致，以及 $lambda\/16$ 的离散误差；
要彻底对上需要另做一组固定网格的实验。

== 参数：`data/infile.pot`

求解器的全部运行参数都在这个 GetPot 文件里，由 `WaveguideParams` 构造函数读取：

#table(
  columns: (auto, auto, 1fr),
  align: (left, left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*键*], [*典型值*], [*含义*]),
  [`mesh_file`], [`cross_waveguide.msh`], [输入网格（Gmsh 格式）],
  [`design_file`], [`design.h5`], [Bézier 控制点、切矢、端口参数、`node_t`],
  [`result_file`], [`result.h5`], [输出文件],
  [`out_mesh_file`], [`out_cross_waveguide.msh`], [求解器写回的网格],
  [`epsilon_wg` / `epsilon_amb`], [`11.56` / `1.0`], [波导芯 / 包层介电常数（$n = 3.4$ 的平方）],
  [`mu_wg` / `mu_amb`], [`1.0` / `1.0`], [磁导率],
  [`freq`], [`1.550000`], [归一化频率（对应 1.55 µm）],
  [`start_point_mode` / `end_point_mode`], [`0` / `0`], [Bézier 端部条件：0 = 给定切矢，1 = 二阶导为零],
  [`save_field`], [`0`], [是否把 `Ez` 自由度写进 `result.h5` 的 `/fields`],
  [`id_wg` / `id_ex_opt` / `id_in_opt`], [`1` / `7` / `8`], [physical id；给了就不必按名字查（名字在格式转换时容易丢）],
)

端口相关的数组（`port_id`、`port_coord`、`port_x_axis`、`incident_on`、
`port_other_params`、`port_is_analytic`）不在 `infile.pot` 里，而是由
`python/design.py` 在生成网格时一并写进 `design.h5` 的 `/ports` 组。
`infile.pot` 里的 `*_varname` 就用来指定这些表在组里的名字。

== 目录与数据流

```text
python/design.py : build_mesh()        data/cross_waveguide.msh
  ├─ 分段三次 Bézier + 尺寸场（gmsh）──►   （节点天生在曲线上）
  ├─ 解色散方程求 β ──────────────────► data/infile.pot
  ├─ 算每个边界节点的参数 t ──────────► data/design.h5
  └─ 写 Bézier 控制点 / 端口表 / 尺寸参数
                                                          │
                                                          ▼
                            waveguide_optimization/waveguide_optimization  (C++)
                                     │  读 .msh / design.h5 / infile.pot
                                     ▼
                            data/result.h5
                                     │  S 参数 + 灵敏度
                                     ▼
python/optimizer.py  ──►  更新控制点，回到第一步（重划网格）
```

== 想改点什么

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*目标*], [*改哪里*]),
  [换目标端口], [`main.cpp` 里的 `unsigned int selected_port = 1;` 与 `QOI_port`],
  [换目标函数（例如串扰）], [`assemble_QOI_rhs()` 里的 $W$；同时改 `python/optimizer.py` 里取 `S_params` 的第几个分量],
  [改变设计变量], [`Bezier_tangent_vector()` 的 `start_point_mode` / `end_point_mode`，以及 `main.cpp` 里 `ex_points` / `in_points` 的计算],
  [加密网格], [`python/design.py` 的 `PORT_SIZE` / `WL_PER_SIZE` / `GRADE_LEN`（尺寸场），或 `infile.pot` 里的频率（内部尺寸由它定）],
  [换优化器], [`python/optimizer.py` 里 `gradient_descent()` 换成 `scipy.optimize.minimize`],
  [改端口模式], [`2D_waveguide_model_params.cpp` 的 `analytic_E_field()`],
)
