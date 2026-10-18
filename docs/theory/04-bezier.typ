#import "../template.typ": project-doc, fig, note, fix, code-ref

#show: project-doc.with(
  title: "Bézier 曲线与网格边界的参数化",
  subtitle: "C² 连续拼接、插值与控制点灵敏度",
  date: "2022-02",
)

优化问题的设计变量就是交叉口边界上的 Bézier 曲线的控制点。这一节把曲线的定义、拼接条件、
插值公式和对控制点的导数（伴随法要用）集中记一遍。实现见 `Bezier.cpp` / `Bezier.h`。

== 单段三次 Bézier 曲线

一段三次 Bézier 曲线由两个端点 $P_0, P_1$ 和两个控制点 $C_0^+, C_1^-$ 决定：

$ B(t) = P_0 (1-t)^3 + 3 C_0^+ t (1-t)^2 + 3 C_1^- t^2 (1-t) + P_1 t^3 $

一阶、二阶导数：

$ B'(t) = 3 C_0^+ (t-1)^2 - 3 P_0 (t-1)^2 + 6 C_0^+ t (t-1) - 6 C_1^- t (t-1) + 3 P_1 t^2 - 3 C_1^- t^2 $

$ B''(t) = 18 C_0^+ t - 12 C_0^+ - 18 C_1^- t + 6 C_1^- - 6 P_0 t + 6 P_0 + 6 P_1 t $

=== 端点处的值

#table(
  columns: (auto, auto),
  align: (center, center),
  stroke: none,
  inset: (x: 1em, y: 0.3em),
  table.header([*$t = 0$*], [*$t = 1$*]),
  [$B(0) = P_0$], [$B(1) = P_1$],
  [$B'(0) = 3(C_0^+ - P_0)$], [$B'(1) = 3(P_1 - C_1^-)$],
  [$B''(0) = -12 C_0^+ + 6 C_1^- + 6 P_0$], [$B''(1) = 6 C_0^+ - 12 C_1^- + 6 P_1$],
)

== 多段拼接：$C^2$ 连续的控制点

把第 $i$ 段记作 $B_i$，它连接 $P_i = B_i(0)$ 与 $P_(i+1) = B_i(1)$，控制点是
$C_i^+, C_i^-$。逐段写出端点条件：

#table(
  columns: (auto, auto),
  align: (center, center),
  stroke: none,
  inset: (x: 1em, y: 0.3em),
  table.header([*$t = 0$*], [*$t = 1$*]),
  [$B_i(0) = P_i$], [$B_i(1) = P_(i+1)$],
  [$B_i'(0) = 3(C_i^+ - P_i)$], [$B_i'(1) = 3(P_(i+1) - C_(i+1)^-)$],
  [$B_i''(0) = -12 C_i^+ + 6 C_(i+1)^- + 6 P_i$], [$B_i''(1) = 6 C_i^+ - 12 C_(i+1)^- + 6 P_(i+1)$],
)

=== 内部节点：导数连续

要求相邻两段在节点处一阶导数相等，把这个公共切矢记作 $V_i$：

$ 3(C_i^+ - P_i) = 3(P_i - C_i^-) =: V_i $

于是控制点完全由节点坐标和切矢决定：

$ C_i^+ = P_i + V_i/3, quad C_i^- = P_i - V_i/3 $

代回二阶导数的表达式，得到用 $P$ 和 $V$ 表示的形式：

$ B_i''(0) = 6 P_(i+1) - 6 P_i - 4 V_i - 2 V_(i+1) $

$ B_i''(1) = 6 P_i - 6 P_(i+1) + 2 V_i + 4 V_(i+1) $

=== 内部节点：二阶导数连续

再要求 $B_(i-1)''(1) = B_i''(0)$：

$ 3 P_(i+1) - 3 P_i - 2 V_i - V_(i+1) = 3 P_(i-1) - 3 P_i + V_(i-1) + 2 V_i $

整理得

$ 3 P_(i+1) - 3 P_(i-1) = V_(i+1) + 4 V_i + V_(i-1) $

这就是一个*三对角*线性系统：对每个内部节点 $i$，系数是 $(1, 4, 1)$，右端是
$3(P_(i+1) - P_(i-1))$。所以给定节点坐标 $P_i$ 后，切矢 $V_i$ 可以 $O(N)$ 地解出来
（`TriDiagonal.h` 里的 `Solve_tridiagonal`，Thomas 算法）。*这就是为什么设计变量取节点
坐标而不是控制点坐标*：控制点由节点唯一确定，且求解是线性的。

#code-ref[`waveguide_optimization/Bezier.cpp` 的 `Bezier_tangent_vector()` 就是这个三对角系统：
  循环 `for(unsigned int i=1; i<N_points-1; i++)` 里 `rhs[i] = 3*P.row(i+1) - 3*P.row(i-1)`、
  `a[i-1] = 1.0; b[i] = 4.0; c[i] = 1.0`，然后把 `a, b, c, rhs` 交给
  `TriDiagonal.h` 的 `Solve_tridiagonal()`（模板函数，`beta`/`rho` 两个数组做前向消元、
  再回代），解出的 `sol[i]` 就是 $V_i$，最后拼成 `Eigen::MatrixXd V`（每行一个节点的切矢）返回。
  注意模板参数 `T` 在 C++ 侧是 `Eigen::VectorXd`（每个节点是一个二维矢量），所以
  $a, b, c$ 用 `double`、右端和解用向量。]

=== 边界条件

首尾两端各差一个方程，由两种模式之一补上：

+ *给定切矢*（`start_point_mode = 0` / `end_point_mode = 0`）：
  $ V_0 = V_"start" $（右端同理 $ V_N = V_"end" $），即两端端口处的切向是外部给定的，
  在代码里就是 `Bezier_tangent_vector()` 的形参 `vec_start` / `vec_end`
  （最后一行取 `b[N_points-1] = 1`、`a[N_points-2] = 0`，右端是 `vec_end`）。
+ *二阶导数为零*（`start_point_mode = 1` / `end_point_mode = 1`）：

  $ 3 P_1 - 3 P_0 = 2 V_0 + V_1 quad "(左端)" $

  $ 3 P_N - 3 P_(N-1) = V_(N-1) + 2 V_N quad "(右端)" $

  代码里对应 `b[0] = 2; c[0] = 1; rhs[0] = 3*(P.row(1) - P.row(0))` 与
  `a[N_points-2] = 1; b[N_points-1] = 2; rhs[N_points-1] = 3*(P.row(N-1) - P.row(N-2))`。

#fix[
  原稿把"给定切矢"这一条写作 $V_i = V_0$（或 $V_N$）：未知量 $V_i$ 和给定的数据用了同一个
  符号 $V_0$，读不出哪边是已知的。按代码的记法应是左端 $V_0 = V_"start"$、右端
  $V_N = V_"end"$（数据即 `vec_start` / `vec_end`）。这一条没有改数学，只把已知量换了个名字。
]

这两种模式在 `Bezier_tangent_vector()` 里分别对应 `start_point_mode`/`end_point_mode`
等于 0 和 1。代码里还有一个 `autoscale` 开关：打开时把给定的切矢按
$|P_1 - P_0| slash |V_"start"|$ 缩放，使切矢的长度与节点间距匹配，避免控制点跑得离曲线太远
（`rhs[0] = vec_start*tan_vec_norm/vec_norm`，右端同理用
$(P_(N-1) - P_(N-2))$ 归一化）。本项目 `data/infile.pot` 里取
`start_point_mode = 0`、`end_point_mode = 0`，即两端都用给定切矢。

== 插值与导数

在反问题里（把网格节点吸附到曲线上）需要频繁计算 $B(t)$ 及其导数。实现里统一用
$C_i^+ = P_i + V_i slash 3$、$C_(i+1)^- = P_(i+1) - V_(i+1) slash 3$ 代换，所以下面
三式分别用 $P, C$（$B$ 本身）、$P, V$（$B'$）和 $C, V$（$B''$）写出：

$ B(t) = P_0 (1-t)^3 + 3 C_0^+ t (1-t)^2 + 3 C_1^- t^2 (1-t) + P_1 t^3 $

$ B'(t) = V_0 (t-1)^2 + (6 P_0 - 6 P_1 + 2 V_0 + 2 V_1) t (t-1) + V_1 t^2 $

$ B''(t) = 2(6 C_0^+ - 6 C_1^- + V_0 + V_1) t - 6 C_0^+ + 6 C_1^- - 2 V_0 $

`Bezier_interp(P, V, t, diff_order)` 按 `diff_order` 返回其中一档。
`Nearest_t()` 用一维搜索求给定点对应的参数 $t$（默认绝对容差 $10^(-2)$，
`Bezier.h` 里的 `double AbsTol=1e-2`）。它现在由 Python 侧的
`design.boundary_param_t()` 调用，给每个边界节点算一次 $t$ 写进 `design.h5`
（C++ 侧不再投影节点——网格是直接在真曲线上生成的）。

顺带一个有用的等价关系：上式就是控制多边形
$[P_i, thin P_i + V_i\/3, thin P_(i+1) - V_(i+1)\/3, thin P_(i+1)]$ 的三次
**Bézier** 段。因此把每段单独建成一条 gmsh 三次 Bézier，得到的几何与
`Bezier_interp()` 逐位相同——这正是"网格节点天生落在曲线上"的依据
（见 `design.bezier_polygons()`）。

#code-ref[`Bezier_interp()` 逐字实现了上面三式。参数 $t$ 先被拆成"第几段 + 段内参数"：
  `id = (int)t(i)`、`loc_t = t(i) - id`，段内用 `Cp = P.row(id) + V.row(id)/3`、
  `Cn = P.row(id+1) - V.row(id+1)/3` 得到两个控制点，然后
  `B.row(i) = P.row(id)*loc_inv_t**3 + 3*Cp*loc_inv_t**2*loc_t + 3*Cn*loc_inv_t*loc_t**2
  + P.row(id+1)*loc_t**3`（`diff_order = 0`）；
  `diff_order = 1` 用 $B'(t)$ 的 $P, V$ 形式，`diff_order = 2` 用 $B''(t)$ 的 $C, V$ 形式。
  整数节点上（`loc_t` 接近 0）它直接用端点表里的值：`diff_order = 2` 时若 `id` 是最后一个点，
  返回 $6P_(id-1) - 6P_(id) + 2V_(id-1) + 4V_(id)$，否则返回
  $6P_(id+1) - 6P_(id) - 4V_(id) - 2V_(id+1)$——正是上一节用 $P, V$ 表示的那两条二阶导数。

  `Nearest_t()` 分两步：先在 `linspace(0, N_calc-1, 100)` 上采样取最近的折线段
  （`distance(curr_row) = (P.row(curr_row) - pt.transpose()).norm()`），再以它为初值做牛顿迭代
  （`T_test(0) += - First_diff/Second_diff`，最多 1000 次，收敛判据是
  $|B'(t) dot (B(t) - p)| < "AbsTol"$）。`python/design.py` 的 `boundary_param_t()` 对
  `ex_opt` / `in_opt` 边界上的每个节点先按象限取绝对值，调用
  `Bezier.Nearest_t(P, V, pt, 1e-8)` 得到 $t$ 写进 `design.h5` 的 `/design/node_t`；
  C++ 侧的 `tag_boundary_param()` 按节点 id 取回 $t$，并用 `Bezier_interp(..., 0)`
  反算曲线点与节点坐标比对（残差 $< 10^(-3)$ 才放行），既是自检也是"节点顺序没错位"的证明。]

== 对控制点的灵敏度

这是伴随法梯度的最后一环。设 $alpha$ 是某个设计变量（节点坐标 $P$ 或切矢 $V$ 的分量），
对 $B(t)$ 求导：

$ (diff B(t))/(diff alpha) = (diff P_i)/(diff alpha) (1-t)^3 + 3 (diff C_i^+)/(diff alpha) t (1-t)^2 + 3 (diff C_(i+1)^-)/(diff alpha) t^2 (1-t) + (diff P_(i+1))/(diff alpha) t^3 $

而由 $C_i^+ = P_i + V_i slash 3$、$C_i^- = P_i - V_i slash 3$：

$ (diff C_i^+)/(diff alpha) = (diff P_i)/(diff alpha) + 1/3 (diff V_i)/(diff alpha) $

$ (diff C_i^-)/(diff alpha) = (diff P_i)/(diff alpha) - 1/3 (diff V_i)/(diff alpha) $

$diff V_i slash diff alpha$ 由三对角系统的解对右端的导数给出：系数矩阵 $(1,4,1)$ 与
$alpha$ 无关，所以 $diff V slash diff alpha$ 是对同一个三对角系统换一个右端再解一次。
代码里没有缓存 LU 分解，而是直接把这个右端交给 `Bezier_tangent_vector()` 重解一遍
（Thomas 算法 $O(N)$，代价可以忽略）；因为求解关于右端是线性的，这样解出来的就是
$diff V_i slash diff alpha$。

由此得到的结论是：*$B(t)$ 对节点坐标 $P$ 的依赖是线性的*（精确地说，在"二阶导数为零"
端条件、或不带 `autoscale` 的"给定切矢"端条件下是线性的），所以那时
$diff B(t) slash diff P$ 与 $P$ 无关，可以在迭代前算好。`Bezier_V_sensitivity()`
就是干这个的（参数 `id_P`、`id_dim` 指定对哪个点、哪个分量求导，
`start_point_mode`/`end_point_mode`/`autoscale` 与前面一致）。

#fix[
  原稿写的是"$diff V_i slash diff alpha$ ……只需要用同一个 LU 分解再解一次"和
  "$diff B(t) slash diff P$ 与 $P$ 无关，可以在迭代前算好"，都没有附加条件。按代码，
  这两句只在端条件取 `mode = 1`（二阶导数为零）时才严格成立：那时三对角系统的右端
  是 $3(P_(i+1) - P_(i-1))$，$V$ 关于 $P$ 是线性的。本项目的 `infile.pot` 用的是
  `start_point_mode = 0`、`end_point_mode = 0` 加 `autoscale = true`，切矢的缩放系数
  里含 $|P_1 - P_0|$，于是 $V$（以及 $B(t)$）对 $P$ 不再是线性的，
  $diff B slash diff P$ 会随 $P$ 变化。`assemble_sensitivity_K()` 也正是这么做的：
  每次调用都重新算一遍 `Bezier_V_sensitivity(...)`，而不是复用迭代前的矩阵。
]

#code-ref[`Bezier_V_sensitivity()` 的构造方式是：先把 $diff P slash diff alpha$ 写成
  `P_sens`（只在 `(id_P, id_dim)` 处置 1），把端条件里对给定切矢的依赖也一并求出来，
  再调用 `Bezier_tangent_vector(P_sens, ..., /*autoscale=*/false)` 解同一个三对角系统，
  返回的就是 $diff V slash diff alpha$。端条件的两种情形在代码里是分开写的：

  - `mode = 0`（给定切矢）：`id_P` 是首点或次点时用
    `vec_start *= (P(id_P,id_dim) - P((id_P+1)%2, id_dim))/vec_norm/tan_vec_norm`，
    正对应 `autoscale` 里 $|P_1 - P_0| slash |V_"start"|$ 那一项对 $P_0, P_1$ 的导数；
    `id_P == N_points` 则是"端口切向角"这个额外设计变量（`cross_waveguide_2` 用
    `(id_P, id_dim) = (N_points, 0)` 表示起点、`(N_points, 1)` 表示终点），
    它的导数是把单位切矢旋转 $90degree$ 再乘上 `tan_vec_norm`。
  - `mode = 1`（二阶导数为零）：`vec_start` / `vec_end` 不参与——`Bezier_tangent_vector()`
    在这个分支里把右端写成 `rhs[0] = 3*(P.row(1) - P.row(0))`、`rhs[N-1] = 3*(P.row(N-1)
    - P.row(N-2))`，于是 $diff V slash diff alpha$ 完全由传进去的 `P_sens` 给出。

  `main.cpp` 的 `assemble_sensitivity_K()` 把它接在边界积分上：先
  `Sens_v = Bezier_V_sensitivity(point_id, dim_id, ..., true)`、`Sens_p(point_id, dim_id) = 1`，
  再用 `Bezier_interp(Sens_p, Sens_v, t_list, 0)` 求出该边界上每个节点的
  $diff B(t) slash diff p$ 与切矢灵敏度 `V_tan`，最后乘上界面两侧折射率差的权重
  （代码里的 `Delta_FEM_beta = -pow(omega,2)*epsilon_wg*mu_wg
  + pow(omega,2)*epsilon_amb*mu_amb`，见 `theory/02-maxwell-2d.typ` 的最后一节）
  和面积变化项
  `volume_elem = average_P_sens x average_V_tan` 组装成 `dKdp`。

  `Bezier.cpp` 在定义了 `EXPORT_PYLIB` 时还会被 pybind11 编译成 Python 扩展模块
  （`Bezier.h` 末尾的 `PYBIND11_MODULE` 导出 `Bezier_tangent_vector`、`Bezier_interp`、
  `Nearest_t`、`Bezier_V_sensitivity` 四个函数），于是 `python/original_mesh.py` 和
  `python/interp_DOF_to_points.py` 里的 `import Bezier` 用的就是和 C++ 完全相同的
  曲线求值/灵敏度代码，不会出现两边公式不一致的问题。本项目用 `make python` 生成它：
  `$(CXX) -O2 -std=gnu++17 -shared -fPIC -DEXPORT_PYLIB ... -o python/Bezier$(PYEXT)`，
  文件名后缀由 `python3-config --extension-suffix` 决定。]

#fix[
  原稿这里写的是 `Bezier.cpython-38-x86_64-linux-gnu.so`。那个文件确实是 3.8 编的
  （在 3.9 下 `import Bezier` 会报版本不匹配），现在的 `Makefile` 由
  `python3-config --extension-suffix` 自动取后缀，本机产出的名字是
  `python/Bezier.cpython-39-x86_64-linux-gnu.so`。
]
