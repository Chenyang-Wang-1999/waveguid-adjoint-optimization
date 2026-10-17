#import "../template.typ": project-doc, fig, fig-row, fig-gallery, note, fix, code-ref

#show: project-doc.with(
  title: "交叉波导伴随优化的结果",
  subtitle: "Adjoint shape optimization of a photonic cross waveguide",
  date: "2022-02",
)

这是本项目的第一项工作，也是后面所有工作（微环阵列的能带与边缘态，见 `02-ring-array.typ`）
的技术底座：把交叉波导交叉口的边界形状用 Bézier 曲线参数化，用有限元解一次正问题、
解一次伴随问题，就拿到目标函数对*全部*控制点坐标的梯度，然后做梯度下降。
优化前的透射率为 $|S_(21)|^2 = 0.712$，优化后达到 $approx 0.96$。

下面按「做什么 → 怎么做 → 结果如何 → 遇到什么问题 → 怎么解决」的顺序叙述。
伴随法灵敏度的完整推导不在这里重复，见 `theory/05-adjoint-sensitivity.typ`。

== 问题与模型

=== 交叉口为什么有损耗

两条单模介质波导以 $90 degree$ 相交。直通波导里的导模走到交叉口时，看到的是一段
横向被"撑开"、边界条件突然改变的波导，于是发生两件事：

- 场在四个角上被反射回输入端口，形成回波 $S_(11)$；
- 场向两侧的垂直波导散射，形成串扰 $S_(31)$，同时有一部分耦合到辐射模上损失掉。

在二维标量近似下，这两条通道之和就解释了为什么一个"什么都没有做"的交叉口
只有 $71.2%$ 的功率能到达直通端。优化的目标就是重新塑造交叉口四个角上的边界，
把被反射和散射掉的功率"引"回直通端。

=== 优化的几何对象

结构沿 $x$ 轴与 $y$ 轴镜像对称，四个端口彼此等价，因此只需要在四个角上各放一条同样形状的
曲线。每条曲线是四次 Bézier 曲线

$ B(t) = P_0 (1-t)^3 + 3 C_0^+ t (1-t)^2 + 3 C_1^- t^2 (1-t) + P_1 t^3 $

其中 $P_0, P_1$ 是曲线与波导边界的连接点，$C_0^+, C_1^-$ 是内部控制点。
曲线的端点条件由物理要求给定：

- 在端口（与均匀波导相接）处：切线方向给定，保证与直波导平滑衔接；
- 在节点（靠近交叉口中心）处：二阶导数为零，保证曲率不出现人为的跳变。

设计变量取这些内部控制点的坐标。交叉口以外是均匀波导段，形状固定、不参与优化。

#fig("assets/optimize/Trial_for_structure_20200223/geom.png",
  caption: [六边形仿真域（$6$ 个端口）与交叉口处的网格。蓝色是均匀波导，
    交叉口附近的深色区域是优化区，其内凹边界用 Bézier 曲线描述；
    网格在优化区与波导内被自适应加密，在包层里明显更粗。],
  width: 70%,
)

=== 参数

模型取 $epsilon_0 = mu_0 = 1$，于是整个系统只用长度单位 $l_"unit"$ 与电流单位
$I_"unit"$ 就能确定，所有几何量都以 $l_"unit"$ 度量。频率、介电常数等运行参数写在一个
GetPot 文本文件里：

#table(
  columns: (auto, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.35em),
  table.header([*参数*], [*取值*]),
  [$w_"wg"$], [$0.08 l_"unit"$（波导宽度）],
  [$nu$], [$1.55 slash t_"unit"$（工作频率）],
  [$L_"opt"$], [$0.3 l_"unit"$（可优化的交叉口区域尺寸）],
  [$epsilon_"wg"$], [$11.56$（波导芯，$n = 3.4$）],
  [$epsilon_"amb"$], [$1$（包层）],
  [$mu_"wg" = mu_"amb"$], [$1$],
)

#code-ref[`data/infile.pot` 是求解器读的全部运行参数，上面表中的材料参数逐项对应
`epsilon_wg = 11.560000`、`epsilon_amb = 1.000000`、`mu_wg = mu_amb = 1.000000`、
`freq = 1.550000`，以及网格/结果文件名。注意这个文件里的 `freq` 是*归一化频率*，
与 $nu = 1.55 slash t_"unit"$ 是同一个数；波导宽度与优化区尺寸则不在 `infile.pot` 里，
而是由 Python 侧的 `optimize.py` 与 `original_mesh.py` 在生成网格时写入几何
（见后文的参数一致性说明）。]

== 方法

=== 目标函数与伴随方程

设计变量向量记作 $bold(p)$（Bézier 控制点坐标），有限元离散后的复线性系统为

$ bold(K)(bold(p)) bold(a)(bold(p)) = bold(b) $

其中 $bold(a) = bold(x) + i bold(y)$ 是复场系数向量，$bold(b)$ 是端口入射模。
目标函数取直通端的透射率

$ W = |S_(21)|^2 $

引入伴随向量 $overline(bold(a))$（一个新的未知量，不是 $bold(a)$ 的共轭），令它满足

$ bold(K)^T overline(bold(a)) = (diff W)/(diff bold(x)) - i (diff W)/(diff bold(y)) $

于是梯度可以用一次"矩阵-向量-矩阵"的缩并写出来：

$ (diff W)/(diff p) = "Re" lr((overline(bold(a)) (diff)/(diff p) bold(K) bold(a))) $

#fix[原稿把梯度公式写成
$ (diff W)/(diff p) = "Re"(overline(bold(a)) (diff)/(diff p) bold(K) bold(a)) $，
*没有负号*。这里保留原稿的写法，但需要说明：按伴随方程的符号约定推下来，
结果应当是 $ (diff W)/(diff p) = - "Re"(overline(bold(a)) (diff)/(diff p) bold(K) bold(a)) $。
两种写法的差别只在伴随向量的符号约定，而这个符号当年确实折腾了很久——git 历史里有一次
提交就叫 “mysterious negative sign”。代码里的记账方式与带负号的版本一致：
`optimize.py` 把 C++ 输出的 `ex_opt_sens` 取负之后才当作 $W$ 的上升方向。
详细讨论见 `theory/05-adjoint-sensitivity.typ` 的「与 2022 年笔记的差异」一节，这里不再重复推导。]

=== 算法

原始记录给出的 6 步流程：

+ 用有限元算出电场 $bold(a)$；
+ 计算 $(diff W)/(diff bold(x))$ 与 $(diff W)/(diff bold(y))$；
+ 解伴随方程，得到与电场对应的伴随场 $overline(bold(a))$；
+ 对每个设计变量组装 $(diff)/(diff p) bold(K)$（*不需要解线性系统*）；
+ 计算 $W$ 对控制点坐标的偏导数；
+ 让控制点沿梯度方向移动。

#code-ref[这 6 步与代码是一一对应的。

- 第 1 步：`waveguide_optimization/main.cpp` 的 `main()`——
  `eq_systems.get_system("waveguide").solve()`。
- 第 2 步：`main.cpp` 的 `assemble_QOI_rhs()` 装配伴随方程的右端项，
  它在被选为目标端口的边界单元上累加
  `(2*std::real(std::conj(S)*dSdA) + 2*Complex(1j)*std::imag(std::conj(S)*dSdA))*JxW_face[qp]`；
  `main()` 里先把 `curr_S` 和 `QOI_port` 塞进 `EquationSystems::parameters`，
  再调用 `solver->adjoint_solve(...)`（libMesh 的转置求解）得到伴随场。
- 第 3 步：同一个 `adjoint_solve()` 调用。
- 第 4 步：`assemble_sensitivity_K()`——只遍历被优化边界
  （`id_ex_opt` / `id_in_opt`）上的单元，用 `Bezier_V_sensitivity()` 取控制点灵敏度，
  把形状导数组装进稀疏矩阵 `dKdp`。
- 第 5 步：`main()` 里的双重循环 `for(dim_id) for(point_id)`，
  `dKdp.vector_mult(middle_vec, solution)` + `middle_vec.conjugate()` +
  `adjoint_solution.dot(middle_vec)`，结果存进 `ex_opt_sens(point_id, dim_id)`
  与 `in_opt_sens(point_id, dim_id)`。
- 第 6 步：`python/optimize.py` 的 `gradient_descent()` 循环，
  沿 `sens` 方向以步长 `alpha` 移动控制点，并做非负截断。

外层循环的完整链路是
`optimize_cross_waveguide2()` → `gradient_descent(get_dof_derivative2, ...)`
→ `get_dof_derivative2()` → `cross_waveguide_2()`（`python/original_mesh.py`，用 Gmsh
重新建几何与网格）→ `run_solver()`（`python/paths.py`，跑 C++ 求解器）
→ 读回 `cross_waveguide_result.mat` 里的 `S_params_0` / `ex_opt_sens` / `in_opt_sens`
→ `points_to_DOF2()` 把两条边界的灵敏度拼回设计变量顺序。]

=== 设计变量的记账方式

#code-ref[设计变量不是"控制点矩阵"，而是一个扁平向量，转换由
`python/optimize.py` 的 `DOF_to_points2()` / `points_to_DOF2()` 完成：
`DOFs` 的前 $4 N_"ctrl" + 4$ 个分量按
$[x_("in,bd"), theta_"in", x_"in,1", y_"in,1", ..., y_("ex,bd"), theta_"ex", x_"ex,1", y_"ex,1", ...]$
排列，即两条边的内点坐标（各 $2 N_"ctrl"$ 个）加上两个边界点的自由坐标
（$x_("in,bd")$、$y_("ex,bd")$）和两个端点切角（$theta_"in"$、$theta_"ex"$）。
`gradient_descent()` 每接受一步就把 `sens`、`S_params`、`DOF` 追加写进
`sens_list.csv` / `S_params_list.csv` / `DOF_list.csv`，所以优化过程的透射率历史
是可以事后从 `S_params_list.csv` 里重建的。]

#note[边界上靠近端口的那一个控制点坐标由端口切向与固定的直波导长度决定，不参与优化，
所以 `get_dof_derivative()` 里写的是 `in_sens[1:,:]`、`ex_sens[1:,:]`（`get_dof_derivative2()`
里是 `in_sens[1:-1,:]` 加上端点切向 `in_sens[-1,1]`）。]

== 结果

=== 优化前后的透射率

在*未优化*的网格上跑一次正问题，四个端口的 S 参数为

$ S_"params,0" = [0.2843 + 0.0574 i, 0.6235 + 0.5686 i, 0.1543 - 0.2834 i, -0.0129 + 0.1721 i] $

因此

$ |S_(21)|^2 = |0.6235 + 0.5686 i|^2 = 0.712 $

这就是基线：一个几何上"干净"的交叉口，只有 $71.2%$ 的功率走直通端
（其余约为 $|S_(11)|^2 = 0.084$、$|S_(31)|^2 = 0.104$、$|S_(41)|^2 = 0.030$）。

优化后 $|S_(21)|^2 approx 0.96$，即透射率从 $0.712$ 提升到 $approx 0.96$（相对提升约 $35%$）。
按功率守恒 $sum_i |S_(i 1)|^2 <= 1$ 反推，反射与串扰之和此时已被压到 $4%$ 以内。

=== 优化后的场分布与宽带行为

#fig-row(
  (
    "assets/optimize/Trial_for_structure_20200223/electric_field_short.jpg",
    "assets/optimize/Trial_for_structure_20200223/transmission2.png",
  ),
  caption: [左：优化后结构中的电场分布（六边形仿真域的俯视图，色标为 $E^z$ 的幅值）。
    右：直通端与相邻端口的透射率随频率的变化（横轴为归一化搜索频率，纵轴为 $|S|^2$）。],
  width: 44%,
)

左侧的场分布说明优化确实起了作用：直通方向的场被连续地"接"了过去，交叉口处没有
明显的驻波斑图（那会是强反射的标志）。右侧的曲线则说明这个结果是*宽带*的而不是
在 $nu = 1.55 slash t_"unit"$ 上凑出来的——高透射率在一个相当宽的频率区间内保持平坦，
而不是一根尖峰；对级联器件来说，这一点比单频数值更重要。

#note[右图的横轴标注为 `f_search`，其数值范围（$1.5$–$1.6$）与 $nu = 1.55 slash t_"unit"$
同一量级，因此可以读作归一化频率；这里不把它解释为"优化迭代次数"。
优化迭代的透射率历史以 `S_params_list.csv` 的形式存在，但那张曲线在原始材料里没有出图。]

== 遇到的问题与改进

=== 原结构的节点区域太小

第一版几何里，交叉口中心附近"可动"的区域非常有限：可优化的只是端口与节点之间很短一段
边界上的 Bézier 控制点（`python/optimize.py` 的 `optimize_cross_waveguide()` 里
$L_"opt"$ 只有几个 $w_"wg"$ 量级），形状能变化的空间本来就小；同时中心附近"放大"的
网格分辨率也不足以分辨形状改变带来的场变化。结果是优化在很低的水平上就停滞了。

#fig("assets/figures/optimize_for_0.jpg",
  caption: [原环间耦合结构的场分布。场主要集中在两端的端口边界附近，
    中心节点区域几乎看不到场，说明可优化的自由度既不足、也没有被场"用到"。],
  width: 62%,
)

=== 增大耦合距离、给最近邻耦合加环路

改进的方向有两个：

- *增大环之间（以及环与波导之间）的耦合距离*，让耦合长度变长、耦合过程更缓变，
  从而给形状优化留出可操作的余地，也降低对网格分辨率的苛求；
- *在结构里加入环路*，用环路的长度去"兑现"最近邻耦合，把原本集中在很短一段上的
  相位积累摊开。

#fig("assets/figures/Long_coupling.svg",
  caption: [改进后的几何示意：中心是由六段波导围成的环路（绿色），
    六个蓝色矩形是输入/输出波导所在的位置，红色圆标注的是环与波导之间的最近邻耦合区。
    环路的作用是把最近邻耦合的相位积累摊到更长的路径上，同时把耦合区从节点处"挪"开。],
  width: 58%,
)

=== 改进结构的场分布

#fig("assets/figures/peak_field.svg",
  caption: [改进结构在谐振附近的场分布（ParaView 导出的 $E^z$ 幅值）。],
  width: 62%,
)

从场分布看，改进后的结构里场确实进入了中心区域、并在环路与耦合区上形成了有结构的
分布，而不是像原结构那样堆在端口边界上。这正是把"节点区域太小"这个问题解决掉之后
应该看到的样子。

=== 优化前后的 S 参数

#fig-row(
  (
    "assets/figures/S_parameters_amplitude_sq.png",
    "assets/figures/S_parameters_phase.png",
  ),
  caption: [改进结构 S 参数的幅度平方（左）与相位（右）随失谐量（Detune）的变化；
    蓝色为 $S_(21)$，橙色为 $S_(31)$。],
  width: 46%,
)

#note[原稿把这两张图并排放在同一页，没有给图注。上面的图注是按图内坐标轴与图例写的。
两张图的横轴都是 `Detune`，因此它们描述的是结构在谐振点附近的行为，
而不是 $|S|^2 = 0.96$ 那一组交叉波导的全景数据；具体哪个结构、哪个端口对应哪条曲线，
原始材料没有说明，这里不做推断。]

=== 后续设想：用低串扰波导测耦合系数

#fig("assets/optimize/Trial_for_structure_20200223/geom_double_rings_coupling.png",
  caption: [测量耦合系数的设想结构：两个环（或环阵列的相邻元胞）之间通过一段低串扰的
    直波导间接耦合，用两端测得的 S 参数反推耦合系数。],
  width: 54%,
)

== 参数与代码的一致性

写这份结果的时候做了一次逐项核对，发现两处需要说明的地方。

#fix[原稿记录的可优化区尺寸是 $L_"opt" = 0.3 l_"unit"$，而 `python/optimize.py` 的
`optimize_cross_waveguide2()` 里写的是 `L_opt = 0.2`（再由
`L_sim = (L_opt + L_fixed) slash cos(pi\/6)` 换算出仿真域尺寸，
`L_fixed = 0.3` 是固定的直波导段长度）。正文表中保留原稿的 $0.3 l_"unit"$，
因为那是本页结果所用的几何；但如果要照现在的代码重跑，得到的会是 $L_"opt" = 0.2$ 的几何。
两者不一致，且原始材料里没有说明哪一个是最终版本，这里不做静默修改。]

#note[`optimize_cross_waveguide()`（第一版几何）里的 `L_fixed = 0.3` 与 $L_"opt" = 0.3$
容易混淆，但两者是不同的量：原稿用 $L_"opt"$ 表示*可优化区*的尺寸，而 `optimize.py` 的
`L_fixed` 表示从端口往里*不可动*的直波导段长度。正文表中已按这个含义标注。]

#note[网格自适应加密的轮数来自 `main()` 里的 `do { ... } while(whether_continue)` 循环，
判据是 `EM_mesh_flagging` 给出的单元标记（`EM_mesh_flagging.{h,cpp}`）。
实测在 `data/cross_waveguide.msh` 上跑 3 轮后不再有单元被标记，于是退出循环；
`README.md` 里 `make run` 的典型输出 `Curr refinement: 1 / 2 / 3` 也是这个意思。
加密之后才做 `mesh_snap_to_Bezier()`：它把优化边界上的每个节点按
`Nearest_t()` 找到对应的曲线参数 $t$，再用 `Bezier_interp()` 把它吸附到曲线上，
并把 $t$ 与象限符号存进节点的 extra datum，供后面 `assemble_sensitivity_K()` 取用。]

#code-ref[灵敏度矩阵的形状是 `ex_opt_sens` / `in_opt_sens` 各 $13 times 2$：
`data/cross_waveguide.mat` 里的 `in_Bezier_P`、`ex_Bezier_P` 都是 $12 times 2$
（$12$ 个控制点，二维坐标），而 `infile.pot` 里 `start_point_mode = 0`、
`end_point_mode = 0`，于是 `main()` 走的是
`ex_points = ex_Bezier_P.rows() + 1` 这一支，得到 $13$ 行；
列数 $2$ 就是坐标的 $x$、$y$ 两个方向。多出来的那一行对应端点切向的自由度。]

== 下一步

这份工作验证了"有限元 + 伴随法 + Bézier 参数化"这条路线是可行的：一次伴随求解就能给出
全部控制点的梯度，优化在单机上可以跑几十轮。下一步是把同样的方法用到
*微环阵列*上——那里耦合系数是最关键的参数，需要先能准确地测出（或算出）它，
再谈优化。这就是 `geom_double_rings_coupling.png` 那张设想图的来由，
也是 `02-ring-array.typ` 的主题。

== 参考文献

- H. Igarashi and K. Watanabe, “Complex adjoint variable method for finite element analysis
  of nonlinear …”, _IEEE Transactions on Magnetics_, *46*(8), 2739–2742, 2010.
- 推导细节：`docs/theory/05-adjoint-sensitivity.typ`；
  Bézier 参数化：`docs/theory/04-bezier.typ`；
  二维麦克斯韦方程与形状导数：`docs/theory/02-maxwell-2d.typ`。
