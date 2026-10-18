#import "../template.typ": project-doc, fig, note, fix, code-ref

#show: project-doc.with(
  title: "伴随法灵敏度分析与 Bézier 控制点梯度",
  subtitle: "Adjoint sensitivity analysis for the discrete FEM system",
  date: "2022-02",
)

这一节是整个项目的核心：怎么在只解一次线性系统的前提下，求出目标函数对*全部*设计变量
（Bézier 控制点坐标）的梯度。方法本身来自 H. Igarashi 和 K. Watanabe,
_IEEE Trans. Magn._ *46*, 2739 (2010)。

== 问题设定

有限元离散之后，电磁场问题化成一个线性代数方程组

$ bold(K)(bold(p)) bold(a)(bold(p)) = bold(b) $

其中：

- $bold(p) in RR^N$ 是设计变量向量（这里是 Bézier 曲线控制点的坐标，$N = 4 N_"ctrl" + 4$，
  取 $N_"ctrl" = 10$ 时 $N = 44$）；
- $bold(K)$ 是复对称（非 Hermite）的刚度矩阵，由边界上的积分组装而成，因此*依赖于形状*，
  也就是依赖于 $bold(p)$；
- $bold(b)$ 是右端项（端口上的入射模），端口本身不参与优化，所以 $bold(b)$ 与 $bold(p)$ 无关；
- $bold(a) = bold(x) + i bold(y)$ 是待求的复场系数向量，$bold(x), bold(y) in RR^n$。

目标函数取直通端的透射率

$ W = |S_(21)|^2 $

$S_(21)$ 是端口上的模式重叠积分（见 `overview.typ`），所以 $W$ 通过 $bold(a)$ 间接依赖 $bold(p)$。

我们来求 $diff W slash diff p$。

== 为什么不能直接求导

最直接的做法是把 $bold(K) bold(a) = bold(b)$ 对 $p$ 求导：

$ bold(K) (diff bold(a))/(diff p) = - (diff bold(K))/(diff p) bold(a) $

然后

$ (diff W)/(diff p) = Re( (diff W)/(diff bold(a)) dot (diff bold(a))/(diff p) ) $

这条路线对*每一个*设计变量都要解一次线性系统。$N = 44$ 时，每次迭代要做 44 次全场求解——
比正问题本身还贵 44 倍。伴随法的意义就在于把这 44 次压到 1 次。

== 伴随方程

引入*伴随向量* $overline(bold(a))$（注意这是一个新的未知量，不是 $bold(a)$ 的共轭），
让它满足

$ bold(K)^T overline(bold(a)) = (diff W)/(diff bold(x)) - i (diff W)/(diff bold(y)) $

这里 $bold(K)^T$ 是*转置*（不是共轭转置），对应的正是 libMesh 里
`LinearSolver::adjoint_solve()` 做的事——它的实现就是把矩阵转置一下再调用普通的 `solve()`
（见 `libmesh/src/solvers/linear_solver.C`）。

关键点是：*伴随方程的右端只跟目标函数有关，跟设计变量的个数无关*。所以不管有 44 个还是
440 个设计变量，伴随方程永远只解一次。

#code-ref[这两件事在 `waveguide_optimization/main.cpp` 的 `main()` 里各占一行：

  ```cpp
  assemble_QOI_rhs(eq_systems, "waveguide");   // 装配伴随方程右端 g = dW/dx - i dW/dy
  LinearSolver<Number> * solver = curr_system.get_linear_solver();
  solver ->adjoint_solve(curr_system.get_system_matrix(),
                         curr_system.add_adjoint_solution(0),
                         curr_system.get_adjoint_rhs(0),
                         solver_params.second, solver_params.first);
  ```

  伴随向量存在 `get_adjoint_solution(0)` 里，右端存在 `get_adjoint_rhs(0)` 里（由
  `assemble_QOI_rhs()` 往里面累加），矩阵参数直接复用正问题的 `get_system_matrix()`；
  `solver_params` 来自 `curr_system.get_linear_solve_parameters()`，它的
  `.first / .second` 分别是最大迭代数和收敛容差，所以上面传参的顺序是
  `(tol, maxits)`。`adjoint_solve()` 解的是转置系统——这正是本文后面用
  `middle_vec.conjugate()` 手工还原的约定（见"算法流程"一节）。]

== 梯度公式

=== 推导

由 $bold(K) bold(a) = bold(b)$ 在 $bold(b)$ 与 $p$ 无关时求导：

$ (diff bold(K))/(diff p) bold(a) = - bold(K) (diff bold(a))/(diff p) $

两边左乘 $overline(bold(a))^T$，并用 $overline(bold(a))^T bold(K) = bold(g)^T$，其中
$bold(g) = diff W slash diff bold(x) - i diff W slash diff bold(y)$，得到

$ overline(bold(a))^T (diff bold(K))/(diff p) bold(a) = - bold(g)^T (diff bold(a))/(diff p) $

而 $bold(a) = bold(x) + i bold(y)$，所以

$ bold(g)^T (diff bold(a))/(diff p) = (diff W)/(diff bold(x))^T (diff bold(x))/(diff p) + (diff W)/(diff bold(y))^T (diff bold(y))/(diff p) + i[:] $

其中方括号里的项是纯虚的，取实部就剩下全微分：

$ Re( bold(g)^T (diff bold(a))/(diff p) ) = (diff W)/(diff p) $

代回去就得到

$ (diff W)/(diff p) = - Re( overline(bold(a))^T (diff bold(K))/(diff p) bold(a) ) $

=== 与 2022 年笔记的差异

#note[
  原始笔记（`PDE_with_libmesh.lyx` / 2022-02-23 周报）把这一条写成
  $ (diff W)/(diff p) = Re(overline(bold(a)) (diff bold(K))/(diff p) bold(a)) $，
  *没有负号*。上面的推导给出的是带负号的版本；两者的差别在于伴随向量（或等价的伴随右端）
  的符号约定。代码里的记账方式与带负号的版本一致：
  `optimize.py` 把 C++ 输出的 `ex_opt_sens` 取负之后才当作 $W$ 的上升方向
  （`get_dof_derivative` 里的 `return(-DOF_sens.flatten(), ...)`，循环里
  `DOF += sens * alpha`）。
  当年这个符号确实折腾了很久——git 历史里有一次提交就叫 “mysterious negative sign”。
  这里保留两种写法，是为了不让整理后的文档显得比当年更确定。
]

=== 右端项的具体形式

目标函数 $W = |S|^2 = S S^*$，其中 $S$ 是端口上的模式重叠积分
$S = integral_Gamma bold(D)_"mode"^* dot bold(E) slash 4 dif Gamma$，
而 $bold(E)$（也就是 $E^z$ 及其梯度）线性依赖于 $bold(a)$。于是

$ (diff W)/(diff x) = 2 Re( S^* (diff S)/(diff x) ), quad (diff W)/(diff y) = 2 Re( S^* (diff S)/(diff y) ) $

代入 $bold(g) = diff W slash diff bold(x) - i diff W slash diff bold(y)$：

$ bold(g) = 2 Re( S^* (diff S)/(diff A) ) + 2 i Im( S^* (diff S)/(diff A) ) = 2 S^* (diff S)/(diff A) $

（用了 $diff S slash diff x = diff S slash diff A$、$diff S slash diff y = i diff S slash diff A$，
其中 $A$ 是场系数；最后一步就是 $z = "Re" z + i "Im" z$，取 $z = S^* diff S slash diff A$。）

#fix[
  原稿这一步写作 $= 2 S (diff S)/(diff A)^*$，也就是 $2 overline(overline(S) (d S)\/(d A))$。
  按前一步，$2 "Re"(S^* S_A) + 2i "Im"(S^* S_A) = 2 S^* S_A$（复共轭标错了位置：
  $overline(overline(S) S_A) = S overline(S_A) != S^* S_A$，两者差一个共轭）。
  代码里 `2*std::real(std::conj(S)*dSdA) + 2*Complex(1j)*std::imag(std::conj(S)*dSdA)`
  算的正是 $2 S^* (d S)\/(d A)$，与改正后的写法一致。
]

代码里 `assemble_QOI_rhs()` 正好就是这一项：

```cpp
Number dSdA = (std::conj(mode_Dz)*phi_face[i][qp]
             + std::conj(mode_Hx)*dBxdA
             + std::conj(mode_Hy)*dBydA)/4;
Fe(i) += (2*std::real(std::conj(S)*dSdA)
        + 2*Complex(1j)*std::imag(std::conj(S)*dSdA))*JxW_face[qp];
```

`2*real(...) + 2i*imag(...)` 就是 $2 S^* (d S)\/(d A)$，和上面的推导逐项对上。
这一段只在被选为目标端口的边界单元上装配。

== 从边界变分到 $diff bold(K) slash diff p$

$diff bold(K) slash diff p$ 不需要真的把矩阵求导。$bold(K)$ 的每个元素都是形如

$ K_(i j) = integral_Omega ( alpha_x (diff N_i)/x (diff N_j)/x + alpha_y (diff N_i)/y (diff N_j)/y + beta N_i N_j ) dif Omega $

的积分。当边界发生 $delta bold(r)$ 的变形时，只有*被优化边界附近的单元*形状改变，
积分区域随之变化；把形状导数写出来就是

$ delta F(phi) = 1/2 integral_("deformed boundary") [ (alpha_x^((2)) - alpha_x^((1))) (diff phi/x)^2 + (alpha_y^((2)) - alpha_y^((1))) (diff phi/y)^2 + (beta^((2)) - beta^((1))) phi^2 ] dif Gamma $

也就是说：*形状导数是一个沿界面的一维积分*，而"沿界面对形状函数积分"得到的正是
$diff bold(K) slash diff p$（推导见 `theory/02-maxwell-2d.typ` 的最后一节）。

实现上，`assemble_sensitivity_K(system, params, dKdp, point_id, dim_id, ...)` 做的事就是：

+ 只遍历被优化边界（`boundary_id`，也就是 `in_opt` / `ex_opt`）上的单元；
+ 用 `Bezier_V_sensitivity(point_id, dim_id, P, ...)` 算出该边界节点坐标对
  第 `(point_id, dim_id)` 个设计变量的导数 $diff B(t) slash diff p$；
+ 把边界积分对 $p$ 的导数组装成一个稀疏矩阵 `dKdp`。

因为 $B(t)$ 对控制点的依赖是线性的（见 `theory/04-bezier.typ`；严格说只在端条件取
"二阶导数为零"时成立，本项目用的是 `start_point_mode = 0` 加 `autoscale`，所以
`assemble_sensitivity_K()` 每次调用都重新算一遍 `Bezier_V_sensitivity(...)`，
而不是在迭代前算好复用），$diff B slash diff p$ 可以随 $p$ 一起更新。

== 算法流程

原文给出的 6 步（2022-02-23 周报）：

+ 用有限元算出电场 $bold(a)$；
+ 计算 $diff W slash diff bold(x)$ 与 $diff W slash diff bold(y)$；
+ 解伴随方程，得到与电场对应的伴随场 $overline(bold(a))$；
+ 对每个设计变量组装 $diff bold(K) slash diff p$（*不解线性系统*）；
+ 计算 $W$ 对控制点坐标的偏导数；
+ 让控制点沿梯度方向移动。

对应到代码（`main.cpp` 的 `main()`，行内注释里的编号就是上面的步号）：

```cpp
// 1. 正问题
eq_systems.get_system("waveguide").solve();

// 2. S 参数：四个端口上的模式重叠积分
Eigen::VectorXcd S_params(4);
std::vector<subdomain_id_type> selected_domains = {3,4,5,6};
integrate_along_side(eq_systems, selected_domains, S_params, mode_decompose);

// 目标端口：output 1；W = |S|^2
unsigned int selected_port = 1;
eq_systems.parameters.set<Complex>("curr_S") = S_params(selected_port);
eq_systems.parameters.set<subdomain_id_type>("QOI_port") = selected_port;

// 3. 伴随右端 + 伴随解
assemble_QOI_rhs(eq_systems, "waveguide");
LinearSolver<Number> * solver = curr_system.get_linear_solver();
solver->adjoint_solve(curr_system.get_system_matrix(),
                      curr_system.add_adjoint_solution(0),
                      curr_system.get_adjoint_rhs(0),
                      solver_params.second, solver_params.first);

// 4-5. 逐设计变量组装 dK/dp，做矩阵-向量乘
for (unsigned int dim_id = 0; dim_id < ex_Bezier_P.cols(); dim_id++)
  for (unsigned int point_id = 0; point_id < ex_points; point_id++) {
    middle_vec.zero();
    curr_system.get_matrix("dKdp").zero();
    assemble_sensitivity_K(eq_systems.get_system("waveguide"), eq_systems.parameters,
                           curr_system.get_matrix("dKdp"), point_id, dim_id,
                           start_point_mode, end_point_mode,
                           id_ex_opt, ex_Bezier_P, ex_Bezier_V);
    curr_system.get_matrix("dKdp").vector_mult(middle_vec,
        * eq_systems.get_system("waveguide").solution);
    middle_vec.conjugate();
    ex_opt_sens(point_id, dim_id) = std::real(
        eq_systems.get_system("waveguide").get_adjoint_solution(0).dot(middle_vec));
  }
// 入口边界 in_opt 用同一段循环再走一遍（dKdp2 / in_opt_sens）
```

第 6 步在 C++ 侧没有：灵敏度写进 `result.h5` 的 `/scalars/ex_opt_sens` /
`/scalars/in_opt_sens` 之后，由 Python 侧的优化循环接手。

#code-ref[6 步与函数的对应关系：

  #table(
    columns: (auto, 1fr),
    align: (left, left),
    stroke: none,
    inset: (x: 0.6em, y: 0.3em),
    [*步*], [*代码位置*],
    [1. 正问题],
    [`main.cpp` 的 `main()`：`eq_systems.get_system("waveguide").solve()`，
      装配函数是 `assemble_waveguide()`],
    [2. $diff W slash diff bold(x) - i diff W slash diff bold(y)$],
    [`assemble_QOI_rhs(eq_systems, "waveguide")`，只遍历 `QOI_port` 指定的那个端口的边界单元],
    [3. 伴随解],
    [`solver->adjoint_solve(...)`，右端是 `get_adjoint_rhs(0)`，解是 `get_adjoint_solution(0)`],
    [4. $diff bold(K) slash diff p$],
    [`assemble_sensitivity_K(...)`，每次只组装一个设计变量对应的稀疏矩阵 `dKdp`
      （入口边界用 `dKdp2`）],
    [5. $diff W slash diff p$],
    [`dKdp.vector_mult(middle_vec, solution)` 做矩阵-向量乘，`middle_vec.conjugate()`
      还原转置（不是共轭转置）约定，再与伴随解点乘取实部],
    [6. 移动控制点],
    [`python/optimize.py`：`gradient_descent(my_der_fun, DOF0, params, alpha)` 里
      `DOF += sens * alpha`；梯度由 `get_dof_derivative()` 返回 `-DOF_sens.flatten()`，
      `points_to_DOF` / `points_to_DOF2` 把 `in_opt_sens` / `ex_opt_sens`
      拼回设计变量顺序],
  )

  两处边界各有一份灵敏度矩阵：出口边界用 `dKdp` + `ex_opt_sens`，入口边界用
  `dKdp2` + `in_opt_sens`，循环体完全一样。]

注意 `middle_vec.conjugate()` 这一步：libMesh 的 `NumericVector::dot(v)` 在复数情形下
对*第二个*参数取共轭，所以先共轭再点乘得到的是 $overline(bold(a))^T (diff bold(K) slash diff p) bold(a)$，
即纯粹的转置形式，与前面的推导一致。

== 代价

设 $n$ 是自由度个数，$N$ 是设计变量个数：

#table(
  columns: (auto, auto, auto),
  align: (left, left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.35em),
  table.header([*做法*], [*线性求解次数*], [*每步额外代价*]),
  [有限差分], [$N + 1$], [——],
  [伴随法], [$2$（正问题 + 伴随）], [$N$ 次 $diff bold(K) slash diff p$ 组装 + $N$ 次矩阵-向量乘],
)

$N = 44$ 时，伴随法把线性求解从 45 次降到 2 次；而多出来的 $N$ 次矩阵-向量乘
（$O(N n)$，稀疏）和 $N$ 次边界组装（只涉及边界单元）相对便宜得多。这就是这个项目
能在单机上跑几十轮优化迭代的原因。

== 验证

`python/test_sensitivity.py` 用有限差分校核伴随梯度：同一个结构算一次 $W$ 和灵敏度，
把设计变量整体随机扰动后再算一次，然后比较方向导数与差分

$ sum_i (diff W)/(diff p_i) Delta p_i quad "与" quad |S_("out") (bold(p) + Delta bold(p))|^2 - |S_("out") (bold(p))|^2 $

```python
DOF_sens = points_to_DOF(in_sens[1:,:],ex_sens[1:,:])
return(DOF_sens, S_params)
...
DOFs += np.random.rand(DOFs.shape[0]) * 0.01
Sens, S_params = get_dof_derivative(DOFs, (0.1,))
...
diff_FD = np.abs(data['all_S_params'][1,1])**2 - np.abs(data['all_S_params'][0,1])**2
diff_sens = data['all_sens'].dot(data['all_DOFs'][1,:] - data['all_DOFs'][0,:])
print(diff_sens)
print(diff_FD)
```

#fix[
  原稿引的代码是"角度分量单独用小步长 $2 pi times 0.01$、其余分量按 $0.01 w_"wg"$ 扰动"那一段
  `delta_DOF`，并把它记在 `test_sensitivity.py` 名下。那段代码实际在
  `python/optimize.py` 的 `test_sensitivity()` 里（它用的是带切向角自由度的
  `cross_waveguide_2` / `get_dof_derivative2`，所以两个角度自由度要用
  `2*pi*0.01` 这种量级）。独立的 `python/test_sensitivity.py` 是它的精简版：
  所有自由度一起加 `np.random.rand(...)*0.01`，再比较
  `diff_sens = sens dot Delta DOF` 与 `diff_FD`（`S_params[1]` 即出口 1，
  也就是 $S_(21)$）。这里按当前脚本改写。
]

#code-ref[校核脚本的调用链：`python/tests/test_fd_sensitivity.py` 先
  `design.make_case(dofs)` 生成网格与 `design.h5`，再
  `solver.run()`（`paths.py` 里拼出 `mpiexec -n 1 <solver> -pc_type lu`）跑一次 C++ 求解器，
  然后从 `result.h5` 里读 `/scalars/S_params`、`ex_opt_sens`、`in_opt_sens`，
  用 `design.sensitivity_to_dofs(in_sens, ex_sens)` 拼回设计变量顺序。
  `optimizer.dof_derivative()` 走的是同一条链，只是多返回一个负号
  （因为 `gradient_descent()` 的更新式是 `DOF += sens * alpha`）。

  注意微扰之后必须**重划网格**（`make_case` / `build_mesh`）：网格是设计变量的函数，
  只改 `design.h5` 不重划，曲线变了而网格没变，物理量根本不会动。]

#note[
  代码里 `points_to_DOF` / `points_to_DOF2` 把 C++ 返回的 `in_opt_sens`、`ex_opt_sens`
  重新拼回设计变量顺序，注意它跳过了第 0 个点（`in_sens[1:,:]`）——因为边界上第一个
  控制点的坐标是固定的（由端口切向与 `L_fixed` 决定），不参与优化。这也是
  `start_point_mode` / `end_point_mode` 在算设计变量个数时要把首点去掉的原因
  （`main.cpp` 里 `ex_points = ex_Bezier_P.rows()` 那一支）。
]

== 参考文献

H. Igarashi and K. Watanabe, “Complex adjoint variable method for finite element analysis
of nonlinear …”, _IEEE Transactions on Magnetics_, vol. 46, no. 8, pp. 2739–2742, 2010.
（当年保存的 PDF 在 `FEM_em_wave/IEEETM46-8_2739-2742_adjoint.pdf`。）
