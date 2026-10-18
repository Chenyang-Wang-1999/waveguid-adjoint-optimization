#import "../template.typ": project-doc, note, fix, code-ref

#show: project-doc.with(
  title: "libMesh 使用笔记",
)

== 编译与链接

=== 编译注意点

+ libMesh 的安装信息统一由 `<libmesh-prefix>/bin/libmesh-config` 给出，本项目
  `Makefile` 就是这么做的（`LIBMESH_CONFIG := $(LIBMESH_DIR)/bin/libmesh-config`）：

  ```make
  CXX      := $(shell $(LIBMESH_CONFIG) --cxx)
  CXXFLAGS := $(shell $(LIBMESH_CONFIG) --cxxflags) -std=gnu++17 -Wall
  CPPFLAGS := $(shell $(LIBMESH_CONFIG) --include)
  LIBS     := $(shell $(LIBMESH_CONFIG) --libs) -lmatio
  ```

  这样切换 libMesh 版本/位置只要改 `LIBMESH_DIR`，不用手写任何 include/lib 路径。
  `libmesh-config` 认识的开关只有 `--cppflags --cxxflags --include --libs --cxx --cc
  --fc --fflags --version --host --ldflags`（不带参数运行会打印用法）。

+ 编译模式（`opt` / `dbg` / `devel` / `prof` / `oprof`，默认 `opt`）由*环境变量*
  `METHOD` 决定：`libmesh-config` 脚本里有一个 `case "$METHOD" in`，按它设置
  `CXXFLAGS` 并把库名后缀换成 `_opt` / `_dbg` / …（例如 `METHOD=dbg libmesh-config
  --libs` 给出 `-lmesh_dbg`）。所以要让 `METHOD` 真的生效，必须让它在
  `libmesh-config` 的环境里，例如 `METHOD=dbg make`，或在 Makefile 里 `export METHOD`。

+ 本机的 libMesh 是 `--enable-complex` 编的（`libmesh_config.h` 里
  `LIBMESH_USE_COMPLEX_NUMBERS` 为 1），`Number` 因此是 `std::complex<double>`。
  本项目求解的是复振幅场，依赖这一点；换用实数版的 libMesh 会直接编不过。

=== 与自己的项目链接

+ *（更新）* libMesh 从 0.8.0 起不再推荐用 `Make.common`：路径与编译规则都改由
  `bin/libmesh-config` 提供（见上一节）。安装包里仍然带着 `<libmesh-prefix>/Make.common`，
  但那是给旧项目用的向后兼容入口，它自己的文件开头就写着
  "For backwards-compatible builds (before v0.8.0) you can include this Make.common file
  in your project. Note this is not the preferred method for using libMesh -- we provide
  libmesh-config as a better solution."（本机 `~/local/Make.common` 第 12–15 行）。
  只有维护很老的 Makefile 时才需要 include 它。

  #fix[
    原稿把 `Make.common` 写成与 libMesh 链接的唯一方式（"libmesh 编译所需的路径及编译规则
    定义在 `<libmesh-prefix>/Make.common` 中……在 Makefile 中引用 `<libmesh-dir>/Make.common`，
    然后使用其中的变量即可"）。这个说法现在过时了，本项目用的是 `libmesh-config`。
    下面两条关于 `obj-suffix` 和 `METHOD` 的知识点仍然成立，保留并补上出处。
  ]

+ *（保留的知识点）* 用 `Make.common` 时，需要由 libMesh 规则编译的源文件，其 object
  文件名必须是 `foo.$(obj-suffix)`——原稿写成"对应的 object 后缀必须是 `.$(obj-suffix)`"。
  `obj-suffix` 在 `Make.common` 里的定义是
  `obj-suffix := x86_64-pc-linux-gnu.$(METHOD).o`（本机 `~/local/Make.common` 第 127 行），
  配套规则是 `%.$(obj-suffix) : %.C`，所以 `foo.C` 会编成
  `foo.x86_64-pc-linux-gnu.opt.o`；这样不同 `METHOD` 的目标文件可以并存在同一个目录里。
  换成 `libmesh-config` 之后没有这个约束，object 名字随便取
  （本项目就是普通的 `main.o`、`Bezier.o` ……）。

+ libmesh 的编译有多种模式，如 `dbg, opt, devel` 等。可供选择的模式由 libmesh 编译时的
  configure 选项决定。详见 libmesh 的 README.md。编译时，定义变量
  `METHOD=<your method>` 即可（用 `Make.common` 时它是 make 变量，用 `libmesh-config`
  时它是环境变量，见"编译注意点"）。

== 基本模块用法

=== Mesh, Element and Nodes

+ Mesh 是整个网格，其中每个点就是一个 Node

+ Element 是有限元单元，用于近似表示一个 Node 邻域上的函数值

=== Context

Context 是 libmesh 用于提取 Element 信息的工具

+ Context 对每个 element 需要初始化，才可以获取数据 `pre_fe_init(const System &, const Elem *)`

+ 数据取自 `system.current_local_solution`

+ Context 工作原理如下：

  + 用户指定需要计算的变量。调用 Context.get\_xxx 函数时，Context 会将对应的标志置 1，然后在后续的初始化中计算

  + pre\_fe\_reinit 与 xx\_fe\_reinit：当用户需要计算某一特定单元中的物理量时，通过 pre\_fe\_reinit 将这个单元传给 Context，再通过 xx\_fe\_reinit 初始化对应的 elem 或 side

  + reinit 函数的工作原理是，通过用户给定的 flag，计算单元内对应的量。默认计算积分点上的权重。

=== 杂项与使用技巧

==== 类继承之间的指针与引用转换

libmesh 使用了大量虚拟类。为了解决类继承时的指针与引用的转换，libmesh 提供了许多转换函数

+ `cast_ref`：用于转换不同类之间的指针，即 `Type1 & type1_ref = cast_ref<Type2 &>(type2_ref)`

==== Point

+ Point 继承自模板类 `TypeVector`

+ 从 libMesh 里取出来的点列（`fe->get_xyz()`、`fe_face->get_xyz()` 等）要用
  `const Point &` 引用，不要整个拷贝一份。原稿由此写成"不能创建 `Point` 类对象，
  只能用 `const Point &` 指向它"，但项目自己的代码就直接构造了 `Point`：
  `N_poav_integrand()` 与 `assemble_sensitivity_K()` 里都有
  `Point averged_normal_vec(0,0,0);`（求边界法向的平均方向），编译运行都正常——
  所以这条限制并不存在，"不要拷贝"才是本意。

  #fix[
    原稿这一条写作"不能创建 `Point` 类对象，只能用 `const Point &` 指向它"。
    与项目代码不符（`main.cpp` 里直接构造了 `Point` 对象），按上面改写。
  ]

==== 初始化

+ Dofmap 等变量需要在系统初始化之后再设置，否则会报段错误。本项目里所有的自由度操作
  （`dof_map.dof_indices()`、`variable_type()`、`constrain_element_matrix_and_vector()`）
  都发生在 `eq_systems.init()` 之后的装配阶段，`init()` 之前只做
  `add_variable()` / `attach_assemble_function()` / `parameters.set()` 这些设置。

== 含时问题求解

=== 含时变量的确定

`libmesh::DifferentiablePhysics` 类有一方法 `time_evolving(unsigned int var, unsigned int order)`，用于指定变量的含时性。若该变量阶数被指定为 $k$，则认为

$ (dif^k u)/(dif t^k) = F(u, v, dots) $

特殊地，当 $k = 0$ 时，认为变量不含时，对应于约束项，或

$ 0 = G(u, v, dots) $

案例见 `examples/fem_system/fem_system_ex1`

=== 求解的方程

内置的迭代算法求解的方程为

$ integral dif bold(x) phi_i (bold(x)) (dif)/(dif t) u(bold(x)) = R_i (bold(u)_j) = integral dif bold(x) phi_i (bold(x)) F(u(bold(x))) \
  = integral dif bold(x) phi_i (bold(x)) F(phi_j (bold(x)) u_j) $

=== 设置初值：

只要设置 `system.solution` 与 `system.current_local_solution` 两个向量即可。注意全局向量设置完要调用 `vec.close()` 进行装配

== 本项目里的用法（代码对应）

上面这些 API 在交叉波导项目里的落点（除特别注明外都在
`waveguide_optimization/main.cpp`）：

#code-ref[

  - *`LinearImplicitSystem`*：`main()` 里 `eq_systems.add_system<LinearImplicitSystem>("waveguide")`，
    `add_variable("Ez", SECOND)`，再 `attach_assemble_function(assemble_waveguide)`；
    装配函数里用 `system.get_system_matrix()` 与 `system.rhs` 累加单元矩阵/右端，
    求解直接用 `eq_systems.get_system("waveguide").solve()`。
  - *`FEMContext`*：`main()` 里只构造了一个
    `FEMContext context(eq_systems.get_system("waveguide"), {Ez_id})`（指定要算的变量），
    真正的取数在 `integrate_along_side()` 里：对每个单元依次
    `context.pre_fe_reinit(sys, elem)`、`context.elem_fe_reinit()`，
    再对每条边 `context.side = side; context.side_fe_reinit()`，然后用
    `context.has_side_boundary_id(...)` 判断边界、`context.side_value(0, qp)` /
    `context.side_gradient(0, qp)` 取积分点上的值和梯度——正是"先声明变量、再 reinit、
    最后按 flag 取量"的用法。`mode_decompose()`、`N_poav_integrand()` 拿到的就是
    这个 context。
  - *网格不再由 libMesh 加密*：早期版本用 `MeshRefinement` + 自定义 flagging
    （`EM_mesh_flagging.{h,cpp}` 的 `EMMeshFlagging`，逐个
    `elem->set_refinement_flag(Elem::REFINE)`）反复二分加密；这一套已经在"网格前移"
    里删掉了，网格改由 Python 侧的 gmsh 尺寸场一次生成。
  - *`allow_renumbering(false)`*：`mesh.read()` 内部会调 `prepare_for_use()`，
    默认按**单元遍历顺序**重排节点 id；而 `design.h5` 里的 `/design/node_t` 是按
    `.msh` 文件里的节点顺序写的，所以读网格之前必须关掉重排（`ReplicatedMesh` 上
    这个开关走的是"只删孤立节点、不动编号"的分支）。
  - *`QGauss`*：`assemble_waveguide()` 里体单元用 `QGauss qrule(dim, FIFTH)`、
    边界用 `QGauss qface(dim-1, FIFTH)`（边界的维数比单元低一维），分别
    `fe->attach_quadrature_rule(&qrule)` / `fe_face->attach_quadrature_rule(&qface)`；
    积分权重取自 `fe->get_JxW()` 与 `fe_face->get_JxW()`，
    `assemble_QOI_rhs()` 与 `assemble_sensitivity_K()` 里也是同样的 `QGauss qface(dim-1, FIFTH)`。
  - *`Mesh` / `Elem` / `Node`*：`mesh.read(curr_params.mesh_file)` 读 gmsh 网格，
    `mesh.get_id_by_name("waveguide")` 与 `mesh_boundary.get_id_by_name("ex_opt")`
    从 `.msh` 的 PhysicalNames 里取子域/边界 id（`infile.pot` 里给了数字 id 时优先用数字）；
    `tag_boundary_param()` 遍历 `mesh.node_ptr_range()`（**所有**节点，不是 local），
    按节点 id 从 `/design/node_t` 取参数 $t$ 存进 `node->set_extra_datum(0, t)`，
    并逐节点校验 $norm(bold(x) - bold(P)(t))$；`assemble_sensitivity_K()` 再用
    `get_extra_datum<double>(0)` 读回来。象限符号不再存 datum，直接用节点坐标的符号。
]

#code-ref[*杂项*那几条在本项目里的情况：

  - 从 libMesh 取出的点：`assemble_waveguide()` 里
    `const std::vector<Point> & qface_point = fe_face->get_xyz();`（引用，不拷贝），
    转到本项目自己的接口时再显式转成 `Eigen::Vector3d pt = {qface_point[qp](0), qface_point[qp](1), 0.0};`
    传给 `WaveguideParams::analytic_E_field()`。
  - `DofMap`：`const DofMap & dof_map = system.get_dof_map();` 之后用
    `dof_map.variable_type(0)` 取 `FEType`、`dof_map.dof_indices(elem, dof_indices)`
    取单元自由度编号、`dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices)`
    在加进全局矩阵前施加悬挂节点约束（多轮 `MeshRefinement` 之后这一步是必需的）。
  - `cast_ref`：本项目没有用到——矩阵和向量都以基类引用传递
    （`SparseMatrix<Number> & dKdp`、`NumericVector<Number> & middle_vec`），
    多态由 libMesh 内部处理。
]
