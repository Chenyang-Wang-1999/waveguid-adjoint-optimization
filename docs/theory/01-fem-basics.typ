#import "../template.typ": project-doc, note, fix, code-ref

#show: project-doc.with(
  title: "有限元方法基础",
)

// 原文用 \eqref 引用了两个编号公式（eq:PDE 与 eq:lemma1），这里给公式编号以便引用。
#set math.equation(numbering: "(1)")

参考书：Claes Johnson, *Numerical Solution of Partial Differential Equations by the
Finite Element Method*, Cambridge University Press, 1987.

这是 Johnson 书上椭圆型问题有限元方法的摘记，同时也是本项目数值部分的背景：后面
`theory/02-maxwell-2d.typ` 里那个标量 $E^z$ 问题，用的就是这里的变分形式加 Galerkin
离散，只不过把 $nabla u dot nabla v$ 换成了带 $omega^2 epsilon mu$ 项的算子。

== 椭圆型问题的有限元方法引论

=== 1D 模型问题的变分形式

考虑 $(D)$

$ -u''(x) = f(x), quad 0 < x < 1 \
  u(0) = u(1) = 0 $ <eq:pde>

将其改写为变分问题：

+ $(M)$ 将试探函数 $v$ 乘到式 (@eq:pde) 两端，得

  $ -(u'', v) = (f, v) $

  以及

  $ -(u'', v) = -integral u'' v dif x \
    = -u'(1) v(1) + u'(0) v(0) + (u', v') \
    = (u', v') $

  $ => $

  $ (u', v') = (f, v), forall v in V $ <eq:lemma1>

+ $(V)$ 证明式 (@eq:lemma1) 与变分问题等价。定义

  $ F(v) = frac(1, 2) (v', v') - "Re"(f, v) $

  则

  $ F(u + w) = frac(1, 2) (u' + w', u' + w') - (f, u + w) \
    = frac(1, 2) (u', u') + frac(1, 2) (u', w') + frac(1, 2) (w', u') \
    + frac(1, 2) (w', w') - (f, u + w) \
    = F(u) + frac(1, 2) (w', w') + "Re"(u', w') - "Re"(f, w) \
    >= F(u) $ <eq:varproblem>

  其中的 “$>=$” 取等号当且仅当 $w' = 0$ 时成立。

+ $(D) => (V) <=> (M)$

#code-ref[这个推导就是 `waveguide_optimization/main.cpp` 里 `assemble_waveguide()` 在做的事。
  双线性型 $(u', v')$ 的二维推广是 $integral_Omega nabla u dot nabla v dif Omega$，对应单元刚度矩阵里的
  `Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp])`（`dphi` 即 $nabla N_i$）；右端 $(f, v)$ 对应
  `Fe(i) += JxW_face[qp]*FEM_q*phi_face[i][qp]`。区别只在于本项目的"源"不是体积内的 $f$，
  而是端口上的入射模（`FEM_q = FEM_gamma*2*E0`），并且多了一项
  `Ke(i,j) += JxW[qp]*FEM_beta_wg*(phi[i][qp]*phi[j][qp])`——那是 $omega^2 epsilon mu$ 项，
  不是这里的模型问题。]

#fix[
  原文这一条是与 $(V)$ 的对比，导出 LyX 时把不等号丢了，写成"其中的 “$=$” 当且仅当
  $w' = 0$ 时成立"。按上一式 $F(u + w) = F(u) + 1/2 (w', w') >= F(u)$，应当是
  *不等号 $>=$ 取等号*当且仅当 $w' = 0$。
]

=== 分片线性函数

用分片线性函数基构造 $V$ 的子空间 $V_h$。

- Galerkin 方法：$(V_h)$
- Ritz 方法：$(M_h)$

#code-ref[实际代码用的不是分片线性函数，而是 libMesh 的二阶（二次）Lagrange 单元：
  `main()` 里 `eq_systems.get_system("waveguide").add_variable("Ez", SECOND)`。
  `assemble_waveguide()` 开头的几个对象正好对应这里的三个概念：`FEType fe_type =
  dof_map.variable_type(0)` 是单元类型（决定试探/检验空间 $V_h$），`FEBase::build(dim,
  fe_type)` 负责在积分点上算出形函数 `phi` 和梯度 `dphi`（即 $N_i$ 与 $nabla N_i$），
  `QGauss qrule(dim, FIFTH)` 给出积分点与权重 `JxW`；`DofMap::dof_indices()` 把单元上的
  局部形函数编号映射到全局自由度，也就是展开系数 $c_i$ 的下标。]

=== 误差估计

+ 误差

  $ ((u - u_h)', v') = 0, forall v in V_h $

  这就是 Galerkin 正交性（误差在能量内积下与整个 $V_h$ 正交），是能量范数误差估计的出发点。

#code-ref[误差估计里的网格尺寸 $h$ 在代码里不是手动指定的，而是由
  `waveguide_optimization/EM_mesh_flagging.cpp` 的 `flag_elements()` 自动决定：
  单元体积超过 $(lambda\/16)^d$（`max_elem_vol = pow(curr_wl/16, dim)`，
  `curr_wl = 1/(freq*index)`，`index` 是该子域的折射率）就置 `Elem::REFINE`；
  端口上的边另有 $lambda\/20$ 的判据（`main()` 里构造 `EMMeshFlagging` 时传入的
  `1.0/20`，对应 `if(side_elem->volume() > pow(port_size, dim-1))`）。
  `main()` 用 `do { my_refinement.flag_elements_by(my_model); } while
  (my_refinement.refine_elements())` 反复加密，直到没有单元需要再细分。]

=== Poisson 方程的有限元方法

Poisson 方程

$ -Delta u = f quad "in " Omega \
  u = 0 quad "on " Gamma $

对 $(u v, 0)$ 与 $(0, u v)$ 应用 Green 公式

$ integral_Omega (diff v)/(diff x_i) w dif x + integral_Omega v (diff w)/(diff x_i) dif x = integral_Gamma v w n_i dif s, quad i = 1, 2 $

$ => $

$ integral_Omega nabla v dot nabla w dif x = integral_Gamma v (diff w)/(diff n) dif s - integral_Omega v Delta w dif x $

==== 变分问题：

+ 定义

  $ a(u, v) = integral_Omega nabla u dot nabla v dif x \
    (f, v) = integral_Omega f v dif x \
    V = { v : v "连续于" Omega, nabla v "分片连续有界，且" v = 0 " on " Gamma } $

  #fix[
    原稿中集合 $V$ 的定义写到 "$nabla v$ is piecewise ..." 就中断了（LyX 导出时丢掉了后半句，
    疑为原文未写完），这里按 Johnson 书里 $H_0^1(Omega)$ 的提法补全为"$nabla v$ 分片连续有界
    且 $v = 0$ 在边界 $Gamma$ 上"。
  ]

+ $(M)$ 求 $u in V$，使得

  $ a(u, v) = (f, v), forall v in V $

+ $(V)$ 极小化

  $ F(v) = frac(1, 2) a(v, v) - "Re"(f, v) $
