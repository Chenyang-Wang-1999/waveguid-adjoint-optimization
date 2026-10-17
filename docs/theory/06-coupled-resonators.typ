#import "../template.typ": project-doc, fig, fig-row, note, fix, code-ref

#show: project-doc.with(
  title: "耦合微环：耦合常数与绝热消除",
  date: "2022-01 ~ 2022-03",
)

微环阵列里，环与环之间隔着一小段连接波导（三环情形里是三角形连接环），耦合实际上是
通过它的模式发生的。要用 Haldane 型的紧束缚有效模型算这个阵列的能带与边缘态，
需要两类输入：

- *数值*：外耦合速率 $kappa_"ex"$、内禀损耗 $kappa_"in"$、环与环之间的等效耦合强度
  $t$（最近邻、次近邻两种几何上分别拟合，记作 $t_1$、$t_2$）以及耦合相位 $phi$；
- *结构*：把连接模式消掉之后，环模式之间剩下的等效耦合长什么样。

本文档的两部分正好对应这两件事。第一部分（@sec-coupling）用 input–output formalism
写出两环结构的 S 参数，给出从 S 参数*拟合*上面这些常数的解析公式；第二部分
（@sec-adiabatic）把连接模式*绝热消除*掉，得到只含环模式的等效哈密顿量 $H_"eff"$。
@sec-issues 记下引入反向散射项时那段走不通的推导，@sec-code 说明这些常数在代码里怎么被用起来。

== input–output formalism <sec-formalism>

系统（环）用算符 $a$ 描述，波导里传播的场用连续算符 $b(omega)$、$c$ 描述，两者的耦合写成
Langevin 方程

$ dot(b)(omega) = -i omega b(omega) + kappa(omega) c \
  dot(a) = -(i)/(planck.reduce) [a, H_"sys"] + integral d omega kappa(omega) { b^dagger(omega) [a, c] - [a, c^dagger] b(omega) } $

在 First Markov 近似下耦合系数与频率无关，

$ kappa(omega) = sqrt(gamma/(2pi)) $

并定义输入/输出场

$ b_"in"(t) := frac(1, sqrt(2pi)) integral d omega e^(-i omega (t - t_0)) b_0(omega) \
  b_"out"(t) := frac(1, sqrt(2pi)) integral d omega e^(-i omega (t - t_1)) b_1(omega) $

于是有 input–output 关系

$ b_"out"(t) = b_"in"(t) + sqrt(gamma) c(t) $

代回 Langevin 方程，运动方程就写成只含输入场的形式：

$ dot(a) = -(i)/(planck.reduce) [a, H_"sys"] - { [a, c^dagger] [frac(gamma, 2) c + sqrt(gamma) b_"in"(t)] - [frac(gamma, 2) c^dagger + sqrt(gamma) b_"in"^dagger(t)] [a, c] } $

多通道（多个耦合点、多条波导）时写成

$ dot(a) = -(i)/(planck.reduce) [a, H_"sys"] - sum_(i,j) { [a, c^((i) dagger)] [frac(gamma_(i j), 2) c^((i)) + sqrt(gamma_(i j)) b_"in"^((j))(t)] - [frac(gamma_(i j), 2) c^((i) dagger) + sqrt(gamma_(i j)) b_"in"^((j) dagger)(t)] [a, c^((i))] } $

#fig("assets/figures/waveguide_modes.png", caption: [波导模式的定义：$A_"in"$、$A_"out"$ 画在波导远端，$b_"in"$、$b_"out"$ 画在环与波导的耦合点，两者相差半段波导的传播相位], width: 62%)

后面一直用这套记号：

$ b_"in"^((i)) = A_"in"^((i)) e^(i beta L/2), quad b_"out"^((i)) = A_"out"^((i)) e^(-i beta L/2) $

其中 $A_"in"$、$A_"out"$ 是波导远端的模式幅值，$L$ 是耦合点到端口的长度。实验上（以及有限元里）
能直接得到的是端口之间的 $S$ 参数 $S = A_"out" slash A_"in"$，所以它里面*含有*传播相位
$e^(-i beta L)$；下面所有 $S$ 参数的表达式都把这个因子显式地写在左边。

== 一环两波导

一个环夹在两条波导之间，环与每条波导的耦合速率都是 $kappa_"ex"$（一个环有两个耦合点），
$kappa_"in"$ 是内禀损耗。设 $H_"sys" = planck.reduce omega_0 a^dagger a$，两个耦合点处
被消掉的系统算符是

$ c^((1)) = a, quad c^((2)) = a e^(i beta_"ring" pi r) $

（第二个耦合点看到的是走过半圈之后的场，$beta_"ring"$ 是环的传播常数，$r$ 是环半径。）
把 $b_"in"^((2)) = 0$ 代入上一节的运动方程，稳态解与两个端口的 S 参数是

$ a = frac(sqrt(kappa_"ex"), i Delta - frac(2kappa_"ex" + kappa_"in", 2)) b_"in"^((1)) $

$ S_21 e^(-i beta L) = 1 + frac(kappa_"ex", i Delta - frac(2kappa_"ex" + kappa_"in", 2)) \
  S_41 e^(-i beta L) e^(-i beta_"ring" pi r) = frac(kappa_"ex", i Delta - frac(2kappa_"ex" + kappa_"in", 2)) $

其中 $Delta = omega - omega_0$ 是相对环谐振频率的失谐。$S_21$ 是直通端，$S_41$ 是下话端；
下话端这一式的左端同时出现两个传播因子，是因为信号既要走过半条波导，也要走过半个环。

#fig("assets/figures/one_ring.png", caption: [一个环两波导：两个耦合点处的环模式分别是 $a$ 与 $a e^(i beta_"ring" pi r)$], width: 40%)

#fix[
  原稿在这一节里先算 $b_"out"^((2))$、再换成 $A_"out"^((2))$ 与 $S_41$，其中
  $b_"out"^((2)) = sqrt(kappa_"ex") a e^(+i beta_"ring" pi r)$ 的*下一个等号*丢掉了相位因子
  $e^(+i beta_"ring" pi r)$：原稿写作 $b_"out"^((2)) = frac(kappa_"ex", i Delta - dots) b_"in"^((1))$。
  按原文自己的约定 $c^((2)) = a e^(+i beta_"ring" pi r)$，这个因子不应该丢；而且原稿紧接着的
  $S_41$ 一式左端出现的 $e^(-i beta_"ring" pi r)$ 正是它的逆因子，可见原稿在同一页里前后并不一致。
  这里按 $c^((2))$ 的约定补回 $e^(+i beta_"ring" pi r)$。
]

#note[
  这一节与下一节的推导结构相同（同一组 Langevin 方程、同样的 $c^((i))$ 取法），区别只在于
  两个耦合点开在同一个环上，所以分母里出现的是 $2 kappa_"ex"$；下一节每个环只有一个耦合点，
  出现的是 $kappa_"ex"$。原稿的逐行推导与上面给出的结果等价，整理时压缩成这几行。
]

== 两环耦合：S 参数 <sec-coupling>

#fig("assets/figures/two_rings.png", caption: [两环结构的模式约定：两个环之间直接耦合，耦合强度 $t_2$、相位 $phi$；两侧波导给出直通端（$S_21$）与下话端（$S_31$）], width: 62%) <fig-two-rings>

设两个环的模式 $a^((1))$、$a^((2))$ 之间直接耦合，哈密顿量取

$ H_"sys" = planck.reduce sum_i omega_0 a^(dagger (i)) a^((i)) + t_2 e^(i phi) a^(dagger (2)) a^((1)) + t_2 e^(-i phi) a^(dagger (1)) a^((2)) $

（原稿把这一对环之间的等效耦合记作 $t_2 e^(i phi)$。同一组公式用在最近邻与次近邻两种几何上，
拟合出来的耦合强度分别记作 $t_1$、$t_2$；单看这一节，$t_2$ 就是这两个环之间的耦合。）

两个耦合点处被消掉的算符仍取

$ c^((1)) = a^((1)), quad c^((2)) = a^((2)) e^(i beta_"ring" pi r), quad gamma_(i j) = kappa_"ex" delta_(i j) $

由 $[a^((i)), a^((j) dagger) a^((k))] = delta_(i j) a^((k))$ 得

$ [a^((1)), H_"sys"] = planck.reduce (omega_0 a^((1)) + t_2 e^(-i phi) a^((2))), quad [a^((2)), H_"sys"] = planck.reduce (omega_0 a^((2)) + t_2 e^(i phi) a^((1))) $

代入 Langevin 方程：

$ dot(a)^((1)) = -i(omega_0 a^((1)) + t_2 e^(-i phi) a^((2))) - frac(kappa_"ex" + kappa_"in", 2) a^((1)) - sqrt(kappa_"ex") b_"in"^((1))(t) \
  dot(a)^((2)) = -i(omega_0 a^((2)) + t_2 e^(i phi) a^((1))) - frac(kappa_"ex" + kappa_"in", 2) a^((2)) - e^(-i beta_"ring" pi r) sqrt(kappa_"ex") b_"in"^((2))(t) $

转到旋转坐标系 $tilde(a)^((i)) = a^((i)) e^(-i omega t)$，稳态下 $dot(tilde(a))^((i)) = 0$，
并记

$ kappa_"eff" := frac(kappa_"ex" + kappa_"in", 2) $

得到线性方程组

$ sqrt(kappa_"ex") mat(tilde(b)_"in"^((1))(t); e^(-i beta_"ring" pi r) tilde(b)_"in"^((2))(t)) = mat(i Delta - kappa_"eff", -i t_2 e^(-i phi); -i t_2 e^(i phi), i Delta - kappa_"eff") mat(tilde(a)^((1)); tilde(a)^((2))) $

只从第一个端口输入（$b_"in"^((2)) = 0$），解出

$ mat(tilde(a)^((1)); tilde(a)^((2))) = frac(sqrt(kappa_"ex") tilde(b)_"in"^((1)), t_2^2 + (i Delta - kappa_"eff")^2) mat(i Delta - kappa_"eff"; i t_2 e^(i phi)) $

再用 output 关系 $tilde(b)_"out"^((1)) = tilde(b)_"in"^((1)) + sqrt(kappa_"ex") tilde(a)^((1))$、
$tilde(b)_"out"^((2)) = sqrt(kappa_"ex") tilde(a)^((2)) e^(i beta_"ring" pi r)$，就得到两个端口的 S 参数：

$ S_21 e^(-i beta L) = 1 + frac(kappa_"ex" (i Delta - kappa_"eff"), t_2^2 + (i Delta - kappa_"eff")^2) \
  S_31 e^(-i beta L) e^(-i beta_"ring" pi r) = frac(kappa_"ex" i t_2 e^(i phi), t_2^2 + (i Delta - kappa_"eff")^2) $

#note[
  端口编号沿用原稿：一环两波导那一节下话端记作 $S_41$，这一节（两环）记作 $S_31$。
  另外原稿把 $kappa_"in"$ 时而写成 $kappa_"in"^((1))$、时而写成 $kappa_"in"^((i))$，
  这里统一写成 $kappa_"in"$（两个环的内禀损耗相同）。
]

== 从 S 参数提取耦合常数

=== 幅值

把上面两式取模方：

$ abs(S_21)^2 = frac((t_2^2 + kappa_"eff"^2 - Delta^2 - kappa_"ex" kappa_"eff")^2 + kappa_"in"^2 Delta^2, (t_2^2 + kappa_"eff"^2 - Delta^2)^2 + 4 Delta^2 kappa_"eff"^2) \
  abs(S_31)^2 = frac(kappa_"ex"^2 t_2^2, (t_2^2 + kappa_"eff"^2 - Delta^2)^2 + 4 Delta^2 kappa_"eff"^2) $

记 $u := Delta^2 - t_2^2 - kappa_"eff"^2$（即 $Delta^2 = u + t_2^2 + kappa_"eff"^2$），把
$abs(S)^2$ 都写成 $u$ 的函数：

$ abs(S_21)^2 = frac((u + kappa_"ex" kappa_"eff")^2 + kappa_"in"^2 (u + t_2^2 + kappa_"eff"^2), u^2 + 4(u + t_2^2 + kappa_"eff"^2) kappa_"eff"^2) \
  abs(S_31)^2 = frac(kappa_"ex"^2 t_2^2, u^2 + 4(u + t_2^2 + kappa_"eff"^2) kappa_"eff"^2) $

这两式就是拟合用的公式：把有限元扫频算出来的 $abs(S)^2$ 随 $Delta$（或随频率）的曲线拿来，
对 $kappa_"ex"$、$kappa_"in"$、$t_2$ 做非线性最小二乘即可。

=== 峰值与零点：闭式解

为了给拟合一个好的初值（也为了交叉检验），可以只盯住 $abs(S_31)^2$ 的峰值和零失谐点。
先定义两个无量纲量

$ alpha := frac(kappa_"in", kappa_"ex"), quad zeta := frac(t_2, kappa_"ex") $

*peak*：$abs(S_31)^2$ 的分母在 $u_"max" = -2 kappa_"eff"^2$ 处最小，也就是

$ Delta^2 = t_2^2 - kappa_"eff"^2 = kappa_"ex"^2 (zeta^2 - ((1 + alpha)/2)^2) $

记 $T_31^"max"$ 为 $abs(S_31)^2$ 在峰值处的值，则

$ T_31^"max" = frac(kappa_"ex"^2 t_2^2, 4 kappa_"eff"^2 t_2^2) = (frac(kappa_"ex", kappa_"ex" + kappa_"in"))^2 = frac(1, (1 + alpha)^2) $

于是

$ alpha = frac(1, abs(S_31)_"max") - 1 = frac(1, sqrt(T_31^"max")) - 1 $

*zero*：在零失谐点 $Delta = 0$ 处，$u = -(t_2^2 + kappa_"eff"^2)$，于是

$ T_31 (0) := abs(S_31)^2 \|_(Delta = 0) = frac(kappa_"ex"^2 t_2^2, (t_2^2 + kappa_"eff"^2)^2) $

取倒数开方：

$ frac(1, sqrt(T_31 (0))) = frac(t_2^2 + kappa_"eff"^2, kappa_"ex" t_2) = zeta + ((1 + alpha)/2)^2 frac(1, zeta) $

乘上 $zeta$ 得到关于 $zeta$ 的二次方程及其解：

$ zeta^2 - frac(1, sqrt(T_31 (0))) zeta + ((1 + alpha)/2)^2 = 0, quad "判别式" = frac(1, T_31 (0)) - (1 + alpha)^2 $

$ zeta = frac(1 plus.minus sqrt(1 - T_31 (0) (1 + alpha)^2), 2 sqrt(T_31 (0))) $

两个根里只有 $+$ 号根满足峰值存在（$Delta$ 为实数）的条件

$ zeta^2 = frac(2 - T_31 (0) (1 + alpha)^2 plus.minus 2 sqrt(1 - T_31 (0) (1 + alpha)^2), 4 T_31 (0)) > ((1 + alpha)/2)^2 $

所以取 $+$ 号根。$alpha$ 与 $zeta$ 定下来之后，$kappa_"in" = alpha kappa_"ex"$、
$t_2 = zeta kappa_"ex"$，再由峰值位置 $Delta^2 = kappa_"ex"^2 (zeta^2 - ((1 + alpha)/2)^2)$
就能定出 $kappa_"ex"$。

#fix[
  $alpha$ 与 $zeta$ 的定义在原稿里出现在第一次使用它们*之后*（$alpha$ 写在
  $Delta^2 = kappa_"ex"^2 (zeta^2 - dots)$ 的下一行，$zeta$ 写在 $= zeta + dots$ 的下一行）。
  这里把两个定义提到本节开头。
]

#fix[
  原稿的 $T_31$ 在 peak 与 zero 两处指的是两个不同的数（峰值处的透射率与零失谐处的透射率），
  却用了同一个符号。这里分别记作 $T_31^"max"$ 与 $T_31 (0)$。
]

#fix[
  原稿在 zero 一节里有两处漏了平方：

  - 原稿写作 $frac(1, sqrt(T_31)) = frac(t_2^2 + kappa_"eff"^2, kappa_"ex" t_2) = zeta + frac(1 + alpha, 2) frac(1, zeta)$，
    但由上一式 $kappa_"eff"^2 slash (kappa_"ex" t_2) = ((1 + alpha)/2)^2 slash zeta$，
    第二项应为 $((1 + alpha)/2)^2 slash zeta$；
  - 同理，二次方程的常数项原稿写作 $frac(1 + alpha, 2)$，应为 $((1 + alpha)/2)^2$。

  只有补上平方，判别式才是 $1 slash T_31 (0) - (1 + alpha)^2$（即原稿下一行给出的结果），
  后面的 $zeta$ 解才和它自洽。
]

#fix[
  原稿把 $zeta^2$ 的分母写成了 $2 sqrt(T_31)$；由量纲（$zeta$ 的分母是 $2 sqrt(T_31)$，平方后分母为 $4 T_31$）应为
  $4 T_31$。这里改为 $4 T_31 (0)$。
]

=== 相位

相位不必单独拟合，直接对复数 S 参数做拟合即可。

原因是 $phi$ 只以 $e^(i phi)$ 的形式出现在 $S_31$ 里，$abs(S_31)^2$ 与 $phi$ 无关；
单独拟合相位要先在每个频点上把相位解缠绕，还得先把 $e^(-i beta L)$、
$e^(-i beta_"ring" pi r)$ 这两个传播相位去掉。直接拟合复数 $S$ 则同时约束幅值与相位，
$phi$ 和其它三个常数一次拟合出来。

== 连接模式的绝热消除 <sec-adiabatic>

=== 两个环与一段连接波导

#fig("assets/adiabatic/Two_rings.png", caption: [两个环与一段连接波导：$a_1, a_2$ 是环模式，$c_L, c_H$ 是连接波导的两个模式，$J_L, J_H$ 是连接模式与环的耦合，$phi_("ring")$、$phi_("link")$ 是环内与连接段上的传播相位], width: 58%)

- $a_1, a_2$：回音壁模式（whispering gallery modes）；
- $c_L, c_H$：连接波导模式，$delta_i = omega_i - omega_0$ 是它与环谐振频率的失谐；
- $phi_("ring")$：环模式的相位变化，$phi_("link")$：连接模式的相位变化。

哈密顿量取

$ hat(H) = planck.reduce lr((omega_0 a_1^dagger a_1 + omega_0 a_2^dagger a_2 + sum_(i=L,H) omega_i c_i^dagger c_i)) + planck.reduce sum_(i=L,H) lr((J_i c_i^dagger a_2 e^(i phi_("ring")) + J_i c_i^dagger a_1 e^(-i phi_("link")^((i))) + "H.c.")) $

Heisenberg 运动方程（$i dif a slash dif t = [a, hat(H)] slash planck.reduce$）给出

$ i (dif)/(dif t) a_1 (t) = omega_0 a_1 + sum_(i=L,H) J_i c_i e^(i phi_("link")^((i))) \
  i (dif)/(dif t) a_2 (t) = omega_0 a_2 + sum_(i=L,H) J_i c_i e^(-i phi_("ring")) \
  i (dif)/(dif t) c_i (t) = omega_i c_i + J_i a_2 e^(i phi_("ring")) + J_i a_1 e^(-i phi_("link")^((i))) $

移到以 $omega_0$ 旋转的坐标系（$a, c arrow.r tilde(a) e^(-i omega_0 t), tilde(c) e^(-i omega_0 t)$）后，

$ i (dif)/(dif t) tilde(a)_1 = sum_(i=L,H) J_i tilde(c)_i e^(i phi_("link")^((i))) \
  i (dif)/(dif t) tilde(a)_2 = sum_(i=L,H) J_i tilde(c)_i e^(-i phi_("ring")) \
  i (dif)/(dif t) tilde(c)_i (t) = delta_i tilde(c)_i + J_i tilde(a)_2 e^(i phi_("ring")) + J_i tilde(a)_1 e^(-i phi_("link")^((i))) $

*绝热消除*：连接模式的失谐 $delta_i$ 比耦合 $J_i$ 大得多，它的演化比环模式快得多，
于是可以取 $dif tilde(c)_i slash dif t = 0$，解出

$ tilde(c)_i = -(1)/(delta_i) lr((J_i tilde(a)_2 e^(i phi_("ring")) + J_i tilde(a)_1 e^(-i phi_("link")^((i))))) $

代回 $a_1$、$a_2$ 的运动方程：

$ i (dif)/(dif t) tilde(a)_1 = -sum_(i=L,H) (J_i^2)/(delta_i) lr((tilde(a)_2 e^(i phi_("ring")) e^(i phi_("link")^((i))) + tilde(a)_1)) \
  i (dif)/(dif t) tilde(a)_2 = -sum_(i=L,H) (J_i^2)/(delta_i) lr((tilde(a)_2 + tilde(a)_1 e^(-i phi_("link")^((i))) e^(-i phi_("ring")))) $

写成矩阵形式

$ i (dif)/(dif t) mat(tilde(a)_1; tilde(a)_2) = -sum_(i=L,H) (J_i^2)/(delta_i) mat(1, e^(i phi_("ring")) e^(i phi_("link")^((i))); e^(-i phi_("link")^((i))) e^(-i phi_("ring")), 1) mat(tilde(a)_1; tilde(a)_2) $

在较小频率范围内可以取 $J_L = J_H = J$，此时等效哈密顿量是

$ H_("eff") = -J^2 mat(lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) )), (1)/(abs(delta_H)) e^(i phi_("ring")) e^(i phi_("link")^((H))) - (1)/(abs(delta_L)) e^(i phi_("ring")) e^(i phi_("link")^((L))); (1)/(abs(delta_H)) e^(-i phi_("link")^((H))) e^(-i phi_("ring")) - (1)/(abs(delta_L)) e^(-i phi_("link")^((L))) e^(-i phi_("ring")), lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) ))) $

对连接波导，对称性保证 $phi_("link")$ 只能取 $pi$ 或 $0$；再用回音壁模式的条件
$phi_("link")^((H)) = phi_("link")^((L)) + pi$（即 $e^(i phi_("link")^((H))) = -e^(i phi_("link")^((L)))$），
非对角元里两项同号相加，矩阵化简为

$ H_("eff") = -J^2 mat(lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) )), lr(( (1)/(abs(delta_H)) + (1)/(abs(delta_L)) )) e^(i phi_("ring")) e^(i phi_("link")^((L))); lr(( (1)/(abs(delta_L)) + (1)/(abs(delta_H)) )) e^(-i phi_("link")^((L))) e^(-i phi_("ring")), lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) ))) $

其中 $e^(i phi_("ring")) = plus.minus 1$、$e^(i phi_("link")^((H))) = plus.minus 1$。
这就是把连接波导消掉之后，两个环之间的等效耦合：对角元是各自频率的移动，非对角元是耦合项。

#fix[
  原稿在化简后的非对角元里写的是 $e^(i phi_("link")^((H)))$。由上面的通式，非对角元为
  $-J^2 [(1 slash delta_H) e^(i phi_("link")^((H))) + (1 slash delta_L) e^(i phi_("link")^((L)))] e^(i phi_("ring"))$；
  按 $H$、$L$ 分别是频率较高、较低的连接模式（$delta_H > 0 > delta_L$，即
  $1 slash delta_H = 1 slash abs(delta_H)$、$1 slash delta_L = -1 slash abs(delta_L)$，这与对角元的写法一致），
  再用 $e^(i phi_("link")^((H))) = -e^(i phi_("link")^((L)))$，结果是
  $+ (1 slash abs(delta_H) + 1 slash abs(delta_L)) e^(i phi_("ring")) e^(i phi_("link")^((L)))$。
  这里把 $e^(i phi_("link")^((H)))$ 改回 $e^(i phi_("link")^((L)))$：两者差一个负号，不是同一个矩阵。
]

#note[
  绝热消除的条件，推导稿写的是 $dif tilde(c)_i slash dif t = 0$，周报写的是
  $dif tilde(c)_i slash dif t approx 0$；这里按推导稿保留等号。
]

周报里紧随其后给出的是两组两环结构的 S 参数（振幅与相位），用来核对上面
$phi_("ring")$、$phi_("link")$ 的取值：

#fig-row(
  ("assets/figures/S_amplitude_t1.png", "assets/figures/S_phase_t1.png"),
  caption: [非对称（non-symmetric）情形：左为振幅，右为相位],
  width: 46%,
)

#fig-row(
  ("assets/figures/S_amplitude_t1_sym.png", "assets/figures/S_phase_t1_sym.png"),
  caption: [对称（symmetric）情形：左为振幅，右为相位],
  width: 46%,
)

=== 三个环与三角形连接环

#fig("assets/adiabatic/Three_rings.png", caption: [三个环与三角形连接环：$a_k$（$k = 1, 2, 3$）是环模式，$c_L, c_H$ 是连接环的模式，箭头标出两个传播方向上的相位 $e^(plus.minus i phi)$], width: 58%)

哈密顿量与上一节同构，只是环从一个变成三个，连接模式与每个环的耦合都记作 $J$：

$ H = planck.reduce lr((sum_(k=1)^3 omega_0 a_k^dagger a_k + sum_(k=L,H) omega_k c_k^dagger c_k + sum_(k=1)^3 sum_(j=L,H) J a_k^dagger e^(-i phi_("ring")^((k))) c_j e^(-i (k-2) phi_("link")^((j))) + "H.c")) $

同样移到旋转坐标系（$delta_k = omega_k - omega_0$），运动方程为

$ i (dif)/(dif t) tilde(a)_k = sum_(j=L,H) J tilde(c)_j e^(-i phi_("ring")^((k))) e^(-i (k-2) phi_("link")^((j))) \
  i (dif)/(dif t) tilde(c)_j = delta_j tilde(c)_j + sum_(k=1)^3 J tilde(a)_k e^(+i phi_("ring")^((k))) e^(+i (k-2) phi_("link")^((j))) $

令 $dif tilde(c)_j slash dif t = 0$ 解出 $tilde(c)_j$ 再代回，就得到只含环模式的等效哈密顿量：

$ (H_("eff"))_(k l) = -sum_(j=L,H) (J^2)/(delta_j) e^(+i phi_("ring")^((l)) - i phi_("ring")^((k))) e^(+i (l-k) phi_("link")^((j))) $

三角形连接环给出

$ phi_("link")^((H)) = phi_("link")^((L)) + (2 pi)/3 $

于是括号里两项可以合并成

$ (H_("eff"))_(k l) = -J^2 e^(+i phi_("ring")^((l)) - i phi_("ring")^((k))) e^(+i (l-k) phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(+i (l-k) frac(2 pi, 3)) - (1)/(abs(delta_L)) )) $

环模式的相位取值（由三个环在三角形上的位置决定，见 @fig-phase）为

$ phi_("ring")^((1)) = (5)/(6) pi, quad phi_("ring")^((2)) = -pi/2, quad phi_("ring")^((3)) = pi/6 $

因此相邻环的相位差固定为

$ phi_("ring")^((l+1)) - phi_("ring")^((l)) = (2 pi)/3 $

#fig("assets/adiabatic/phase.png", caption: [环模式相位差的记账：相邻两个环之间是 $plus.minus 2 pi slash 3$（图中 $150 degree + 90 degree = 240 degree equiv -120 degree$ 等）], width: 52%) <fig-phase>

把上面的相位代进通式，就得到三环的等效哈密顿量

$ H_("eff") = -J^2 mat(
  lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) )),
  e^(+i frac(2 pi, 3)) e^(+i phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(+i frac(2 pi, 3)) - (1)/(abs(delta_L)) )),
  e^(-i frac(2 pi, 3)) e^(+i 2 phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(-i frac(2 pi, 3)) - (1)/(abs(delta_L)) ));
  e^(-i frac(2 pi, 3)) e^(-i phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(-i frac(2 pi, 3)) - (1)/(abs(delta_L)) )),
  lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) )),
  e^(+i frac(2 pi, 3)) e^(+i phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(+i frac(2 pi, 3)) - (1)/(abs(delta_L)) ));
  e^(+i frac(2 pi, 3)) e^(-i 2 phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(+i frac(2 pi, 3)) - (1)/(abs(delta_L)) )),
  e^(-i frac(2 pi, 3)) e^(-i phi_("link")^((L))) lr(( (1)/(abs(delta_H)) e^(-i frac(2 pi, 3)) - (1)/(abs(delta_L)) )),
  lr(( (1)/(abs(delta_H)) - (1)/(abs(delta_L)) ))
) $

其中 $phi_("link")^((L)) = 0, plus.minus 2 pi slash 3$。周报里把同一个矩阵以图片形式给出：

#fig("assets/figures/Heff.svg", caption: [周报里以图片形式给出的三环 $H_("eff")$], width: 100%)

#fix[
  通式里非对角元的连接相位是 $e^(+i (l-k) phi_("link")^((L)))$。原稿在矩阵的
  $(1,3)$、$(3,1)$ 两个元素（也就是 $l - k = plus.minus 2$ 的那一对）里把它写成了
  $e^(-i phi_("link")^((L)))$ 与 $e^(+i phi_("link")^((L)))$，指数上少了因子 $l - k = 2$，符号也反了；
  正确的写法是 $(1,3)$ 取 $e^(+i 2 phi_("link")^((L)))$、$(3,1)$ 取 $e^(-i 2 phi_("link")^((L)))$
  （其余四个非对角元 $l - k = plus.minus 1$，原稿写对了）。这里已按通式改正，改过的矩阵仍然 Hermite。
]

#note[
  上面那张周报图片里的矩阵与推导稿的矩阵是同一个，同样带着 $(1,3)$、$(3,1)$ 元素的连接相位问题；
  图片前面写的整体因子是 $+J^2$，而正文推导给的是 $-J^2$，两者相差一个整体符号。
]

== 已知问题：引入反向散射项之后的推导 <sec-issues>

上面两节里的环模式没有区分传播方向。原稿最后一节想把两个传播方向
（赝自旋 $sigma = plus.minus$）通过反向散射耦合起来，即在哈密顿量里加上一项

$ H_"back" = sum_(j=L,H) t c_j^{(+) dagger} c_j^{(-)} + "H.c." $

然后再做一次绝热消除。这一节推导到一半，作者自己写了 "The derivation above may be not correct."；
核对之后我们确认它目前不自洽，所以这里只记下问题，不采用它的结论（后面的能带与边缘态计算
用的都是 @sec-adiabatic 里不含反向散射的 $H_("eff")$）。至少有三处具体的矛盾：

+ 反向散射系数在原文里有三种互不等价的写法：哈密顿量里是常数 $t$，运动方程里变成
  $t^(plus.minus)$（并补充定义 $t^+ = t$、$t^- = t^*$），绝热消除之后又变成 $r u_j$；
  而 $u$ 先被定义成六项指数之和（不带下标），随后又以 $u_j$（$j = L, H$）的形式带下标使用。
  被消掉的到底是什么项，无法确定。
+ 绝热消除后的运动方程指标对不上：写成
  $i dif/(dif t) (a_l^((+)); a_l^((-))) = (H_(++), H_(+-); H_(-+), H_(--)) (a_l^((+)); a_l^((-)))$ 时，
  左端只有下标 $l$，右端却仍对 $k$ 求和；而且 $H_(++)$ 的 $(1,2)$ 元素里把算符 $a_k^((+))$
  留在了矩阵元中，同一个矩阵的其它元素却都是数。
+ $H_(+-)$ 的每个矩阵元里多乘了一个 $J^2$：原稿的系数写成
  $J^2 t sum_j (J^2)/(delta_j^2 - abs(t)^2)$，而同一节的 $H_(++)$ 用的系数是
  $J^2 delta_j slash (delta_j^2 - abs(t)^2)$。按推导 $H_(+-)$ 的系数应当是
  $J^2 t sum_j (1)/(delta_j^2 - abs(t)^2)$；原稿多出来的那个 $J^2$ 使它与 $H_(++)$
  连量纲都对不上。

此外，这一节末了的假设 $phi_("link")^((H)) = phi_("link")^((L)) + pi slash 3$ 与上一节的
$phi_("link")^((H)) = phi_("link")^((L)) + 2 pi slash 3$ 差一个因子 $2$；$H_(++)$ 与上一节的三环
$H_("eff")$ 是同一个矩阵，因此也带着 $(1,3)$、$(3,1)$ 元素连接相位的那处问题。

#fig("assets/adiabatic/backscattering.svg", caption: [反向散射的几何：两个传播方向（$sigma = plus.minus$，红/蓝）的环模式与连接模式，连接段上的相位为 $e^(plus.minus i phi_("link"))$], width: 76%)

== 这些常数在代码里怎么用 <sec-code>

#code-ref[
  这套推导的目的，是从有限元算出来的 S 参数里拟合出 $kappa_"ex"$、$kappa_"in"$、
  $t_1$、$t_2$、$phi$，供 Haldane 有效模型使用。流程是：先用有限元（当年这批耦合结构
  的扫描还有一部分是 COMSOL 做的，见 `docs/assets/comsol/`）扫出一组几何/波长下的
  S 参数，再在 Python 侧用 `scipy` 的非线性最小二乘按上面 $abs(S_21)^2$、$abs(S_31)^2$
  与复数 $S$ 的表达式拟合。拟合前要先把传播相位去掉，这需要传播常数 $beta$：
  `python/propagation_constant0.py` 里的 `get_beta_from_k0(n0, n1, h, k0_sweep)`
  用牛顿迭代解色散方程给出 $beta(k_0)$，这正是 $S_21 e^(-i beta L)$、
  $S_31 e^(-i beta L) e^(-i beta_"ring" pi r)$ 这两个记号的来源（$beta$ 与 $k pi r$
  的数值见当年的参数表）。拟合脚本与扫描数据当时留在 Python 侧，没有进这个仓库
  （原稿 `adiabatic_elimination.lyx` 开头写的是 "see scripts"，指的就是这些脚本）；
  拟合出来的数值结果见 `docs/results/02-ring-array.typ`。
]

#code-ref[
  有效哈密顿量 $H_("eff")$ 把连接模式绝热消除掉之后，只留下环模式之间的等效耦合：
  对角元 $lr((1 slash abs(delta_H) - 1 slash abs(delta_L)))$ 是频率移动，
  非对角元 $lr((1 slash abs(delta_H) + 1 slash abs(delta_L))) e^(i phi_("ring")) e^(i phi_("link"))$
  就是环与环之间的等效耦合项（它的模与相位分别对应拟合出来的 $t$ 和 $phi$）。
  这正是后面用 Haldane 有效模型算能带与边缘态的出发点，见
  `docs/theory/07-bands-and-edge-states.typ`；三个环那一节的 $H_("eff")$ 给出的就是
  三角形上三个环之间的等效耦合矩阵。
]
