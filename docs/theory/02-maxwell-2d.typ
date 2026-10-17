#import "../template.typ": project-doc, note, fix, code-ref

#show: project-doc.with(
  title: "二维麦克斯韦方程与波导的有限元求解",
  date: "2022",
)

== 教程：Poisson 方程

=== Poisson 方程的弱形式

带 Dirichlet 条件的 Poisson 方程

$ -Delta u|_Omega = f \
  u|_(diff Omega) = phi.alt $

首先，第一个方程可以改写为

$ -nabla dot nabla u = f $

$ arrow.r.double $

$ integral_Omega f v^* dif omega = integral_Omega -v^* nabla dot nabla u dif omega $

注意到

$ nabla dot (v^* nabla u) = nabla v^* dot nabla u + v^* nabla^2 u $

$ arrow.r.double $

$ integral_Omega f v^* dif omega = integral_Omega -v^* nabla dot nabla u dif omega \
  = integral_Omega nabla v^* dot nabla u dif omega - integral_(diff Omega) v^* bold(n) dot nabla u dif omega $

$ integral_(diff Omega) v^* u dif sigma = integral_(diff Omega) v^* phi.alt dif sigma $

FEM：设

$ u = sum_i c_i phi_i \
  v = phi_j $

则

$ integral_Omega phi_j^* f dif omega = sum_i c_i integral_Omega nabla phi_j^* dot nabla phi_i dif omega - sum_i c_i integral_(diff Omega) phi_j^* bold(n) dot nabla phi_i dif omega $

$ c_i integral_(diff Omega) phi_j^* phi_i dif sigma = integral_(diff Omega) phi_j^* phi.alt dif sigma $

其中，当 $j$ 是内点时

$ F_e^((j)) = integral_Omega phi_j^* f dif omega \
  K_e^((j,i)) = integral_Omega nabla phi_j^* dot nabla phi_i dif omega $

当 $j$ 是边界点时

$ F_e^((j)) = integral_Omega phi_j^* f dif omega + integral_(diff Omega) phi_j^* phi.alt dif sigma \
  K_e^((j,i)) = integral_Omega nabla phi_j^* dot nabla phi_i dif omega + integral_(diff Omega) phi_j^* phi_i dif sigma - integral_(diff Omega) phi_j^* bold(n) dot nabla phi_i dif omega $

=== 算例

$ nabla^2 u = 0, quad (x,y) in [-1,1] times [-1,1] \
  u|_(diff Omega) = cases(1 & y = 1, 0 & "others") $

==== 精确解

设

$ u = sum_(n=1)^infinity c_n (y) sin((n pi (x+1))/2) $

那么由 $nabla^2 u = 0$ 可得

$ sum_(n=1)^infinity c_n'' (y) sin((n pi (x+1))/2) - sum_(n=1)^infinity (n pi/2)^2 c_n (y) sin((n pi (x+1))/2) = 0 $

这就要求

$ c_n'' (y) = (n pi/2)^2 c_n (y) $

$ arrow.r.double $

$ c_n (y) = C_(1,n) sinh(n pi/2 (y+1)) + C_(2,n) cosh(n pi/2 (y+1)) $

由边界条件

$ sum_(n=1)^infinity C_(2,n) sin((n pi (x+1))/2) = 0 \
  sum_(n=1)^infinity (C_(1,n) sinh(n pi) + C_(2,n) cosh(n pi)) sin((n pi (x+1))/2) = 1 $

$ arrow.r.double $

$ C_(2,n) = 0 \
  sum_(n=1)^infinity C_(1,n) sinh(n pi) sin((n pi (x+1))/2) = 1 $

正交性

$ integral_(-1)^1 sin((n pi (x+1))/2) sin((m pi (x+1))/2) dif x \
  = 1/2 integral_(-1)^1 (cos(((n-m) pi (x+1))/2) - cos(((n+m) pi (x+1))/2)) dif x \
  = delta_(m n) $

因此

$ D_n = integral_(-1)^1 1 * sin((n pi (x+1))/2) dif x \
  = 2/(n pi) cos((n pi (x+1))/2) |_1^(-1) \
  = 2/(n pi) (1 - (-1)^n) $

与前面的展开式比较

$ 1 = sum_(n=1)^infinity C_(1,n) sinh(n pi) sin((n pi (x+1))/2) \
  = sum_(n=1)^infinity D_n sin((n pi (x+1))/2) $

$ arrow.r.double $

$ C_(1,n) = D_n/(sinh n pi) = 2/(n pi) (1 - (-1)^n)/(sinh n pi) $

$ u = sum_(n=1)^infinity c_n (y) sin((n pi (x+1))/2) \
  = sum_(n=1)^infinity (2(1 - (-1)^n))/(n pi) (sinh(n pi/2 (y+1)))/(sinh n pi) sin((n pi (x+1))/2) $

== 单位制

取

$ epsilon_0 = mu_0 = 1 $

于是所有单位都具有长度的量纲。

一些单位的推导：

+ 时间

  $ 1/[c]^2 = [epsilon_0][mu_0] = [T]^2/[L]^2 $

  $ arrow.r.double $

  $ [T]^2 = [epsilon_0][mu_0][L]^2 $

+ 电荷与质量

  $ [F] = [q]^2/([epsilon_0][L]^2) = ([M][L])/[T]^2 = [M]/([epsilon_0][mu_0][L]) $

  $ arrow.r.double $

  $ [q]^2 = ([M][L])/[mu_0] $

  质量

  $ [S] = [M]/[T]^3 \
    [M] = [S][T]^3 = [S]([epsilon_0][mu_0])^(3/2)[L]^3 \
    [I]_("current") = [q]/[T] = ([M]/([epsilon_0][mu_0]^2[L]))^(1/2) $

+ $H, E, B, D$

  $ [H] = [q]/([L][T]) = ([S][epsilon_0]^(1/2)/[mu_0]^(1/2))^(1/2) = ([S]/[Z_0])^(1/2) \
    [E] = ([S][Z_0])^(1/2) \
    [B] = [mu_0][H] = (([mu_0][S])/[c])^(1/2) \
    [D] = [epsilon_0][E] = (([epsilon_0])/[c][S])^(1/2) $

== 二维麦克斯韦方程

=== 二维麦克斯韦方程组

（$z$ 方向均匀）

$ nabla dot bold(D) = rho \
  nabla dot bold(B) = 0 \
  nabla times bold(E) = -(diff)/(diff t) bold(B) = i omega bold(B) \
  nabla times bold(H) = bold(J) + (diff)/(diff t) bold(D) = bold(J) - i omega bold(D) $

其中 $(diff)/(diff z) = i beta_z$。

$ diff_i (epsilon_(i j) E^j) + i beta_z epsilon_z E^z = rho \
  diff_i (mu_(i j) H^j) + i beta_z mu_z H^z = 0 $

$ i beta_z E^2 = diff_2 E^z - i omega mu_(1 j) H^j \
  i beta_z E^1 = diff_1 E^z + i omega mu_(2 j) H^j \
  i beta_z H^2 = i omega epsilon_(1 j) E^j + diff_2 H^z - J^1 \
  i beta_z H^1 = -i omega epsilon_(2 j) E^j + diff_1 H^z + J^2 $

$ arrow.r.double $

$ i beta_z E^i = diff_i E^z + i omega (delta_1^i delta_2^j - delta_2^i delta_1^j) mu_(j k) H^k \
  i beta_z H^i = diff_i H^z + (delta_1^i delta_2^j - delta_2^i delta_1^j)(J^j - i omega epsilon_(j k) E^k) $

$ diff_1 E^2 - diff_2 E^1 = i omega mu_(33) H^z \
  diff_1 H^2 - diff_2 H^1 = J^z - i omega epsilon_(33) E^z $

$ i beta_z bold(E)_t = nabla_t E^z + i omega mat(0, 1; -1, 0) bold(mu) bold(H)_t \
  i beta_z bold(H)_t = nabla_t H^z + mat(0, 1; -1, 0)(bold(J)_t - i omega bold(epsilon) bold(E)_t) \
  nabla_t dot mat(0, 1; -1, 0) bold(E)_t = i omega mu_z H^z \
  nabla_t dot mat(0, 1; -1, 0) bold(H)_t = J^z - i omega epsilon_z E^z \
  nabla_t dot (bold(epsilon) bold(E)_t) + i beta_z epsilon_z E^z = rho \
  nabla_t dot (bold(mu) bold(H)_t) + i beta_z mu_z H^z = 0 $

==== TE mode（$E_z = 0$）

$ beta_z bold(E)_t = omega mat(0, 1; -1, 0) bold(mu) bold(H)_t \
  i beta_z bold(H)_t = nabla_t H^z + mat(0, 1; -1, 0)(bold(J)_t - i omega bold(epsilon) bold(E)_t) \
  nabla_t dot mat(0, 1; -1, 0) bold(E)_t = i omega mu_z H^z \
  nabla_t dot mat(0, 1; -1, 0) bold(H)_t = J^z \
  nabla_t dot (bold(epsilon) bold(E)_t) = rho \
  nabla_t dot (bold(mu) bold(H)_t) + i beta_z mu_z H^z = 0 $

+ $beta_z = 0$

  $ 0 = omega mat(0, 1; -1, 0) bold(mu) bold(H)_t \
    0 = nabla_t H^z + mat(0, 1; -1, 0)(bold(J)_t - i omega bold(epsilon) bold(E)_t) \
    nabla_t dot mat(0, 1; -1, 0) bold(E)_t = i omega mu_z H^z \
    nabla_t dot mat(0, 1; -1, 0) bold(H)_t = J^z \
    nabla_t dot (bold(epsilon) bold(E)_t) = rho \
    nabla_t dot (bold(mu) bold(H)_t) = 0 $

  $ arrow.r.double $

  $ 0 = bold(H)_t \
    i omega mat(0, 1; -1, 0) bold(E)_t = bold(epsilon)^(-1)(nabla_t H^z + mat(0, 1; -1, 0) bold(J)_t) \
    nabla_t dot (bold(epsilon)^(-1) nabla_t H^z) + omega^2 mu_z H^z = -nabla_t dot bold(epsilon)^(-1) mat(0, 1; -1, 0) bold(J)_t \
    0 = J^z \
    nabla_t dot (bold(epsilon) mat(0, -1; 1, 0) bold(epsilon)^(-1)(nabla_t H^z + mat(0, 1; -1, 0) bold(J)_t)) = i omega rho $

==== TM mode（$H_z = 0$）

$ i beta_z bold(E)_t = nabla_t E^z + i omega mat(0, 1; -1, 0) bold(mu) bold(H)_t \
  i beta_z bold(H)_t = mat(0, 1; -1, 0)(bold(J)_t - i omega bold(epsilon) bold(E)_t) \
  nabla_t dot mat(0, 1; -1, 0) bold(E)_t = 0 \
  nabla_t dot mat(0, 1; -1, 0) bold(H)_t = J^z - i omega epsilon_z E^z \
  nabla_t dot (bold(epsilon) bold(E)_t) + i beta_z epsilon_z E^z = rho \
  nabla_t dot (bold(mu) bold(H)_t) = 0 $

+ $beta_z = 0$

  $ i omega mat(0, 1; -1, 0) bold(H)_t = -bold(mu)^(-1) nabla_t E^z \
    i omega bold(epsilon) bold(E)_t = bold(J)_t \
    nabla_t dot mat(0, 1; -1, 0) bold(E)_t = 0 \
    nabla_t dot (bold(mu)^(-1) nabla_t E^z) + omega^2 epsilon_z E^z = -i omega J^z \
    nabla_t dot (bold(mu) bold(H)_t) = 0 $

#code-ref[本项目只解 TM 这一支（$H_z = 0$、未知量是标量 $E^z$）：`main()` 里
  `eq_systems.get_system("waveguide").add_variable("Ez", SECOND)`，在
  `assemble_waveguide()` 中按 $integral (nabla E^z dot nabla v - omega^2 epsilon mu E^z v) dif Omega$
  装配——这正是上面 $nabla_t dot (bold(mu)^(-1) nabla_t E^z) + omega^2 epsilon_z E^z = -i omega J^z$
  乘以检验函数分部积分的结果（$mu = 1$）：
  `Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);` 与
  `Ke(i,j) += JxW[qp]*FEM_beta_wg*(phi[i][qp]*phi[j][qp]);`，
  其中 `FEM_beta_wg = -pow(omega,2)*epsilon_wg*mu_wg`。
  注意代码里 `FEM_beta_*` 这个名字指的是弱形式里 $omega^2 epsilon mu$ 那个系数（带负号），
  与下文作为传播常数的 $beta$ 无关。TE 支（$E^z = 0$、未知量是 $H^z$）在代码里没有实现。]

=== 边值问题

==== 问题

问题

$ (diff)/(diff x)(alpha_x (diff phi)/(diff x)) + (diff)/(diff y)(alpha_y (diff phi)/(diff y)) + beta phi = f \
  phi = p, quad "on" Gamma_1 \
  (alpha_x (diff phi)/(diff x) hat(x) + alpha_y (diff phi)/(diff y) hat(y)) dot hat(n) + gamma phi = q, quad "on" Gamma_2 $

不连续条件

$ phi^+ = phi^- quad "on" Gamma_d \
  (alpha_x^+ (diff phi^+)/(diff x) hat(x) + alpha_y^+ (diff phi^+)/(diff y) hat(y)) dot hat(n) = (alpha_x^- (diff phi^-)/(diff x) hat(x) + alpha_y^- (diff phi^-)/(diff y) hat(y)) dot hat(n) $

==== 变分形式

$ cases(delta F(phi) = 0, phi = p & "on" Gamma_1) $

其中

$ F(phi) = 1/2 integral.double_Omega lr([alpha_x ((diff phi)/(diff x))^2 + alpha_y ((diff phi)/(diff y))^2 + beta phi^2]) dif Omega \
  + integral_(Gamma_2) (gamma/2 phi^2 - q phi) dif Gamma - integral.double_Omega f phi dif Omega $

对于不连续的情形，把每一部分按连续的 $alpha$ 拆开后，同样可以证明该变分形式。

==== von Neumann 边界

$ phi^e (x,y) = sum_(j=1)^(N_p) N_j^e (x,y) phi_j^e \
  F(phi) = sum_(e=1)^M F^e (phi^e) $

则

$ (diff F^e)/(diff phi_i^e) = sum_(j=1)^(N_p) phi_j^e integral.double_(Omega^e) (alpha_x (diff N_i^e)/(diff x) (diff N_j^e)/(diff x) + alpha_y (diff N_i^e)/(diff y) (diff N_j^e)/(diff y) + beta N_i^e N_j^e) dif Omega \
  - integral.double_(Omega^e) f N_i^e dif Omega $

==== 第三类边界条件

$ F_b (phi) = integral_(Gamma_2) (gamma/2 phi^2 - q phi) dif Gamma \
  F_b (phi) = sum_(s=1)^(M_s) F_b^s (phi^s) \
  phi^s = sum_(j=1)^(N_s) N_j^s phi_j^s $

则

$ F_b (phi) = integral_(Gamma_2) sum_i (gamma/2 sum_j N_j^s phi_j^s N_i^s phi_i^s - q phi_i^s N_i^s) dif Gamma \
  (diff F_b)/(diff phi_i^s) = integral_(Gamma_2) (gamma/2 sum_j N_j^s phi_j^s N_i^s - q phi_i^s N_i^s) dif Gamma \
  = sum_j phi_j^s integral_(Gamma_2) gamma N_j^s N_i^s dif Gamma - integral_(Gamma_2) q N_i^s dif Gamma $

==== Dirichlet 边界条件

Dirichlet 条件 $phi = p$ 在有限元里用惩罚法（penalty method）处理：把约束以惩罚因子 $P$
（一个大数，代码注释里取 $P = 10^20$）乘进泛函，得到的单元贡献就是下面"有限元形式"
一节里的 $K_(i j) = P integral_(Gamma_1) N_j^s N_i^s dif Gamma$ 和
$F_i = P integral_(Gamma_1) N_i^s p dif Gamma$。$P$ 越大约束越硬，代价是矩阵条件数变差。

=== 金属波导

$ nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z mu_(x y) E^z = 0 \
  E^z = 0, quad "side boundaries" \
  E^z = E_0 e^(i k_0 x) + R E_0 e^(-i k_0 x), quad "input" \
  E^z = T E_0 e^(i k_0 x), quad "output" $

==== 本征模 $E_0$ 的计算

设

$ E^z = E_0 (x) e^(i beta y) $

代入体方程

$ nabla_t dot (e^(i beta y) (dif)/(dif x) E_0 (x) hat(e)_x + E_0 (x) i beta e^(i beta y) hat(e)_y) + omega^2 epsilon_z mu E_0 (x) e^(i beta y) = 0 $

$ ((dif)/(dif x))^2 E_0 (x) = -(omega^2/c^2 epsilon_z mu_(x y) - beta^2) E_0 (x) $

解为

$ E_0^((m)) (x) = sin((m pi)/w x) \
  ((dif)/(dif x))^2 E_0^((m)) (x) = -((m pi)/w)^2 E_0^((m)) (x) $

因此

$ beta^2 = omega^2/c^2 epsilon_z mu_(x y) - ((m pi)/w)^2 $

==== 端口边界条件

+ 输入边界

  $ E^z = E_0 e^(i beta x) + R E_0 e^(-i beta x) \
    (diff E^z)/(diff x) = i beta E_0 e^(i beta x) - i beta R E_0 e^(-i beta x) $

  $ arrow.r.double $

  $ E^z (0) = E_0 + R E_0 \
    (diff E^z)/(diff x)(0) = i beta E_0 - i beta R E_0 $

  $ arrow.r.double $

  $ R = (E^z (0) - E_0)/E_0 \
    (diff E^z)/(diff x)(0) = i beta E_0 - i beta R E_0 = 2 i beta E_0 - i beta E^z (0) $

  在端口上 $E^z (0) = E^z$，把它代回上式，并把 $diff E^z slash diff x (0)$ 写成
  端口上的 $diff E^z slash diff x$：

  $ (diff E^z)/(diff x) + i beta E^z = 2 i beta E_0 $

  输入端口在 $x = 0$、求解域在 $x > 0$ 一侧，外法向是 $-hat(x)$，即
  $diff slash diff n = - diff slash diff x$，两端乘 $-1$ 得

  $ (diff E^z)/(diff n) - i beta E^z = -2 i beta E_0 $

  #fix[
    原稿写作
    $ R = (E^z (0) - E_0)/E_0 \
      (diff E^z)/(diff x) = i beta E_0 - i beta E^z (0) + i beta E_0 \
      (diff E^z)/(diff x) + i beta E^z = 2 i beta E_0 $
    中间一步把 $R E_0 = E^z (0) - E_0$ 代进去之后，左边是*在 $x = 0$ 处*取值的
    $diff E^z slash diff x (0)$，右边却留着 $E^z (0)$，最后一行又直接用 $E^z$ 替换掉了它——
    也就是漏写了 $E^z (0) = E^z$（在端口上端口值就是场值）这一步；从
    $diff E^z slash diff x$ 到 $diff E^z slash diff n$ 的外法向符号也省略了。这里补上这两步，
    结论不变。
  ]

+ 输出边界

  $ E^z = T E_0 e^(i beta x) \
    (diff E^z)/(diff x) = i beta T E_0 e^(i beta x) = i beta E^z $

  $ arrow.r.double $

  $ (diff E^z)/(diff n) - i beta E^z = 0 \
    T = E^z/E_0 $

==== 有限元形式

+ 原型方程

  $ (diff)/(diff x)(alpha_x (diff phi)/(diff x)) + (diff)/(diff y)(alpha_y (diff phi)/(diff y)) + beta phi = f \
    phi = p, quad "on" Gamma_1 \
    (alpha_x (diff phi)/(diff x) hat(x) + alpha_y (diff phi)/(diff y) hat(y)) dot hat(n) + gamma phi = q, quad "on" Gamma_2 $

  我的方程

  $ nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z mu_(x y) E^z = 0 \
    E^z = 0, quad "side boundaries" \
    (diff E^z)/(diff n) - i beta E^z = -2 i beta E_0, quad "input" \
    (diff E^z)/(diff n) - i beta E^z = 0, quad "output" $

  其中

  $ alpha = I \
    f = 0 \
    beta = omega^2/c^2 epsilon_z mu_(x y) \
    p = 0 \
    gamma = -i beta \
    q = -2 i beta E_0, quad "input" $

+ Bulk

  $ (diff F^e)/(diff phi_i^e) = sum_(j=1)^(N_p) phi_j^e integral.double_(Omega^e) (alpha_x (diff N_i^e)/(diff x) (diff N_j^e)/(diff x) + alpha_y (diff N_i^e)/(diff y) (diff N_j^e)/(diff y) + beta N_i^e N_j^e) dif Omega $

+ Third（第三类边界）

  $ (diff F_b)/(diff phi_i^s) = integral_(Gamma_2) (gamma/2 sum_j N_j^s phi_j^s N_i^s - q phi_i^s N_i^s) dif Gamma \
    = sum_j phi_j^s integral_(Gamma_2) gamma N_j^s N_i^s dif Gamma - integral_(Gamma_2) q N_i^s dif Gamma $

+ Dirichlet

  $ K_(i j) = P integral_(Gamma_1) N_j^s N_i^s dif Gamma \
    F_i = P integral_(Gamma_1) N_i^s p dif Gamma $

#code-ref[这些边界项就是 `waveguide_optimization/main.cpp` 的 `assemble_waveguide()` 里
  遍历单元各条边（`for (auto side : elem->side_index_range())`，只处理
  `elem->neighbor_ptr(side) == nullptr` 的真边界）时装配的东西：

  - 端口（第三类边界）：`FEM_gamma = - 1j*beta[port_index];` 对应 $gamma = -i beta$，
    矩阵部分是 `Ke(i,j) += JxW_face[qp]*FEM_gamma*phi_face[i][qp]*phi_face[j][qp]`，
    即 $integral_(Gamma_2) gamma N_i N_j dif Gamma$；只有 `port_incidence_on(port_index)`
    为真的端口才加右端 `FEM_q = FEM_gamma*2*E0;` 与
    `Fe(i) += JxW_face[qp]*FEM_q*phi_face[i][qp]`，对应
    $q = -2 i beta E_0$ 的 $integral_(Gamma_2) q N_i dif Gamma$。端口上的 $beta$ 与
    $E_0$ 来自 `WaveguideParams::analytic_E_field()`（见下文平板波导一节）。
  - 非端口边界：`FEM_gamma = - Complex(1j)*sqrt(epsilon_amb*mu_amb)*omega;`，即
    $gamma = -i k$、$q = 0$，这就是下面"散射边界条件"的一阶透明边界。
  - Dirichlet/PEC 那一支在代码里是注释掉的（`const Real penalty = 1.e20;` 所在的那段），
    当前模型的四个外边界全是端口，用不到。
]

==== 散射边界条件

$ bold(n) times (nabla times bold(E)) + i k bold(n) times (bold(E) times bold(n)) = 0 $

$ epsilon_(i j k) epsilon_(k l m) (n^j diff_l E^m + i k n^j E^l n^m) = 0 \
  (delta_(i l) delta_(j m) - delta_(i m) delta_(j l))(n^j diff_l E^m + i k n^j E^l n^m) = 0 $

+ $E$ 偏振

  $ (delta_(i l) delta_(j 3) - delta_(i 3) delta_(j l)) n^j diff_l E^z + (delta_(i 3) delta_(j m) - delta_(i m) delta_(j 3)) i k n^j E^z n^m = 0 \
    -delta_(i 3) delta_(j l) n^j diff_l E^z + delta_(i 3) delta_(j m) i k n^j E^z n^m = 0 \
    -n dot nabla E^z + i k E^z = 0 $

==== 后处理

+ 由 $E$ 求 $H$

  $ nabla times bold(E) = -(diff)/(diff t) bold(B) = i omega bold(B) $

  $ arrow.r.double $

  $ i omega bold(mu) bold(H) = nabla times bold(E) \
    = mat(diff_y E^z - diff_z E^y; diff_z E^x - diff_x E^z; diff_x E^y - diff_y E^x) $

  $ arrow.r.double $

  $ bold(H) = -i/omega bold(mu)^(-1) mat(diff_y E^z - diff_z E^y; diff_z E^x - diff_x E^z; diff_x E^y - diff_y E^x) $

+ Poynting vector

  $ angle.l S angle.r = 1/2 "Re"(bold(E)^* times bold(H)) $

  $ arrow.r.double $

  $ S_x = 1/2 "Re"(E^(y*) H^z - E^(z*) H^y) \
    S_y = 1/2 "Re"(E^(z*) H^x - E^(x*) H^z) \
    S_z = 1/2 "Re"(E^(x*) H^y - E^(y*) H^x) $

#code-ref[后处理在代码里有两个入口，公式与上面完全一致。`main.cpp` 的 `mode_decompose()`（S 参数
  的模式分解）在端口积分点上做

  ```cpp
  Number Hx_qp = - (Complex(1j)/(omega*mu)) * Ez_grad(1);   // Hx = -i/(omega*mu) * dEz/dy
  Number Hy_qp =   (Complex(1j)/(omega*mu)) * Ez_grad(0);   // Hy = +i/(omega*mu) * dEz/dx
  result += (std::conj(mode_Dz_qp)*Ez_qp + std::conj(mode_Bx_qp)*Hx_qp
             + std::conj(mode_By_qp)*Hy_qp)*JxW[qp]/4;
  ```

  其中 `Ez_grad = context.side_gradient(0, qp)` 就是 $nabla_t E^z$；模式场的 $D_z, B_x, B_y$
  由 `WaveguideParams::analytic_E_field()` 给出（`B_x = mu H_x`、`B_y = mu H_y`，见该函数的
  `case WaveguideModes::Bx / By`），所以被积函数正是
  $1/4 (bold(D)_"mode"^* dot bold(E) + bold(H)_"mode"^* dot bold(B))$。
  `N_poav_integrand()` 是同一套 $H_x, H_y$ 用在另一种调用形式里：它算
  `Poavx = -(conj(Ez)*Hy).real()/2`、`Poavy = (conj(Ez)*Hx).real()/2`，再投影到（平均）法向
  得到 $angle.l S angle.r dot hat(n)$，即上面 Poynting 矢量的法向分量。]

=== 平板波导

==== $E_z$ 偏振

$ cases(nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z^("wg") mu_(x y)^("wg") E^z = 0 & "waveguide", nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z^("amb1") mu_(x y)^("amb1") E^z = 0 & "ambient 1", nabla_t dot (nabla_t E^z) + omega^2/c^2 epsilon_z^("amb2") mu_(x y)^("amb2") E^z = 0 & "ambient 2") $

$ H^y = i/(omega mu) diff_x E^z = -beta/(omega mu) E^z \
  H^x = -i/(omega mu) diff_y E^z $

令 $n_("wg")^2 = epsilon_z^("wg") mu_(x y)^("wg")$，$n_("amb1")^2 = epsilon_z^("amb1") mu_(x y)^("amb1")$，$n_("amb2")^2 = epsilon_z^("amb2") mu_(x y)^("amb2")$，则

$ ((diff)/(diff y))^2 E^z = (beta^2 - ((n_i omega)/c)^2) E^z, quad "domain" i $

设

$ E^z = cases(A_+ exp(-kappa_+ y) & y > w/2, C_1 sin(k y) + C_2 cos(k y) & abs(y) < w/2, A_- exp(kappa_- y) & y < -w/2) $

则

$ (dif E^z)/(dif y) = cases(-kappa_+ A_+ exp(-kappa_+ y) & y > w/2, k C_1 cos(k y) - k C_2 sin(k y) & abs(y) < w/2, kappa_- A_- exp(kappa_- y) & y < -w/2) $

连续性条件

$ A_+ exp(-(kappa_+ w)/2) - C_1 sin((k w)/2) - C_2 cos((k w)/2) = 0 \
  A_- exp(-(kappa_- w)/2) + C_1 sin((k w)/2) - C_2 cos((k w)/2) = 0 \
  -kappa_+ A_+ exp(-(kappa_+ w)/2) - k C_1 cos((k w)/2) + k C_2 sin((k w)/2) = 0 \
  kappa_- A_- exp(-(kappa_- w)/2) - k C_1 cos((k w)/2) - k C_2 sin((k w)/2) = 0 $

$ exp(-(kappa_+ w)/2) exp(-(kappa_- w)/2) mat(delim: "|", 1, 0, -sin((k w)/2), -cos((k w)/2); 0, 1, +sin((k w)/2), -cos((k w)/2); -kappa_+, 0, -k cos((k w)/2), +k sin((k w)/2); 0, kappa_-, -k cos((k w)/2), -k sin((k w)/2)) = 0 $

$ arrow.r.double $

$ -cos((k w)/2)^2 kappa_+ k - cos((k w)/2)^2 kappa_- k - 2 cos((k w)/2) kappa_+ kappa_- sin((k w)/2) \
  + 2 cos((k w)/2) sin((k w)/2) k^2 + kappa_+ sin((k w)/2)^2 k + kappa_- sin((k w)/2)^2 k = 0 $

$ arrow.r.double $

$ -cos^2((k w)/2)(kappa_+ + kappa_-) k + sin^2((k w)/2)(kappa_+ + kappa_-) k \
  + [k^2 - (k sin((k w)/2) - k cos((k w)/2))^2] \
  + (cos((k w)/2) kappa_+ - sin((k w)/2) kappa_-)^2 - cos^2((k w)/2) kappa_+^2 - sin^2((k w)/2) kappa_-^2 = 0 $

#fix[
  原稿这一行写作
  $ -cos^2((k w)/2)(kappa_+ + kappa_-) k + sin^2((k w)/2)(kappa_+ + kappa_-) k + k^2 \
    + (cos((k w)/2) kappa_+ - sin((k w)/2) kappa_-)^2 - cos^2((k w)/2) kappa_+^2 - sin^2((k w)/2) kappa_-^2 \
    - (k sin((k w)/2) - k cos((k w)/2))^2 = 0 $
  其中的 $+k^2$ 孤立地写在式子中间，看上去像是多出来的一项。它其实是后面
  $-(k sin((k w)/2) - k cos((k w)/2))^2$ 的"配平方"搭档：
  $k^2 - (k sin - k cos)^2 = k^2 [1 - (sin - cos)^2] = 2 k^2 sin cos$（用 $sin^2 + cos^2 = 1$），
  正好等于前一步里的 $2 cos((k w)/2) sin((k w)/2) k^2$；而
  $(cos kappa_+ - sin kappa_-)^2 - cos^2 kappa_+^2 - sin^2 kappa_-^2 = -2 cos sin kappa_+ kappa_-$
  给出交叉项。所以原式在代数上是对的，只是把一对量拆开写了。这里把 $k^2$ 与那个平方项括在一起，
  数值完全不变（两条式子都已用符号计算逐项核对过）。
]

对于 $kappa_+ = kappa_-$

$ -2 cos((k w)/2)^2 kappa k - 2 cos((k w)/2) kappa^2 sin((k w)/2) + 2 cos((k w)/2) sin((k w)/2) k^2 + 2 kappa sin((k w)/2)^2 k = 0 \
  (k sin((k w)/2) - kappa cos((k w)/2))(k cos((k w)/2) + kappa sin((k w)/2)) = 0 $

$ arrow.r.double $

$ k = kappa cot((k w)/2) $

或

$ -k = kappa tan((k w)/2) $

$ A_+ = C_1 sin((k w)/2) + C_2 cos((k w)/2) \
  A_- = -C_1 sin((k w)/2) + C_2 cos((k w)/2) \
  -A_+ kappa = k C_1 cos((k w)/2) - k C_2 sin((k w)/2) \
  A_- kappa = k C_1 cos((k w)/2) + k C_2 sin((k w)/2) $

+ $k = kappa cot((k w)/2)$

  $ A_+ = C_1 sin((k w)/2) + C_2 cos((k w)/2) \
    A_- = -C_1 sin((k w)/2) + C_2 cos((k w)/2) \
    -A_+ = C_1 cot((k w)/2) cos((k w)/2) - C_2 cos((k w)/2) \
    A_- = C_1 cot((k w)/2) cos((k w)/2) + C_2 cos((k w)/2) $

  $ arrow.r.double $

  $ A_- = C cos((k w)/2) \
    A_+ = C cos((k w)/2) $

+ $-k = kappa tan((k w)/2)$

  $ A_+ = C_1 sin((k w)/2) + C_2 cos((k w)/2) \
    A_- = -C_1 sin((k w)/2) + C_2 cos((k w)/2) \
    A_+ = C_1 sin((k w)/2) - C_2 tan((k w)/2) sin((k w)/2) \
    -A_- = C_1 sin((k w)/2) + C_2 tan((k w)/2) sin((k w)/2) $

  $ arrow.r.double $

  $ A_+ = C sin((k w)/2) \
    A_- = -C sin((k w)/2) $

+ 归一化

  $ H^y = i/omega bold(mu)^(-1) diff_x E^z = -beta/omega bold(mu)^(-1) E^z \
    S^x = -1/2 E^z H^y = 1/2 beta/(omega mu) E^z E^z \
    w = 1/2 bold(E) dot bold(D) + 1/2 bold(H) dot bold(B) \
    = 1/4 "Re"(bold(E)^* dot bold(D) + bold(H)^* dot bold(B)) \
    = 1/4 (epsilon + beta^2/(omega^2 mu))(E^z)^2 $

  奇模

  $ E^z = cases(C sin((k w)/2) exp(-kappa y) & y > w/2, C sin(k y) & abs(y) < w/2, -C sin((k w)/2) exp(kappa y) & y < -w/2) $

  $ angle.l E, E angle.r = C^2 sin^2((k w)/2) 1/4 (epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) integral_(-infinity)^(-w/2) exp(2 kappa y) dif y \
    + 1/4 (epsilon_("wg") + beta^2/(omega^2 mu_("wg"))) C^2 integral_(-w/2)^(w/2) sin^2(k y) dif y \
    + 1/4 (epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) C^2 sin^2((k w)/2) integral_(w/2)^infinity exp(-2 kappa y) dif y \
    = C^2 1/2 sin^2((k w)/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) integral_(-infinity)^(-w/2) exp(2 kappa y) dif y \
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg"))) integral_(-w/2)^(w/2) (1 - cos(2 k y)) dif y \
    = C^2 1/2 sin^2((k w)/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) 1/(2 kappa) exp(-kappa w) \
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg")))(w - 1/(2k) sin(2 k y)|_(-w/2)^(w/2)) \
    = C^2 1/(4 kappa) sin^2((k w)/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) exp(-kappa w) \
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg")))(w - 1/k sin(k w)) $

  偶模

  $ E^z = cases(C cos((k w)/2) exp(-kappa y) & y > w/2, C cos(k y) & abs(y) < w/2, C cos((k w)/2) exp(kappa y) & y < -w/2) $

  $ angle.l E, E angle.r = C^2 cos^2((k w)/2) 1/4 (epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) integral_(-infinity)^(-w/2) exp(2 kappa y) dif y \
    + 1/4 (epsilon_("wg") + beta^2/(omega^2 mu_("wg"))) C^2 integral_(-w/2)^(w/2) cos^2(k y) dif y \
    + 1/4 (epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) C^2 cos^2((k w)/2) integral_(w/2)^infinity exp(-2 kappa y) dif y \
    = C^2 1/2 cos^2((k w)/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) integral_(-infinity)^(-w/2) exp(2 kappa y) dif y \
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg"))) integral_(-w/2)^(w/2) (1 + cos(2 k y)) dif y \
    = C^2 1/2 cos^2((k w)/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) 1/(2 kappa) exp(-kappa w) \
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg")))(w + 1/(2k) sin(2 k y)|_(-w/2)^(w/2)) \
    = C^2 1/(4 kappa) cos^2((k w)/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) exp(-kappa w) \
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg")))(w + 1/k sin(k w)) $

+ 功率输出

  + 偶模

    $ S_x = 1/2 beta/(omega mu) E^z E^z \
      = 1/2 beta/(omega mu_("amb")) integral_(w/2)^infinity C^2 cos^2((k w)/2) exp(-2 kappa y) dif y \
      + 1/2 beta/(omega mu_("amb")) integral_(-infinity)^(-w/2) C^2 cos^2((k w)/2) exp(2 kappa y) dif y \
      + 1/2 beta/(omega mu_("wg")) integral_(-w/2)^(w/2) C^2 cos^2(k y) dif y \
      = beta/(omega mu_("amb")) cos^2((k w)/2) 1/(2 kappa) C^2 exp(-kappa w) \
      + 1/4 beta/(omega mu_("wg")) integral_(-w/2)^(w/2) C^2 (1 + cos(2 k y)) dif y \
      = beta/(omega mu_("amb")) cos^2((k w)/2) 1/(2 kappa) C^2 exp(-kappa w) \
      + 1/4 beta/(omega mu_("wg")) C^2 (w + 1/k sin(k w)) $

  + 奇模

    $ S_x = 1/2 beta/(omega mu) E^z E^z \
      = 1/2 beta/(omega mu_("amb")) integral_(w/2)^infinity C^2 sin^2((k w)/2) exp(-2 kappa y) dif y \
      + 1/2 beta/(omega mu_("amb")) integral_(-infinity)^(-w/2) C^2 sin^2((k w)/2) exp(2 kappa y) dif y \
      + 1/2 beta/(omega mu_("wg")) integral_(-w/2)^(w/2) C^2 sin^2(k y) dif y \
      = beta/(omega mu_("amb")) sin^2((k w)/2) 1/(2 kappa) C^2 exp(-kappa w) \
      + 1/4 beta/(omega mu_("wg")) integral_(-w/2)^(w/2) C^2 (1 - cos(2 k y)) dif y \
      = beta/(omega mu_("amb")) sin^2((k w)/2) 1/(2 kappa) C^2 exp(-kappa w) \
      + 1/4 beta/(omega mu_("wg")) C^2 (w - 1/k sin(k w)) $

#code-ref[这一节的偶模/奇模场、归一化积分 $angle.l E, E angle.r$ 和功率 $S_x$，就是
  `waveguide_optimization/2D_waveguide_model_params.cpp` 的 `analytic_E_field()` 里逐行写下
  来的东西：

  - `mode_index % 2` 选模支：奇数取 `sin`（奇模），偶数取 `cos`（偶模，基模 `mode_index = 0`）；
    芯层里的 $k$ 与包层里的 $kappa$ 分别是
    `k = sqrt(omega*omega*epsilon_wg*mu_wg - beta*beta)`、
    `kappa = sqrt(beta*beta - omega*omega*epsilon_amb*mu_amb)`，正是上面的 $k$ 与 $kappa$。
  - 归一化积分：`mode_energy` 的两项逐字对应上面奇模的
    $C^2 1/(4 kappa) sin^2(k w\/2)(epsilon_("amb") + beta^2/(omega^2 mu_("amb"))) exp(-kappa w)
    + C^2 1/8 (epsilon_("wg") + beta^2/(omega^2 mu_("wg")))(w - 1/k sin(k w))$（偶模把
    $sin^2 -> cos^2$、$w - -> w +$）；`normalize = true` 时取 `C_factor = 1/sqrt(mode_energy)`。
  - 功率：`mode_flux` 对应上面的 $S_x$ 积分结果，`normalize = false` 时
    `C_factor = amplitude/sqrt(mode_flux)`（按给定功率定标）。
  - 场分量：芯层 $abs(x) < w\/2$ 用 `Ez = C_factor*cos(k*x)`（偶模）/`sin(k*x)`（奇模），
    包层用 `cos(k*w/2)*exp(-+kappa*x)`；`Hy_loc = - beta/(omega*mu)*Ez` 与上面的
    $H^y = -beta/(omega mu) E^z$ 一致，`Hx_loc` 由 $diff_x E^z$ 给出。
  - 注意代码里的横向坐标叫 `x`（`pt` 沿 `port_x_axis` 投影），本文这一节里叫 $y$，两者是同一个
    东西；`Dz = epsilon*Ez`、`Bx = mu*Hx`、`By = mu*Hy` 由该函数末尾的 `switch` 乘上局部
    $epsilon, mu$ 得到，而 `Hx_loc` / `Hy_loc` 是端口局部坐标下的分量，返回前要按
    `x_axis`（横向轴 `port_x_axis` 旋转 $-90degree$）与 `y_axis`（即 `port_x_axis`）
    投影成全局的 $H_x, H_y$，才能和 `mode_decompose()` 里由全局梯度算出的 $H_x, H_y$ 逐点对上。]

== 基于 Bézier 曲线的拓扑优化

=== 边界的变分

$ F(phi) = 1/2 integral.double_Omega lr([alpha_x ((diff phi)/(diff x))^2 + alpha_y ((diff phi)/(diff y))^2 + beta phi^2]) dif Omega \
  + integral_(Gamma_2) (gamma/2 phi^2 - q phi) dif Gamma - integral.double_Omega f phi dif Omega $

$ delta F(phi) = 1/2 integral.double_("deformed boundary") lr([(alpha_x^((2)) - alpha_x^((1)))((diff phi)/(diff x))^2 + (alpha_y^((2)) - alpha_y^((1)))((diff phi)/(diff y))^2 + (beta^((2)) - beta^((1))) phi^2]) dif Gamma $

沿界面对形状函数积分，就得到 $K$ 矩阵的导数。

伴随法与 Bézier 控制点梯度的完整推导见 `theory/05-adjoint-sensitivity.typ`。
