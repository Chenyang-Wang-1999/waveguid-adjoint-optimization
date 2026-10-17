#import "../template.typ": project-doc, fig, note, fix, code-ref

#show: project-doc.with(
  title: "平板波导传播常数的计算",
  subtitle: "色散方程的牛顿迭代与数值稳定性",
  date: "2022-01",
)

参考书：《导波光学》，ISBN 9787030189226。

这一节解决的问题是：给定三层平板波导（衬底 $n_0$ / 芯层 $n_1$ / 覆盖层 $n_2$、芯层厚度 $h$）
和工作频率 $k_0 = omega\/c$，求出导模的传播常数 $beta$。它是后面所有端口边界条件
（$beta$ 必须出现在一阶透明边界里）和模式归一化的输入。

== 记号与色散方程

设横向衰减/振荡常数为

$ kappa_1 = sqrt(k_0^2 n_1^2 - beta^2) $

$ p_0 = sqrt(beta^2 - k_0^2 n_0^2) $

$ p_2 = sqrt(beta^2 - k_0^2 n_2^2) $

其中 $kappa_1$ 是芯层里的横向波数（导模要求 $n_0 k_0 < beta < n_1 k_0$，因此 $kappa_1$ 和
$p_0, p_2$ 都是实的），$p_0, p_2$ 是两侧包层里的衰减常数。

色散方程（相位匹配条件）为

$ kappa_1 h = m pi + phi_(10) + phi_(12) $

$ phi_(1 i) = arctan(p_i/kappa_1) $

$m$ 是模阶数，$phi_(1 i)$ 是芯层-包层界面上的反射相位。

当两侧包层对称（$n_0 = n_2$，于是 $p_0 = p_2$）时，方程化简为

$ tan(kappa_1 h/2) = tan(phi_(10)) = p_0/kappa_1 $

即

$ kappa_1 tan(kappa_1 h/2) = p_0 $

对称情形的导模有两支，上面这一支是*偶模*（芯层里 $E^z prop cos(kappa_1 y)$）：由
$kappa_1 h = m pi + 2 phi_(10)$ 取 $m = 0$ 即得。另一支是*奇模*
（$E^z prop sin(kappa_1 y)$），满足

$ -kappa_1 cot(kappa_1 h/2) = p_0 $

它对应 $m$ 为奇数；两支合起来就是 $tan(kappa_1 h/2 - m pi/2) = p_0\/kappa_1$。
`theory/02-maxwell-2d.typ` 里从场分布和导数在界面上的连续性条件直接解出了这两支：
$k = kappa cot(k w\/2)$（偶模）与 $-k = kappa tan(k w\/2)$（奇模），其中那里的 $k$ 是芯层
横向波数（即这里的 $kappa_1$）、$kappa$ 是包层衰减常数（即这里的 $p_0$），横向坐标和芯层
厚度则记作 $y$ 和 $w$。下面的牛顿迭代只针对偶模那一支（也是代码里用到的那一支，
模式序号 `mode_index = 0`）。

== 牛顿迭代法

把色散方程写成求根问题 $f(beta) = 0$：

$ f(beta) = kappa_1 tan(kappa_1 h/2) - p_0 $

$ f'(beta) = (diff kappa_1)/(diff beta) (tan(kappa_1 h/2) + kappa_1 h/2 sec^2(kappa_1 h/2)) - (diff p_0)/(diff beta) $

其中

$ (diff kappa_1)/(diff beta) = -beta/kappa_1, quad (diff p_0)/(diff beta) = beta/p_0 $

于是可以套牛顿迭代

$ beta_(n+1) = beta_n - f(beta_n)/f'(beta_n) $

== 数值不稳定与 arctan 形式

直接用上面的 $f(beta)$ 做牛顿迭代会碰到数值不稳定：`tan` 在 $kappa_1 h/2$ 接近
$pi/2$ 的整数倍时会发散，$f$ 的零点附近斜率极大，迭代容易跳出去。

解决办法是把方程写成 `arctan` 形式，把 `tan` 解掉。由
$arctan(tan x) = x - floor(x/pi) pi$，取主值分支后得到

$ g(beta) = kappa_1 h/2 - arctan(p_0/kappa_1) $

$ g'(beta) = -h/2 dot beta/kappa_1 - beta/(p_0 kappa_1) $

`arctan` 是有界的（值域在 $(-pi/2, pi/2)$），因此 $g$ 在整个可行区间上都是缓变的，
牛顿迭代收敛稳定得多。

#note[
  导数 $g'(beta)$ 可以直接对 $g$ 求导验证：
  $g' = (h/2) kappa_1' - [1/(1 + (p_0\/kappa_1)^2)] (p_0\/kappa_1)'$，
  代入 $kappa_1' = -beta\/kappa_1$、$p_0' = beta\/p_0$ 化简后正好得到上式。
]

#code-ref[`python/propagation_constant0.py` 就是这三条公式的直译：

  ```python
  def cost_function_2(beta, h, n0, n1, k0):
      # cost function g(beta) = kappa1*h/2 - arctan(p0/kappa1)
      kappa1 = np.sqrt((k0**2) * (n1**2) - beta**2)
      p0 = np.sqrt(beta**2 - (k0**2)*(n0**2))
      return (kappa1*h/2 - np.arctan(p0/kappa1))
  ```

  `cost_function_2` 返回的正是 $g(beta)$（函数上方的注释写成
  `f(beta) = kappa1*h - arctan(p0/kappa1)`，比函数体少了一个 $1\/2$，以函数体为准），
  `diff_cost_function_2` 返回上式的
  $g'(beta) = -h\/2 dot beta\/kappa_1 - beta\/(p_0 kappa_1)$，逐项一致；
  $kappa_1, p_0$ 分别由 `sqrt(k0**2*n1**2 - beta**2)`、`sqrt(beta**2 - k0**2*n0**2)`
  算出。C++ 侧的 `analytic_E_field()` 用的是同一对量，只是换了名字：
  `k = sqrt(omega*omega*epsilon_wg*mu_wg - beta*beta)`（对应 $kappa_1$）、
  `kappa = sqrt(beta*beta - omega*omega*epsilon_amb*mu_amb)`（对应 $p_0$），
  并以 `mode_index % 2` 在偶模/奇模之间选择。]

== 代码实现

这一段推导在两种语言里都实现过，用的是同一个 $g(beta)$：

- Python：`propagation_constant0.py` 的 `get_beta_from_k0()` 里，
  `cost_function_2` 对应 $g(beta)$、`diff_cost_function_2` 对应 $g'(beta)$，
  用 `scipy.optimize.newton` 迭代，并且把上一次的解作为下一个频率点的初值
  （频率扫描时 $beta(k_0)$ 是连续的，这样可以保证不跳支）：

  ```python
  beta = opt.newton(cost_fun, beta_approx, diff_cost_fun)
  beta_sweep.append(beta)
  beta_approx = beta        # 用上一个频率点的解做初值
  ```

  初值在进入循环前取 $beta = (n_0 + 8 n_1) k_0 slash 9$，落在导模区间
  $(n_0 k_0, n_1 k_0)$ 内部。`python/original_mesh.py` 调用它时只传了当前
  $k_0 = 2 pi f$（`k0_sweep=[freq*2*pi]`），算出 $beta$ 后写进
  `port_other_params` 的第 4 行，C++ 侧的 `WaveguideParams` 再从这个矩阵里读 $beta$。

- C++：`2D_waveguide_model_params.cpp` 里的 `analytic_E_field()`，为端口提供解析的
  $E_0$、$beta$ 以及由此算出的 $H_x$、$H_y$、$B_x$、$B_y$（它不自己解色散方程，
  $beta$ 是 `port_other_params(3, curr_port)` 读进来的）。

#note[
  原文的代码片段里有一处保护性检查 `if(n1*k0 < beta or n0*k0 > beta): print("Error! ...")`，
  用来在迭代跑出导模区间（$n_0 k_0 < beta < n_1 k_0$）时报警。原文这一行是裸的
  `if`，实际上应当包在每次迭代之后。
]
