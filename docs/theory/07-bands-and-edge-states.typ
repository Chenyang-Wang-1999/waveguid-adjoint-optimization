#import "../template.typ": project-doc, fig, fig-row, fig-gallery, note, fix, code-ref

#show: project-doc.with(
  title: "微环阵列的能带与边缘态：模型",
  date: "2022-03 ~ 2022-04",
)

// 摘要：模板的 `abstract:` 参数在 docs/template.typ 里漏了 `#`（第 108、118 行），
// 直接传参会把源码当正文印出来，这里手工排出同样的样式。
#block(
  width: 100%,
  inset: (x: 1em, y: 0.8em),
  fill: luma(96%),
  radius: 3pt,
)[
  #set text(size: 9.5pt)
  *摘要。* 微环阵列是芯片上的“光子量子霍尔效应”平台：环与环之间的耦合系数扮演跃迁振幅，
  连接环带来的绕行相位扮演规范场。本笔记整理 2022 年 3–4 月建立的模型——
  从两个环的 4×4 哈密顿量、次近邻（三角连接）耦合，到 Haldane 有效模型与几何参数，
  最后说明这些系数在本项目的求解器里对应哪些步骤。
  具体的拟合数值、能带与边缘态计算结果见 `docs/results/02-ring-array.typ`。
]
#v(1em)

== 两谐振环的耦合

=== 哈密顿量与 4×4 矩阵

把每个微环的两种模式（记作奇 $o$、偶 $e$；物理上就是环内驻波的两种横向宇称）
当作格点轨道，左侧环的两个模式记作 $A_("Lo")$、$A_("Le")$，右侧环记作 $A_("Ro")$、$A_("Re")$。
两个环的线性耦合可以写成紧束缚形式的哈密顿量

$ H = sum_(m,n in {A_("Lo"), A_("Le"), A_("Ro"), A_("Re")}) h_(m n) a_m^dagger a_n $

其中 $a_m^dagger$、$a_n$ 是模式 $m$、$n$ 的产生、湮灭算符，$h_(m n)$ 是耦合矩阵元。在以
$(A_("Lo"), A_("Le"), A_("Ro"), A_("Re"))$ 为基的表示下，$h$ 是一个实对称的 4×4 矩阵：

$ h mat(A_("Lo"); A_("Le"); A_("Ro"); A_("Re")) = mat(
    Omega, , t_("oo"), t_("oe");
    , Omega, t_("oe"), t_("ee");
    t_("oo"), t_("oe"), Omega, ;
    t_("oe"), t_("ee"), , Omega
  ) mat(A_("Lo"); A_("Le"); A_("Ro"); A_("Re")) $

对角元 $Omega$ 是单个环的谐振频率（两个环相同），非对角元：
$t_("oo")$ 是左右环奇模之间的耦合，$t_("ee")$ 是左右环偶模之间的耦合，
$t_("oe")$ 是奇模与偶模之间的耦合。空白项为 0——同一个环内的两个模式正交，不直接耦合。

#fix[原稿只写了“求解 $h - omega_0 I$ 本征值”。这里补全为本征方程
$(h - omega_0 I) v = omega v$；符号约定与原稿一致，即 $omega$ 从参考频率 $omega_0$ 量起。]

=== 本征值与本征矢

令 $delta Omega = Omega - omega_0$ 为单环谐振频率相对参考频率的失谐，解本征方程得四个本征值

$ omega = 1/2 (plus.minus (t_("ee") + t_("oo")) plus.minus sqrt((t_("ee") - t_("oo"))^2 + 4 t_("oe")^2)) + delta Omega $

其中

$ Delta = sqrt((t_("ee") - t_("oo"))^2 + 4 t_("oe")^2) $

也就是说，把 4×4 矩阵按左右环分块写成 $h = Omega I + mat(0, B; B, 0)$（$B$ 是以 $(o,e)$ 为基的
2×2 对称矩阵）之后，四个本征值是 $plus.minus lambda_+$ 与 $plus.minus lambda_-$，
$lambda_\pm = 1/2 (t_("ee") + t_("oo") plus.minus Delta)$，对应的本征矢是 $(u; plus.minus u)$ 形式。

取第 4 个分量归一化为 1，本征矢矩阵为

$ V = mat(
    -frac(Delta - t_("ee") + t_("oo"), 2 t_("oe")), -frac(Delta + t_("ee") - t_("oo"), 2 t_("oe")), -frac(-Delta - t_("ee") + t_("oo"), 2 t_("oe")), -frac(-Delta + t_("ee") - t_("oo"), 2 t_("oe"));
    -1, 1, -1, 1;
    -frac(-t_("oo") Delta + t_("ee") t_("oo") - 2 t_("oe")^2 - t_("oo")^2, t_("oe") (Delta + t_("ee") + t_("oo"))), -frac(-t_("oo") Delta - t_("ee") t_("oo") + 2 t_("oe")^2 + t_("oo")^2, t_("oe") (Delta - t_("ee") - t_("oo"))), -frac(-t_("oo") Delta - t_("ee") t_("oo") + 2 t_("oe")^2 + t_("oo")^2, t_("oe") (Delta - t_("ee") - t_("oo"))), -frac(-t_("oo") Delta + t_("ee") t_("oo") - 2 t_("oe")^2 - t_("oo")^2, t_("oe") (Delta + t_("ee") + t_("oo")));
    1, 1, 1, 1
  ) $

$V$ 的第 $i$ 列是第 $i$ 个本征值的本征矢（未归一化）。

=== $t_("oe") = 0$ 的退化情形

两环并排、结构关于两环连线的中垂线（以及连线本身）对称时，奇模与偶模不混合，$t_("oe") = 0$。
这时 $Delta = |t_("ee") - t_("oo")|$，本征值退化为

$ omega = {-t_("ee"), t_("ee"), -t_("oo"), t_("oo")} + delta Omega $

本征矢也退化为纯的宇称组合：

$ V = mat(0, 0, -1, 1; -1, 1, 0, 0; 0, 0, 1, 1; 1, 1, 0, 0) $

即四个模式分别具有确定的“左右宇称”与“奇偶宇称”。这正是把环的两种模式当作两个子格
（赝自旋）来用的前提：$t_("ee")$ 与 $t_("oo")$ 给出两套独立的耦合，$t_("oe")$ 是破坏这一图像的微扰。

=== 谐振频率与 FSR

取环波导宽度 $l_("ring") = 0.1 r_("ex")$、波导折射率 $n_("wg") = 3.4$，单环的谐振频率
（以 $c\/r_("ex")$ 为单位，$m$ 为角向量子数）为

$ m = 13 &: 0.92521 \
  m = 14 &: 0.97525 \
  m = 15 &: 1.0248 $

相邻纵模间隔即自由光谱范围

$ "FSR" = 0.050 c/r_("ex") $

与上表相邻两行的差值（$0.05004$、$0.04955$）一致。若把环半径取为 $r_("ex")$，
则 $2 pi r_("ex") n_("eff") nu \/ c = m$，由 $m = 15$ 可反解出 $n_("eff") = m \/ (2 pi nu) approx 2.33$ ——
这是把无量纲频率换算回有效折射率的关系。

=== 两环耦合的结果：$t_("ee")$、$t_("oo")$ 与 $delta Omega$

直接对两环结构做本征频率计算（$m = 15$，$omega_0 = 1.024776 c\/r_("ex")$），
把四个本征频率相对 $omega_0$ 的变化量按模式宇称列出来：

#figure(
  {
    set text(size: 8.5pt)
    table(
      columns: 5,
      align: center,
      stroke: 0.4pt + luma(70%),
      inset: (x: 0.45em, y: 0.32em),
      table.header(
        [*距离（单位 $r_("ex")$）*],
        [*x 偶 y 偶*],
        [*x 偶 y 奇*],
        [*x 奇 y 奇*],
        [*x 奇 y 偶*],
      ),
      [$0$], [$-1.3126726212 times 10^(-2)$], [$-1.3122551529 times 10^(-2)$], [$2.0999619319 times 10^(-2)$], [$2.1192387116 times 10^(-2)$],
      [$0.2$], [$-1.131629 times 10^(-3)$], [$-1.131686 times 10^(-3)$], [$1.297629 times 10^(-3)$], [$1.298016 times 10^(-3)$],
      [$0.4$], [$-1.240000 times 10^(-4)$], [$-1.240000 times 10^(-4)$], [$1.280000 times 10^(-4)$], [$1.280000 times 10^(-4)$],
      [$0.6$], [$-1.580000 times 10^(-5)$], [$-1.620000 times 10^(-5)$], [$1.760000 times 10^(-5)$], [$1.800000 times 10^(-5)$],
    )
  },
  caption: [本征频率变化量 $omega - omega_0$ vs 耦合距离，$m = 15$，$omega_0 = 1.024776 c\/r_("ex")$。],
  kind: table,
  supplement: "表",
  placement: none,
) <tab-eigen-shift>

把上表按 $omega = {-t_("ee"), t_("ee"), -t_("oo"), t_("oo")} + delta Omega$ 反解，
得到耦合系数随距离的变化：

#figure(
  {
    set text(size: 8.5pt)
    table(
      columns: 4,
      align: center,
      stroke: 0.4pt + luma(70%),
      inset: (x: 0.45em, y: 0.32em),
      table.header(
        [*距离（单位 $r_("ex")$）*],
        [*$t_("ee")$*],
        [*$t_("oo")$*],
        [*$delta Omega$*],
      ),
      [$0$], [$-1.71596 times 10^(-2)$], [$-1.70611 times 10^(-2)$], [$0.004$],
      [$0.2$], [$-1.21482 times 10^(-3)$], [$-1.21466 times 10^(-3)$], [$8.3 times 10^(-5)$],
      [$0.4$], [$-1.26000 times 10^(-4)$], [$-1.26000 times 10^(-4)$], [$2 times 10^(-6)$],
      [$0.6$], [$-1.69000 times 10^(-5)$], [$-1.69000 times 10^(-5)$], [$1 times 10^(-6)$],
    )
  },
  caption: [耦合系数，$m = 15$，$omega_0 = 1.024776 c\/r_("ex")$。],
  kind: table,
  supplement: "表",
  placement: none,
) <tab-coupling-coeffs>

两张表描述同一组计算，但只能做量级校核：@tab-coupling-coeffs 给出的
$delta Omega = 0.004$ 与 @tab-eigen-shift 中四个值的中点 $3.99 times 10^(-3)$ 相符。
若进一步把 $(t_("ee"), t_("oo"))$ 代入 $omega = {-t_("ee"), t_("ee"), -t_("oo"), t_("oo")} + delta Omega$
反推，得到的四个值与 @tab-eigen-shift 相差在 $6 times 10^(-5)$ 以内（约为数值本身的 0.4%）。
两组数都做了舍入（$delta Omega$ 只保留一位有效数字），$t_("oe")$ 又未在表中给出，
无法判断残差来自拟合还是来自 $t_("oe")$，因此两组数值一律照原样保留。
耦合强度随距离从 $1.7 times 10^(-2)$（$d = 0$）衰减到 $1.7 times 10^(-5)$（$d = 0.6$），
大致每 $0.2 r_("ex")$ 掉一个数量级——这也是把环间距离当作可调参数来改变带隙的依据之一。

#note[原表（$omega - omega_0$）的列序为“x 偶 y 偶、x 偶 y 奇、x 奇 y 奇、x 奇 y 偶”，
而配套的场分布图注按“左上 x 偶 y 偶、右上 x 偶 y 奇、左下 x 奇 y 偶、右下 x 奇 y 奇”描述，
最后两列的顺序不一致。原稿也没有说明这四个宇称标记与 $(t_("ee"), t_("oo"), t_("oe"))$ 的对应关系；
缺少原始数据无法判定，故表头与数值一律照原样保留。]

== 次近邻耦合：三角连接

最近邻的两个环可以直接耦合，也可以各自通过一段连接波导耦合；*次近邻*的两个环之间
则隔着一个由三段连接波导拼成的“三角”。把这个三角（连同它的两个连接模式）绝热消除掉，
两个环之间就出现一个等效的次近邻跃迁

$ t_2 e^(i phi) a_2^dagger a_1 + "h.c." $

其中振幅 $t_2$ 由连接环与环的耦合 $J$ 和连接环的失谐 $delta_i$ 决定，相位 $phi$ 来自
连接环的绕行相位 $phi_("ring")$ 与连接波导的相位 $phi_("link")$ 之和。完整的推导
（两环 + 一个连接环、三环 + 三角形连接环、以及反向散射项）见
`docs/theory/06-coupled-resonators.typ`；这里只记录该结构下拟合出来的耦合系数。

几何参数：$n = 3.4$，环波导宽度 $l_("ring") = 0.1 r_("ex")$，
三角连接波导宽度 $l_("tri") = 0.01 a = 0.02 sqrt(3) r_("ex")$（$a$ 是环阵列的晶格常数）。
三角顶点耦合（上标 $(2)$ 表示二阶、即次近邻过程）的拟合结果：

$ t_("ee")^((2)) = -4.17 times 10^(-5) \
  t_("oo")^((2)) = -1.796 times 10^(-4) \
  delta Omega_("ee")^((2)) = -0.004189 \
  delta Omega_("oo")^((2)) = -0.0047394 $

与最近邻耦合相比，$t_("oo")^((2))$ 比 $t_("ee")^((2))$ 大 4 倍多——三角连接对奇模的耦合明显更强。
连接三角同时带来一个远大于耦合强度的失谐 $delta Omega^((2)) approx -4.2 times 10^(-3)$
（比 $t^((2))$ 大一个量级以上，符号也与两环直接耦合时的 $delta Omega$ 相反）：
它相当于给单环谐振频率 $Omega$ 加了一个平移，本身不进 Haldane 模型，
但会明显改变环阵列的工作频率。

=== 不同连接组态的能级

为了挑出合适的连接方式，逐一比较了几种组态（把连接环放在三角的顶点或边、单个或双三角、
以及两环/三环/六环的排布）下、本征频率相对 $omega_0$ 的变化量：

#figure(
  {
    set text(size: 8.5pt)
    table(
      columns: 7,
      align: (left, center, center, center, center, center, center),
      stroke: 0.4pt + luma(70%),
      inset: (x: 0.45em, y: 0.32em),
      [三角], [$-0.027212$], [$-0.027396$], [$-0.027212$], [$0.027079$], [], [],
      [双环，顶角], [$-0.0045345$], [$-0.0045017$], [$0.0032438$], [$0.0033359$], [], [],
      [双环，边], [$-3 times 10^(-5)$], [$-3 times 10^(-5)$], [$-3 times 10^(-5)$], [$-3 times 10^(-5)$], [], [],
      [双环，双三角], [$-0.0075662$], [$-0.0058602$], [$-2.2732 times 10^(-4)$], [$-1.6533 times 10^(-4)$], [], [],
      [三环，顶角], [$-0.0061529$], [$-0.0061523$], [$-0.0013098$], [$-0.0013094$], [$0.0065624$], [$0.0067247$],
      [三环，双三角], [$-0.008215$], [$-0.0082147$], [$-0.0030433$], [$-0.003043$], [$0.0019159$], [$0.0039945$],
    )
  },
  caption: [各种连接组态下的本征频率变化量（相对 $omega_0$）。第一列是组态，后面是各模式的数值。],
  kind: table,
  supplement: "表",
  placement: none,
) <tab-config-levels>

#note[原表没有表头，也没有说明各列对应哪个模式；前四行只填了 4 个数，后两行填了 6 个数。
数值照录。表中“三角”一行第 1、3 格完全相同（$-0.027212$），与
“同一组态下不同模式不应严格简并”的预期不符，疑为抄写错误；因无法确定正确值，保持原样。]

从这张表里最终选中了“三角”这一组态（@tab-config-levels 的第一行）：它的本征频率移动最大
（$2.7 times 10^(-2)$ 量级），也就是连接最强；随后按“选择三角”一节取
$d = 0.14 r_("ex")$、$rho = 0.1 a$。

=== 工作点的其他参数

#figure(
  {
    set text(size: 8.5pt)
    table(
      columns: (auto, 1fr),
      align: (left, left),
      stroke: 0.4pt + luma(70%),
      inset: (x: 0.45em, y: 0.32em),
      table.header([*符号*], [*取值 / 含义*]),
      [$n$], [芯层折射率 $3.4$],
      [$l_("ring")$], [环波导宽度 $0.1 r_("ex")$],
      [$l_("tri")$], [三角连接波导宽度 $0.01 a = 0.02 sqrt(3) r_("ex")$],
      [$d_("ring")$], [两个环之间的距离（两环耦合计算中取 $0.80$–$1.10 r_("ex")$；阵列结构中取 $0.1 r_("ex")$）],
      [$d_("tri")$], [三角连接环到环的距离 $0.08 r_("ex")$],
      [$rho_("tri")$], [三角连接环的半径参数 $0.085 a$],
      [$d$, $rho$], [“选择三角”一节的取法：$d = 0.14 r_("ex")$，$rho = 0.1 a$，$Delta nu = -0.029568$ 与 $0.025241$],
      [$nu$, $m$], [工作点 $nu = 1.552169$（单位 $c\/r_("ex")$），角向量子数 $m = 26$],
      [圆环横向高阶模], [约 $1.78 c\/r_("ex")$，工作频率必须选在它以下],
    )
  },
  caption: [微环阵列的几何与工作参数。],
  kind: table,
  supplement: "表",
  placement: none,
) <tab-geometry-params>

#note[$d_("ring")$ 在两处含义不同：两环耦合的拟合里它是两个环之间的距离（$0.80$–$1.10 r_("ex")$，
连接环半径随之取 $r = (d_("ring") - w_("ring")) \/ 2$）；而在 $m = 26$ 的阵列结构里
它被记作 $0.1 r_("ex")$，量级差一个数量级。原稿没有说明两者的关系，此处按出现位置分别列出，
不做统一。]

== Haldane 模型与能带拟合

=== 有效哈密顿量

把环阵列的两种环模式（奇/偶）看成两个独立的赝自旋，每个赝自旋对应一个 Haldane 模型：
格点分成 $A$、$B$ 两个子格（对应环阵列里两套不等价的环位），只保留最近邻跃迁与
带相位的次近邻跃迁，有效哈密顿量为

$ H = sum_(<<i j>>) t_1 a_i^dagger b_j + sum_(<<<i j>>>) t_2 e^(i phi_(i j)) a_i^dagger a_j + t_2 e^(-i phi_(i j)) b_i^dagger b_j + h.c. + delta nu sum_i (a_i^dagger a_i - b_i^dagger b_i) $

其中 $angle.l i j angle.r$ 表示最近邻、$angle.l angle.l i j angle.r angle.r$ 表示次近邻，
$phi_(i j) = plus.minus phi$ 在 $A$、$B$ 两个子格上反号（量子霍尔效应里就是这一项提供了
等效磁场）。写成 $k$ 空间的 2×2 形式：

$ h(k) = mat(d_0(k) + m(k), f(k); overline(f(k)), d_0(k) - m(k)) $

$ d_0(k) = 2 t_2 cos phi sum_i cos(k dot b_i), quad
  m(k) = delta nu - 2 t_2 sin phi sum_i sin(k dot b_i), quad
  f(k) = t_1 sum_i e^(-i k dot a_i) $

这里 $a_i$ 是最近邻矢量、$b_i$ 是次近邻矢量：$k$ 空间对角项里出现 $sin phi$ 的那一项就是
Haldane 的“质量项”，它使对角质量 $m(k)$ 在 $K$、$K'$ 两点取相反的符号，从而给出拓扑非平凡的
带隙——当 $delta nu = plus.minus 3 sqrt(3) t_2 sin phi$ 时带隙在 $K$（或 $K'$）点闭合，
两侧分别是陈数 $cal(C) = plus.minus 1$ 的相。拟合参数里的 $t_1$、$t_2$、$phi$、$delta nu$、
$cal(C)$ 五个量就对应这里的跃迁振幅、相位、失谐项和拓扑符号。

#note[$delta nu$ 的物理身份原稿没有写清楚：它既可能是 Haldane 模型的子格质量项 $M$，
也可能像两环模型里的 $delta Omega$ 那样只是一个整体频率偏移。两种读法给出的拓扑性质完全相反。
按 2022-03-30 那组参数，相变点在 $3 sqrt(3) t_2 abs(sin phi) = 2.0 times 10^(-4)$：
若 $delta nu = 9.82 times 10^(-4)$ 是质量项，模型已远在拓扑相之外，是平凡绝缘体、不该有边缘态；
若它只是整体偏移（$M approx 0$），则带隙宽度应为 $2 times 3 sqrt(3) t_2 abs(sin phi)
= 4.0 times 10^(-4)$，与 2022-04-13 记录的 $0.0009$–$0.0013$（宽 $4 times 10^(-4)$）相符。
后一种读法能同时解释带隙大小和边缘态的存在，但原稿没有明说，这里保留为待确认的问题。]

#fig("assets/figures/Haldane_theory.svg", width: 62%, caption: [Haldane 模型给出的能带（理论）。上、下两条体态能带之间是拓扑带隙，带隙宽度由 $t_2 sin phi$ 与 $delta nu$ 的竞争决定。]) <fig-haldane-theory>

=== 用 FEM/COMSOL 数据拟合能带

模型参数的标定思路是“先算能带、再拟合模型”：

+ 用有限元（本工作里是 COMSOL）算出环阵列（或有限宽条带）在第一布里渊区
  $Gamma - M - K - Gamma$ 路径上的本征频率；
+ 对每个 $k$ 点把上面的 2×2 哈密顿量对角化，得到理论上的两条能带
  $omega_\pm (k) = d_0(k) plus.minus sqrt(m(k)^2 + |f(k)|^2)$；
+ 以 $(t_1, t_2, phi, delta nu, cal(C))$ 为拟合参数（$cal(C) = plus.minus 1$ 只用来
  挑选带隙的拓扑符号），对本征频率做最小二乘拟合。

拟合出来的具体数值、以及拟合曲线与 COMSOL 数据的对比，见
`docs/results/02-ring-array.typ` 的“能带与带隙”一节。

=== 相关工作：Hafezi 2013

这一路设计直接来自 Hafezi 等人的硅基微环阵列（Nat. Photon. 7, 1003 (2013)）：
他们用 SOI 工艺做高 $Q$ 环形谐振腔（$Q > 1 times 10^4$），衬底是 220 nm 厚的 Si
加 2 µm 厚的埋氧（BOX）；波导截面取 510 nm $times$ 220 nm 是为了保证 TE 光在
$lambda approx 1.55$ µm 处*单模*传输；环与连接波导之间的空气间隙取 180 nm 与 200 nm；
转角半径取 6 µm 以压低弯曲损耗。该工作证明了光子体系里可以构造出等效规范场，
并把环阵列用作量子霍尔边缘态的平台；本项目的参数（$n = 3.4$、TE 单模、
以 $1.55$ µm 附近的本征频率为工作点）沿用同一套取值。

#fix[原稿把“转角半径”写成 “radius of the cornor”；此处改为 “corner”。]

#note[原稿这一节是阅读笔记的草稿，条目多为“Photons 好”“操”之类的片段；
已重写为上面这段正常综述，只保留其中的工艺与几何参数。]

== 模型与代码的对应

#code-ref[
  上面这些模型参数不是凭空给的，它们对应本项目求解器里的一条固定链路：
  “生成结构 → 解麦克斯韦方程 → 取本征频率/S 参数 → 用 Python 拟合”。逐段对应关系如下。

  *结构与网格*：`python/original_mesh.py` 里的 `cross_waveguide()`、`cross_waveguide_2()`
  用 `pygmsh` + `gmsh` 生成几何与二阶网格，边界由 Bézier 控制点给出
  （`Bezier.Bezier_interp()`、`Bezier.Bezier_tangent_vector()` 求曲线上的点与端点切矢），
  并写出 `data/cross_waveguide.msh`、`data/cross_waveguide.mat`
  （端口编号/坐标/朝向/入射端口，以及曲线的控制点与切矢）和 `data/infile.pot`（运行参数）。
  端口处的传播常数 $beta$（上面的相位 $k pi r$ 就来自它）由
  `python/propagation_constant0.py` 的 `get_beta_from_k0()` 用牛顿迭代解色散方程得到。

  *求解与 S 参数*：`waveguide_optimization/main.cpp` 读入上述文件，
  `mesh_snap_to_Bezier()` 把网格节点贴到 Bézier 曲线上，`assemble_waveguide()` 装配二维
  麦克斯韦系统，端口的解析模式场在 `waveguide_optimization/2D_waveguide_model_params.{h,cpp}`
  的 `analytic_E_field()` 里，求解后由 `integrate_along_side(..., mode_decompose)` 把
  每个端口的场投影到模式上，得到 $S$ 参数并写进 `cross_waveguide_result.mat` 的 `S_params_0`。

  *驱动与拟合*：`python/optimize.py` 的 `get_dof_derivative()` 负责“重建几何 → 调求解器 →
  读回 `S_params_0`”这一循环，`get_S_params()` / `gradient_descent()` 再把 $S$ 参数、
  设计变量与灵敏度写成 CSV；
  `docs/theory/06-coupled-resonators.typ` 给出把这些 $S$ 参数换算成
  $(t, kappa_("ex"), kappa_("in"), phi)$ 的公式。环阵列的 $t_1$、$t_2$ 与 Haldane 参数
  就是在同样的“$S$ 参数/本征频率 → 最小二乘”这一步上标定的。
]

#note[需要说清楚的是：2022 年 3–4 月这一批微环阵列的本征频率与 $S$ 参数实际是用 COMSOL
算出来的（本文与结果文档里带 COMSOL 水印的几何/数据图都是），Haldane 模型的拟合脚本也未随
本仓库入仓。本仓库提供的是*同构*的一条链路——`python/original_mesh.py` 建网格 →
`waveguide_optimization` 解场与取 $S$ 参数 → `python/optimize.py` 在 Python 侧拟合/优化。
把环阵列接进这条链路，只需要在 `original_mesh.py` 里换掉几何构造函数，
其余（`.mat`/`infile.pot` 接口、`mode_decompose()`、结果 `.mat` 的字段名）都不用动。]

== 与其他文档的关系

+ `docs/theory/06-coupled-resonators.typ`：两环耦合的 input–output 形式，
  $S$ 参数与 $(t, kappa_("ex"), kappa_("in"), phi)$ 的换算公式；
+ `docs/theory/06-coupled-resonators.typ`：连接环/三角连接环的绝热消除，
  即本节 $t_2 e^(i phi)$ 的完整推导；
+ `docs/results/02-ring-array.typ`：本文模型在 2022 年 3–4 月的全部数值结果。
