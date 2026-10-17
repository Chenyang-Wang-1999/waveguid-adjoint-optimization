#import "template.typ": project-doc, note, fix, code-ref

#show: project-doc.with(
  title: "文档索引",
  subtitle: "交叉波导的伴随法形状优化",
  date: "2026 年整理",
)

== 这个项目是什么

用有限元 + 伴随法优化光子交叉波导的交叉口形状，把直通端透射率从
$|S_(21)|^2 = 0.712$ 提升到 $approx 0.96$。设计变量是交叉口边界上三次 Bézier 曲线的
控制点坐标，梯度由一次伴随求解得到。代码写在 2022 年 2 月，2022 年 3–4 月又用同一套
求解器做了微环阵列的能带与边缘态计算。

代码结构与运行方式看 `overview.typ` 和 `build-and-run.typ`；仓库首页是根目录的
`README.md`。

== 文档地图

=== 总览

#table(
  columns: (3.6cm, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*内容*]),
  [`overview.typ`], [项目总览：物理问题、单位制、几何参数化、数值方法、代码结构、数据流、依赖],
  [`build-and-run.typ`], [依赖与配置、`make` 各目标、`infile.pot` 参数表、想改点什么时改哪里],
)

=== 理论与方法（`theory/`）

#table(
  columns: (3.6cm, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*内容*]),
  [`01-fem-basics.typ`], [有限元方法基础：1D 模型问题的变分形式、分片线性函数、误差估计、Poisson 方程],
  [`02-maxwell-2d.typ`], [二维麦克斯韦方程、TE/TM 分解、三类边界条件、金属波导的端口透明边界与 FEM 形式、平板波导解析解],
  [`03-slab-waveguide.typ`], [平板波导传播常数：色散方程、牛顿迭代、把 `tan` 换成 `arctan` 的数值稳定性处理],
  [`04-bezier.typ`], [Bézier 边界参数化：$C^2$ 连续的拼接条件、端点条件、插值、控制点灵敏度],
  [`05-adjoint-sensitivity.typ`], [伴随法灵敏度分析：伴随方程、梯度公式、右端项、代价对比、有限差分验证],
  [`06-coupled-resonators.typ`], [耦合微环：input–output formalism 推出 S 参数、耦合系数的拟合公式、连接环模式的绝热消除],
  [`07-bands-and-edge-states.typ`], [微环阵列的能带与边缘态模型：环间耦合矩阵、次近邻耦合、Haldane 有效模型],
)

=== 工程笔记（`notes/`）

#table(
  columns: (3.6cm, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*内容*]),
  [`libmesh-notes.typ`], [libMesh 的编译链接、Mesh/Element/`FEMContext`、含时问题求解，以及本项目用到这些 API 的地方],
)

=== 结果（`results/`）

#table(
  columns: (3.6cm, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*文件*], [*内容*]),
  [`01-cross-waveguide.typ`], [交叉波导优化：模型参数、优化前后对比、遇到的问题与改进、场分布与 S 参数],
  [`02-ring-array.typ`], [微环阵列：$t_1$/$t_2$ 与 $kappa$ 的拟合值、相位校核、能带与带隙、边缘态、两个赝自旋的耦合],
)

== 编译

```bash
./docs/build.sh                              # 全部编译到 docs/build/
typst compile --root . docs/overview.typ out.pdf    # 单篇
typst watch   --root . docs/theory/05-adjoint-sensitivity.typ
```

> *必须加 `--root .`*：Typst 默认把“输入文件所在目录”当作项目根，而各篇文档都用
> `../template.typ` 引用上层的模板，不加这个参数会被判成越界访问。

排版全部集中在 `docs/template.typ` 里。正文只依赖它导出的接口：

```typst
#import "../template.typ": project-doc, fig, fig-row, fig-gallery, note, fix, code-ref
#show: project-doc.with(title: "标题", date: "2022-03")
```

想换版式时，只要保持这些名字和签名不变、替换 `template.typ` 即可，正文一行都不用动。

除了排版，模板还导出三个写作辅助：

#table(
  columns: (3.6cm, 1fr),
  align: (left, left),
  stroke: none,
  inset: (x: 0.5em, y: 0.4em),
  table.header([*接口*], [*用途*]),
  [`note[...]`], [编者注：需要额外说明、但不属于正文的地方],
  [`fix[...]`], [修正说明：正文写的是修正后的内容，这里注明原稿写法与修改理由],
  [`code-ref[...]`], [代码对应：把这一节的公式或结论映射到具体文件与函数],
)

#note[
  *图片路径的写法。* 所有插图路径都写成相对 `docs/` 的形式（例如 `assets/figures/x.png`），
  不写 `../assets/...`。原因是 `fig()` / `fig-row()` / `fig-gallery()` 定义在
  `docs/template.typ` 里，Typst 是按“路径字符串出现的那个文件”来解析相对路径的，
  也就是按 `docs/` 解析。自己写 `image()` 时请沿用同一条约定。
]

== 图片

原稿（LyX 写的笔记与周报）里的插图散落在项目之外的 `CouplingRings/` 目录下。整理时把
被正文引用到的 55 张图复制到了 `docs/assets/`，保持了原有的相对目录结构：

```text
docs/assets/
├── figures/      交叉波导与微环实验的通用图（结构、S 参数、拟合、能带）
├── mode-m26/     m = 26 那组结构参数下的能带图
├── adiabatic/    绝热消除用的几何示意图
├── py/           Python 侧画出来的图
└── optimize/     优化过程中的结构与场分布
```

少数体积过大的位图（> 200 kB）等比缩小到最长边 1200 px 并重新编码：线稿类存成调色板 PNG，
连续场云图存成 JPEG；矢量图（SVG）原样保留。重写文档时又筛掉了只起过程记录作用、
对结论没有贡献的图，所以实际保留的比原稿引用的少得多。

== 约定

- 全文采用 *$epsilon_0 = mu_0 = 1$* 的单位制，因此所有物理量都量纲化为长度
  （详见 `theory/02-maxwell-2d.typ`）。
- 时间依赖取 $e^(-i omega t)$，因此 $diff/bold(D) = -i omega bold(D)$。
- $tilde(a)$ 表示旋转坐标系（rotating frame）中的算符，$a$ 表示实验室系算符。
- 频率、长度等数值都是归一化单位下的值，不带量纲。

#note[
  *整理说明。* 这些文档不是原稿的影印，而是重写：只保留了对理解项目和代码有帮助的内容，
  删掉了走不通的推导、重复的公式和过程性的中间结果。原稿里明显抄错或不自洽的地方做了修正，
  每一处都用 `#fix[...]` 注明了原稿的写法与修改理由，方便对照。原始 LyX 文件没有随仓库发布。
]
