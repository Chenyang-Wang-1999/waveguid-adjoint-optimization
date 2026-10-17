// ============================================================================
//  docs/template.typ —— 交叉波导伴随优化项目 / 文档统一模板
// ----------------------------------------------------------------------------
//  这个文件只管“排版”，所有正文内容都在各个 .typ 文件里。
//  如果你有自己的模板，只需要保证导出同样的接口即可整体替换：
//
//      #import "../template.typ": project-doc, fig, note
//      #show: project-doc.with(title: "标题", date: "2022-03-09")
//
//  用法：
//      typst compile docs/theory/01-fem-basics.typ
//      typst watch   docs/theory/01-fem-basics.typ
// ============================================================================

#let project-doc(
  title: [],
  subtitle: none,
  date: none,
  author: "王晨阳 (Chenyang Wang)",
  abstract: none,
  body,
) = {
  set document(title: title, author: author)

  set page(
    paper: "a4",
    margin: (x: 2.4cm, top: 2.4cm, bottom: 2.4cm),
    numbering: "1",
    number-align: center,
    header: context {
      let page = counter(page).get().first()
      if page > 1 {
        set text(size: 8.5pt, fill: luma(45%))
        grid(
          columns: (1fr, auto),
          align: (left, right),
          title,
          [#page],
        )
        v(-0.6em)
        line(length: 100%, stroke: 0.4pt + luma(75%))
      }
    },
  )

  // 中文字体。本机可用：Noto Serif CJK SC（思源宋体）/ Noto Sans CJK SC（思源黑体）。
  // 换成你自己装的字体即可，例如 ("Source Han Serif SC", "Latin Modern Roman")。
  set text(
    font: ("Noto Serif CJK SC", "New Computer Modern"),
    size: 10.5pt,
    lang: "zh",
    region: "cn",
  )
  set par(justify: true, leading: 0.78em, spacing: 1.05em)

  // 编号：文档标题用自定义标题块（不占 heading 层级），因此正文从 `==` 开始。
  // 直接写 numbering: "1.1" 会让第一层显示成 “0.1”，这里把前导的 0 去掉。
  set heading(numbering: (..nums) => {
    let n = nums.pos()
    let out = ()
    let started = false
    for v in n {
      if v != 0 { started = true }
      if started { out.push(v) }
    }
    if out.len() == 0 { return none }
    numbering("1.1", ..out)
  })
  show heading.where(level: 1): it => block(above: 1.6em, below: 0.9em)[
    #set text(size: 15pt, weight: "bold")
    #it
  ]
  show heading.where(level: 2): it => block(above: 1.3em, below: 0.7em)[
    #set text(size: 12.5pt, weight: "bold")
    #it
  ]
  show heading.where(level: 3): it => block(above: 1.1em, below: 0.6em)[
    #set text(size: 11pt, weight: "bold")
    #it
  ]
  show raw.where(block: true): it => block(
    width: 100%,
    inset: 0.6em,
    radius: 2pt,
    fill: luma(97%),
    stroke: 0.4pt + luma(85%),
    it,
  )

  // ---- 标题区 -------------------------------------------------------------
  block(above: 0em, below: 1.8em)[
    #text(size: 19pt, weight: "bold")[#title]
    #if subtitle != none [
      #par(leading: 0.9em)[#text(size: 12pt, fill: luma(35%))[#subtitle]]
    ]
    #v(0.5em)
    #line(length: 100%, stroke: 0.8pt)
    #v(0.35em)
    #grid(
      columns: (1fr, auto),
      align: (left, right),
      text(size: 9.5pt, fill: luma(35%))[#author],
      if date != none { text(size: 9.5pt, fill: luma(35%))[#date] },
    )
  ]

  if abstract != none [
    #block(
      width: 100%,
      inset: (x: 1em, y: 0.8em),
      fill: luma(96%),
      radius: 3pt,
    )[
      #text(size: 9.5pt)[*摘要。* #abstract]
    ]
    #v(1em)
  ]

  body
}

// 插图：`#fig("../assets/figures/x.png", caption: [...], width: 65%)`
#let fig(path, caption: none, width: 70%, height: auto) = figure(
  image(path, width: width, height: height),
  caption: caption,
  kind: image,
  supplement: "图",
)

// 并排插图
#let fig-row(paths, caption: none, width: 46%) = figure(
  grid(
    columns: paths.len(),
    column-gutter: 1.2em,
    ..paths.map(p => image(p, width: width)),
  ),
  caption: caption,
  kind: image,
  supplement: "图",
)

// 图组：把一组图按每行 `per-row` 张排成多行，整组共用一个图注、占一个图号。
// 用于原文一个 Float 里放多张图的情形（此时用多个 `fig-row` 会散成好几张图）。
#let fig-gallery(paths, per-row: 2, caption: none, width: 92%) = figure(
  stack(
    spacing: 0.7em,
    ..paths.chunks(per-row).map(row => grid(
      columns: row.len(),
      column-gutter: 1.2em,
      ..row.map(p => image(p, width: width)),
    )),
  ),
  caption: caption,
  kind: image,
  supplement: "图",
)

// 编者注 / 待办：整理文档时留下的说明，排版时可整体删掉。
#let note(body) = block(
  width: 100%,
  inset: (x: 0.9em, y: 0.6em),
  stroke: (left: 2pt + rgb("#c8a24a")),
  fill: rgb("#fdf8ec"),
  radius: 2pt,
)[
  #set text(size: 9.5pt)
  #body
]

// 修正说明：正文里写的是*修正后*的内容，这里注明原稿的写法与修改理由。
// 用法：`#fix[原稿写作 $x = y$，此处按 … 修正。]`
#let fix(body) = block(
  width: 100%,
  inset: (x: 0.9em, y: 0.6em),
  stroke: (left: 2pt + rgb("#4a7fb5")),
  fill: rgb("#eef4fb"),
  radius: 2pt,
)[
  #set text(size: 9pt)
  *修正。* #body
]

// 代码对应：把这一节的公式/结论映射到具体文件与函数。
// 用法：`#code-ref[(waveguide_optimization/main.cpp)[assemble_waveguide()]]`
#let code-ref(body) = block(
  width: 100%,
  inset: (x: 0.9em, y: 0.6em),
  stroke: (left: 2pt + luma(60%)),
  fill: luma(97%),
  radius: 2pt,
)[
  #set text(size: 9.5pt)
  *代码对应。* #body
]
