#import "@preview/cetz:0.3.4": canvas, draw

// 共享 preamble 模板：三篇郭玲笔记各自 include
#let boxed(content, fill: rgb("#f4f8ff"), stroke: rgb("#c3d5ef"), title: none) = block(
  width: 100%,
  inset: 10pt,
  radius: 4pt,
  fill: fill,
  stroke: (left: 2.6pt + stroke),
)[
  #if title != none [#text(weight: "bold")[#title] #h(0.4em)]
  #content
]

#let keypoint(title, body) = boxed(body, fill: rgb("#fff9ec"), stroke: rgb("#ecc76a"), title: title)
#let warn(title, body) = boxed(body, fill: rgb("#fdf1f1"), stroke: rgb("#e0a3a3"), title: title)
#let insight(title, body) = boxed(body, fill: rgb("#f2fbf5"), stroke: rgb("#8dc9a3"), title: title)

#let definition(name, body) = block(
  width: 100%,
  inset: (x: 11pt, y: 8pt),
  radius: 3pt,
  stroke: (left: 2.6pt + rgb("#7a7a7a")),
  fill: rgb("#fafafa"),
)[
  *#name* #h(0.3em) #body
]

#let proposition(name, body) = block(
  width: 100%,
  inset: (x: 11pt, y: 8pt),
  radius: 3pt,
  stroke: (left: 3pt + rgb("#1a4f9a")),
  fill: rgb("#f4f8ff"),
)[
  *#name* #h(0.3em) #body
]

#let lemma(name, body) = block(
  width: 100%,
  inset: (x: 11pt, y: 7pt),
  radius: 3pt,
  stroke: (left: 2.4pt + rgb("#6a4c93")),
  fill: rgb("#faf7ff"),
)[
  *#name* #h(0.3em) #body
]

#let proof(body) = block(
  width: 100%,
  inset: (x: 11pt, y: 6pt),
  radius: 2pt,
  fill: rgb("#fcfcfc"),
  stroke: (left: 1pt + luma(190)),
)[
  #text(style: "italic")[证明.] #body
]

#let term(en, zh) = [#text(fill: rgb("#1a4f9a"))[#en]（#zh）]

#let figure-caption(body) = block(width: 100%, inset: (top: 5pt), text(size: 9pt, fill: luma(70))[#body])

// CeTZ 用的节点盒
#let cnode(w, h, c, body, size: 8.5pt) = rect(
  width: w, height: h, radius: 3pt, stroke: 1pt + c, fill: c.lighten(88%),
)[#align(center + horizon, text(size: size)[#body])]

// 简洁表格
#let ttable(header, rows, widths: auto, align-x: left) = table(
  columns: widths,
  inset: 6pt,
  align: align-x,
  stroke: 0.4pt + luma(190),
  table.header(..header.map(h => { set par(justify: false); h })),
  ..rows.map(r => r.map(c => { set par(justify: false); c })).flatten(),
)
