// ============================================================
//  Information Shapes Koopman Representation 详尽学习笔记
//  基于 Cheng, Yuan, Yang, Zhang, Cheng, He, Sun (ICLR 2026)
//  arXiv:2510.13025
// ============================================================

#import "@preview/cetz:0.3.4": canvas, draw

#set page(
  paper: "a4",
  margin: (x: 2.1cm, y: 2.3cm),
  numbering: "1",
  header: context {
    let h = counter(heading).get()
    if h.len() > 0 {
      align(right, text(size: 8.5pt, fill: luma(110))[
        《Information Shapes Koopman Representation》精读笔记
      ])
    }
  },
)

#set text(font: ("New Computer Modern", "Noto Serif SC"), size: 10.3pt, lang: "zh", region: "cn")
#show math.equation: set text(font: ("New Computer Modern Math", "Noto Serif SC"))
#set par(justify: true, leading: 0.78em, first-line-indent: 1.6em)
#set heading(numbering: "1.1")

#import "common.typ": *

#align(center)[
  #v(1.2cm)
  #text(size: 21pt, weight: "bold")[信息如何塑造 Koopman 表示]
  #v(0.2cm)
  #text(size: 13pt)[Information Shapes Koopman Representation]
  #v(0.7cm)
  #text(size: 10.5pt)[
    Xiaoyuan Cheng, Wenxuan Yuan, Yiming Yang, Yuanzhao Zhang,\
    Sibo Cheng, Yi He, Zhuo Sun
  ]
  #v(0.15cm)
  #text(size: 9.5pt, fill: luma(90))[ICLR 2026 · arXiv:2510.13025]
  #v(1.1cm)
]

#align(center)[
  #box(width: 78%)[
    #text(size: 9.8pt, style: "italic")[
      一份面向"有机器学习基础、但从未接触过 Koopman 算子或信息瓶颈"的读者的\
      自足讲解。目标不是复述论文，而是把论文的推理链拆开、补全，\
      让读者读完能自己判断这篇工作在做什么、站在哪里、下一步该往哪走。
    ]
  ]
]

#v(0.8cm)

#boxed(title: "怎么读这份笔记")[
  这份笔记刻意不按"摘要—引言—方法—实验—结论"的顺序走。论文的论证主线是一条
  *信息在潜在空间中如何流动、流失、被塑造* 的线索，笔记就沿着这条线索走：

  #v(0.3em)
  - *第 0 章* 预备知识。默认读者只学过数学分析、线性代数、概率论与常微分方程；
    信息论、动力系统、自编码器这些工具在第 0 章按同一套路补齐
    （先讲它回答什么问题，再给数值例子，再给定义）。已学过对应内容的读者可以跳过，
    只读末尾的阅读地图与自测。
  - *第一部分（§1–§2）* 两条线索的谱系。先把问题立起来：为什么需要 Koopman 线性化，
    为什么它天然困难，深度学习接手后具体坏在哪里；再把信息瓶颈这条线的工具备齐。
    读完你会知道这*不是一个工程调参问题*。
  - *第二部分（§3–§5）* 是核心。沿着论文自己的提问链
    「信息会丢吗？→ 丢多少？→ 丢的是什么？→ 怎么优化？→ 怎么防止副作用？」
    逐步证明五个命题。每个命题都给完整的推导，不跳步。
  - *第三部分（§6）* 把理论翻译成可运行的损失函数与两份算法。
  - *第四部分（§7）* 实验如何反过来验证理论。
  - *第五部分（§8）* 横向坐标：这条线索相对于 Koopman 谱系和 IB 谱系的位置。
  - *第六部分（§9–§10）* 局限、下一步，以及记号速查与参考文献。

  数学密度较高，但每一个符号第一次出现时都会先解释*它为什么存在*。遇到
  #text(fill: rgb("#b8860b"))[黄色框] 时，那是最需要记住的一句话；
  遇到 #text(fill: rgb("#1a7a45"))[绿色框] 是直觉总结；
  #text(fill: rgb("#a33"))[红色框] 是容易出错的细节。
]

#v(0.5cm)

#boxed(title: "术语地图（第一次遇到不懂的词，先回来查这里）", fill: rgb("#f7f7fb"), stroke: rgb("#b9b9d6"))[
  下面这些词在正文里第一次出现时都会展开解释，这里先给一句最简版，方便随时回查。

  #v(0.3em)
  #table(
    columns: (auto, 1fr),
    inset: (x: 6pt, y: 3pt),
    align: (left, left),
    stroke: none,
    [*可观测量*], [状态的函数，比如“温度场”“某个区域的平均速度”。Koopman 理论作用在它们上面。],
    [*谱 / 特征值*], [线性算子“各个方向上的缩放倍数”。特征值模长 $<1$ 就衰减，$=1$ 就持续，$>1$ 就爆炸。],
    [*不变子空间*], [一个函数集合，被算子作用之后仍然留在自己里面。只有在这种集合上，算子才能写成有限矩阵。],
    [*潜变量 / 潜在空间*], [把高维状态压缩之后得到的低维向量，以及它所在的空间。],
    [*互信息*], [两个变量“共享了多少信息”。用它可以量化“$z$ 里还留着多少关于未来 $x$ 的线索”。],
    [*率失真*], [“码率”与“误差”之间的最优兑换关系。这里用来把“信息量”换算成“误差有多大”。],
    [*注水*], [一种资源分配方式：水（方差预算）先流向门槛最低的容器（增益最大的方向）。低增益方向一滴都拿不到。],
    [*密度矩阵*], [半正定、迹为 1 的矩阵。它的特征值构成一个概率分布，可以理解成“方差在各方向的占比”。],
    [*von Neumann 熵*], [上面那个特征值分布的香农熵。它只关心“分布铺得开不开”，不关心总量。],
    [*有效维度*], [$exp("熵")$。回答“这些方差实际上用到了几个方向”，取值在 1 到 $d$ 之间。],
    [*拉格朗日乘子*], [把带约束的优化改写成无约束优化的辅助变量，负责“守住约束”。],
    [*变分界 / ELBO*], [真实目标算不出来时，用一个能算的量从上方或下方夹住它。深度生成模型的标准工具。],
  )
]

#v(0.6cm)
#outline(title: [目录], indent: 1.4em, depth: 2)

#pagebreak()

#include "parts/00_foundation.typ"
#include "parts/01_problem.typ"
#include "parts/02_tools.typ"
#include "parts/03_q12.typ"
#include "parts/04_q3.typ"
#include "parts/05_q45.typ"
#include "parts/06_algorithm.typ"
#include "parts/07_experiments.typ"
#include "parts/08_related.typ"
#include "parts/09_future.typ"
#include "parts/10_appendix.typ"
