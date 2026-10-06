// ============================================================
//  稀疏多项式逼近与压缩感知配点法：郭玲团队方向 A 精读笔记
// ============================================================

#set page(
  paper: "a4",
  margin: (x: 2.1cm, y: 2.3cm),
  numbering: "1",
  header: context {
    let h = counter(heading).get()
    if h.len() > 0 {
      align(right, text(size: 8.5pt, fill: luma(110))[
        稀疏多项式逼近与压缩感知配点法 · 精读笔记
      ])
    }
  },
)

#set text(font: ("New Computer Modern", "Noto Serif SC"), size: 10.3pt, lang: "zh", region: "cn")
#show math.equation: set text(font: ("New Computer Modern Math", "Noto Serif SC"))
#set par(justify: true, leading: 0.78em, first-line-indent: 1.6em)
#set heading(numbering: "1.1")

#import "../preamble.typ": *

#align(center)[
  #v(1.2cm)
  #text(size: 21pt, weight: "bold")[稀疏多项式逼近与压缩感知配点法]
  #v(0.2cm)
  #text(size: 13pt)[一条由 Christoffel 函数贯穿的研究线]
  #v(0.7cm)
  #text(size: 10.5pt)[
    精读对象：Guo, Narayan, Zhou, Chen (2017)；Guo, Narayan, Yan, Zhou (2017)；\
    Guo, Narayan, Zhou (2018)；Guo, Li, Liu (2018)；Guo, Liu, Zhou (2018)；\
    Guo, Narayan, Liu, Zhou (2020)
  ]
  #v(0.15cm)
  #text(size: 9.5pt, fill: luma(90))[上海师范大学 郭玲 · 不确定性量化方向 · arXiv 1602.00995 / 1708.01296 / 1802.08837 / 1805.05416 / 1805.10893 / 2008.10121]
  #v(1.1cm)
]

#align(center)[
  #box(width: 78%)[
    #text(size: 9.8pt, style: "italic")[
      本笔记面向已修完数学分析、线性代数、概率论与常微分方程、
      没接触过不确定性量化的读者。
      问题只有一个：一个 PDE 的系数是随机的，解的统计量怎么算才不必把随机维数
      的网格铺满。六篇论文给出同一个答案的六个版本：利用稀疏性，
      用远少于基函数个数的样本恢复展开系数。贯穿全部六篇的数学对象是
      Christoffel 函数，它在每篇里换一个身份出现：求积权、点列加权、
      列归一、预条件、采样测度。
    ]
  ]
]

#v(0.8cm)

#boxed(title: "怎么读这份笔记")[
  - *第 0 章* 预备知识：正交多项式、条件数、最小二乘、稀疏性与采样分布
    这些工具集中在一章里补齐，每个概念从它回答什么问题讲起；
    正文遇到陌生符号时先去该章查，章末有阅读地图与自测问题。
  - *第 1 到 5 章* 问题与工具：随机方程的维数灾难，正交多项式与 Gauss 求积，
    稀疏恢复的几何图像（基追踪、RIP、互不相干），以及三条主定理的证明思路。
    Christoffel 函数在第 1 章末尾定义，它是整条线的公共主轴。
  - *第 6 到 8 章* 是前半条线：怎么采样。Gauss 求积网格的随机子采样（第 6 章），
    确定性的加权 Fekete 点（第 7 章），第 8 章补齐两篇的实验与一维最优性推导。
  - *第 9 到 12 章* 加信息与换罚函数：把导数塞进测量矩阵（第 9 章），
    非凸 TL1（第 10 章），证明链与实现细节（第 11、12 章）。
  - *第 13 到 17 章* 换数据假设：分布未知时的数据驱动基（第 13 章），
    稀疏恢复与数据驱动的合流 induced sampling（第 14 章），
    两篇各自的实现台（第 15 到 17 章）。
  - *第 18 到 22 章* 用法与坐标：最小算例走查（第 18 章），符号对照（第 19 章），
    六篇算例与选型规则（第 20 章），前人路线谱系（第 21 章），
    这条线的立场与未做之事（第 22 章）。第 23、24 章是附录。

  每个术语第一次出现时先说明它解决什么困难。公式后面的“读法”段落
  把符号翻译回问题本身。
]

#v(0.5cm)

#boxed(title: "术语地图", fill: rgb("#f7f7fb"), stroke: rgb("#b9b9d6"))[
  #table(
    columns: (auto, 1fr),
    inset: (x: 6pt, y: 3pt),
    align: (left, left),
    stroke: none,
    [*gPC*], [广义多项式混沌：把随机系统的解展开成随机变量的正交多项式。],
    [*随机配点*], [在随机空间选一批点，每点解一次确定性 PDE，再把解投影回多项式基。],
    [*维数灾难*], [张量积配点数随维数指数增长；本文所有方法的对手。],
    [*稀疏性*], [真实解的 gPC 系数大多接近零，只有少数几个大系数。],
    [*基追踪*], [在数据约束下最小化系数向量的 $ell_1$ 范数，用来逼出最稀疏解。],
    [*RIP*], [有限等距性质：采样矩阵近似保持稀疏向量的长度，是恢复保证的载体。],
    [*互不相干*], [测量矩阵列与列的归一化内积的最大值；越小恢复越容易。],
    [*Christoffel 函数*], [多项式平方和的倒数。它衡量“在某个点附近多项式空间有多挤”。],
    [*求积权*], [数值积分公式里每个节点的权重。],
    [*平衡测度*], [位势论意义下让多项式“均匀分布”的采样密度；有界域上是 Chebyshev 测度。],
    [*Fekete 点*], [让 Vandermonde 行列式最大的点集，插值意义下最“张得开”的点。],
    [*诱导测度*], [由“原密度除以 Christoffel 函数”定义的采样分布，只依赖已有样本。],
  )
]

#v(0.6cm)
#outline(title: [目录], indent: 1.4em, depth: 2)

#pagebreak()

#include "parts/a00_foundation.typ"
#include "parts/a01_problem.typ"
#include "parts/a01b_prep.typ"
#include "parts/a02_cs.typ"
#include "parts/a02c_theory.typ"
#include "parts/a03_p1_p2.typ"
#include "parts/a03b_more.typ"
#include "parts/a04_p3_p4.typ"
#include "parts/a04b_more.typ"
#include "parts/a03c_impl.typ"
#include "parts/a05_p5_p6.typ"
#include "parts/a05b_more.typ"
#include "parts/a05d_p5_lab.typ"
#include "parts/a05e_p6_lab.typ"
#include "parts/a05c_practice.typ"
#include "parts/a06b_case.typ"
#include "parts/a06_related.typ"
#include "parts/a07_appendix.typ"
#include "parts/a07_extra.typ"
