#import "@preview/cetz:0.3.4": canvas, draw

// ============================================================
//  用神经网络解随机微分方程 —— 郭玲团队方向 B 精读笔记
//  论文线：NN-aPC (2019) → NN-DO/BO (2020) → NFF (2022)
//          → MC-PINN (2022) → BO-fPINN (2023) → sPI-GeM (2025)
// ============================================================

#set page(
  paper: "a4",
  margin: (x: 2.1cm, y: 2.3cm),
  numbering: "1",
  header: context {
    let h = counter(heading).get()
    if h.len() > 0 {
      align(right, text(size: 8.5pt, fill: luma(110))[
        用神经网络解随机微分方程 · 精读笔记
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
  #text(size: 21pt, weight: "bold")[用神经网络解随机微分方程]
  #v(0.2cm)
  #text(size: 13pt)[从物理信息网络到深度生成模型的一条研究线]
  #v(0.7cm)
  #text(size: 10.5pt)[
    精读对象：Zhang, Lu, Guo, Karniadakis (2019)；Zhang, Guo, Karniadakis (2020)；\
    Guo, Wu, Zhou (2022)；Guo, Wu, Yu, Zhou (2022)；Ma, Li, Zeng, Guo, Karniadakis (2023)；\
    Zhou, You, Guo, Meng (2025)
  ]
  #v(0.15cm)
  #text(size: 9.5pt, fill: luma(90))[上海师范大学 郭玲 · 科学机器学习方向 · arXiv 1809.08327 / 1905.01205 / 2108.12956 / 2203.08501 / 2303.10913 / 2503.18012]
  #v(1.1cm)
]

#align(center)[
  #box(width: 78%)[
    #text(size: 9.8pt, style: "italic")[
      本笔记面向“有机器学习基础、解过常微分方程、但没接触过不确定性量化
      与物理信息网络”的读者。它沿一条线索走：神经网络怎么被用来解
      带随机性的偏微分方程，每一步是为了克服上一步的哪个具体困难。\
      六篇论文以失败模式串成一条链：每篇都以解决前一篇暴露的困难为出发点。
    ]
  ]
]

#v(0.8cm)

#boxed(title: "怎么读这份笔记")[
  这份笔记不按“摘要、方法、实验”的顺序逐篇复述，而是按论文自身的递进逻辑走：

  #v(0.3em)
  - *第 1 到 3 章* 备齐基础：随机方程为什么难算，物理信息网络（PINN）
    的基本构造，以及 PINN 类损失的三个训练实践问题。没有这三块，后面的每一篇都读不动。
  - *第 4 到 9 章* 是前半条主线。NN-aPC 把随机性第一次放进 PINN 的损失（第 4、5 章）；
    NN-DO/BO 把“模态空间”的约束写成弱形式损失，克服长时间积分的困难（第 6、7 章）；
    NFF 用可逆变换给出严格的概率密度，替代近似式的模态展开（第 8、9 章）。
  - *第 10 到 13 章* 是后半条主线。MC-PINN 用随机采样估计分数阶算子，
    解锁高维分数阶方程（第 10、11 章）；BO-fPINN 把模态弱形式与分数阶离散拼在一起（第 12 章）；
    sPI-GeM 把生成模型降到“系数的分布”上，把可扩展性推到空间二十维（第 13 章）。
  - *第 14 到 19 章* 横向坐标与下一步：六方法总对比与选型（第 14 章），
    实验画廊（第 15、18 章），谱系与立场（第 16、17 章），附录（第 19 章）。

  每个符号第一次出现时都会先说明它为什么存在。
  #text(fill: rgb("#b8860b"))[黄色框] 是需要记住的一句话；
  #text(fill: rgb("#1a7a45"))[绿色框] 是直觉总结；
  #text(fill: rgb("#a33"))[红色框] 是容易出错的细节。
]

#v(0.5cm)

#boxed(title: "术语地图（第一次遇到不懂的词，先回来查这里）", fill: rgb("#f7f7fb"), stroke: rgb("#b9b9d6"))[
  下面这些词在正文里第一次出现时都会展开解释，这里先给一句最简版。

  #v(0.3em)
  #table(
    columns: (auto, 1fr),
    inset: (x: 6pt, y: 3pt),
    align: (left, left),
    stroke: none,
    [*随机场*], [定义在空间点上、取值随机的函数，比如“渗透率场 $k(x, omega)$”。],
    [*维数灾难*], [需要的计算量随随机维数指数增长的现象。几乎所有新方法的动机都是绕开它。],
    [*KL 展开*], [把随机场按“方差从大到小”分解成不相关随机变量加权的空间模态之和。],
    [*gPC / aPC*], [把解展开成随机变量的正交多项式。gPC 需要已知分布，aPC 只需要数据。],
    [*PINN*], [物理信息网络：把“方程残差在配点上接近零”写进神经网络的损失函数。],
    [*自动微分*], [框架（PyTorch 等）对网络输出求精确导数的能力。分数阶导数是它够不着的。],
    [*弱形式*], [把方程两边乘测试函数再积分，降低对光滑性的要求，也降低自动微分的阶数。],
    [*DO / BO*], [让时变模态展开“不冗余”的两类正交约束；经典算法要显式解演化方程，NN 版本把它们写进损失。],
    [*特征值交叉*], [两个模态方差曲线相撞。经典 BO 方法在此失效，是 NN 版本的直接动机。],
    [*normalizing flow*], [可逆变换链：把简单分布精确地搬成复杂分布，密度可以逐层算出来。],
    [*分数阶导数*], [积分形式的求导，描述记忆与长程相互作用；非局部，无法自动微分。],
    [*MC 估计*], [用随机采样近似一个积分或期望，样本平均在期望意义下等于真值（无偏）。],
    [*WGAN*], [生成对抗网络的一种，训练目标是让生成分布与真实分布的 Wasserstein 距离变小。],
    [*置信区间校准*], [预测的不确定度与真实误差相匹配：95% 的区间应该盖住约 95% 的真值。],
  )
]

#v(0.6cm)
#outline(title: [目录], indent: 1.4em, depth: 2)

#pagebreak()

#include "parts/b01_background.typ"
#include "parts/b02_pinn.typ"
#include "parts/b02b_training.typ"
#include "parts/b03_nnapc.typ"
#include "parts/b03b_nnapc_more.typ"
#include "parts/b04_dobo.typ"
#include "parts/b04b_dobo_more.typ"
#include "parts/b05_nff.typ"
#include "parts/b05b_nff_more.typ"
#include "parts/b06_mcpinn.typ"
#include "parts/b06b_mc_more.typ"
#include "parts/b07_bofpinn.typ"
#include "parts/b08_spigem.typ"
#include "parts/b08b_compare.typ"
#include "parts/b08c_gallery.typ"
#include "parts/b09_related.typ"
#include "parts/b09c_final.typ"
#include "parts/b10_appendix.typ"
