// ============================================================
//  科学机器学习的不确定性量化 —— 郭玲团队方向 C 精读笔记
// ============================================================

#set page(
  paper: "a4",
  margin: (x: 2.1cm, y: 2.3cm),
  numbering: "1",
  header: context {
    let h = counter(heading).get()
    if h.len() > 0 {
      align(right, text(size: 8.5pt, fill: luma(110))[
        科学机器学习的不确定性量化 · 精读笔记
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
  #text(size: 21pt, weight: "bold")[科学机器学习的不确定性量化]
  #v(0.2cm)
  #text(size: 13pt)[从一个综述到一个信息瓶颈框架，再到四个后继]
  #v(0.7cm)
  #text(size: 10.5pt)[
    精读对象：Psaros, Meng, Zou, Guo, Karniadakis (SIAM Review 2023)；\
    Guo, Wu, Zhou, Wang, Zhou (2023)；Wang, Guo, Wu, Zhou (2024)；\
    Zhang, Guo, Wu (2024)；Ma, Guo, Wu, Zhou (2025)；Feng, Guo, Wan, Wu, Zhou, Zhou (2025)
  ]
  #v(0.15cm)
  #text(size: 9.5pt, fill: luma(90))[上海师范大学 郭玲 · SciML 不确定性量化方向 · arXiv 2201.07766 / 2302.03271 / 2401.02080 / 2408.15981 / 2509.25646 / 2507.22493]
  #v(1.1cm)
]

#align(center)[
  #box(width: 78%)[
    #text(size: 9.8pt, style: "italic")[
      神经网络开始替代数值求解器之后，出现一个新问题：网络给出的预测
      要不要信、信多少。这份笔记沿这条线索走：先把“不确定性”本身拆开
      （它有几种来源、怎么度量、怎么校准），再进入方法主线，
      信息瓶颈如何同时给出预测均值与可信的不确定度，以及它的四个演化方向
      （采样、降维理论、置换不变算子、物理约束的潜变量模型）。
    ]
  ]
]

#v(0.8cm)

#boxed(title: "怎么读这份笔记")[
  - *概念准备* 备齐起点：不确定性从哪来（数据、模型、网络三层的分解），
    贝叶斯模型平均与 MC 采样怎么把参数不确定性传播成预测分布，
    以及“校准”这个容易被忽略的维度。
  - *综述篇* 到 *评测协议实操* 四章是全景：方法分类、经典方法的逐个精讲、
    指标体系与大规模比较的结论。读它们建立坐标。
  - *方法核心* 与 *IB-UQ 深入* 两章是主线：信息瓶颈怎么被改造成 UQ 工具。
    变分下界的推导完整保留；它与信息瓶颈谱系（Tishby、VIB）的关系逐条对照。
  - *后继一* 到 *后继四* 是四个演化方向：扩散生成器采样 Boltzmann 分布（EDG）、
    反应坐标的最优降维（FMRC）、置换不变的算子学习 UQ（UQ-SONet）、
    物理约束的潜变量模型（LVM-GP）。每个后继后面紧跟一章“实验台／实现台”，
    把论文正文与附录里的算法、数据配置、超参数和全部数值摊开核对。
  - *谱系* 与 *坐标* 两章把这条线放回更大的文献网络，并说明它现在停在何处。
    *附录* 是符号表与复现清单。
]

#v(0.5cm)

#boxed(title: "术语地图", fill: rgb("#f7f7fb"), stroke: rgb("#b9b9d6"))[
  #table(
    columns: (auto, 1fr),
    inset: (x: 6pt, y: 3pt),
    align: (left, left),
    stroke: none,
    [*aleatoric / epistemic*], [数据噪声带来的不可约不确定性 / 有限知识带来的可约不确定性。],
    [*后验分布*], [看到数据之后参数（或解）的条件分布；UQ 的中心对象。],
    [*贝叶斯模型平均*], [按后验加权平均各参数取值下的预测，得到预测分布。],
    [*校准*], [预测的置信区间与真实误差频率一致：95% 区间盖住约 95% 真值。],
    [*过置信*], [区间太窄：声称的置信度高于实际覆盖率。深度网络最常见。],
    [*OOD*], [分布外（out-of-distribution）：输入落在训练数据没覆盖的区域。],
    [*信息瓶颈*], [在压缩输入与保留目标信息之间取权衡的表示学习框架。],
    [*变分下界*], [把算不出的量变成可优化的上/下界，深度生成模型的标准工具。],
    [*normalizing flow / GIN*], [可逆变换链；保体积版本（GIN）用来给输入分布“加温”。],
    [*配分函数*], [Boltzmann 分布的归一化常数，积分通常算不出，采样方法必须绕开它。],
    [*反应坐标*], [高维动力系统的少数关键变量，长期行为由它们决定。],
    [*置换不变*], [输入集合换顺序输出不变；处理“变长传感器集合”的性质。],
    [*潜变量模型*], [用隐随机变量解释观测相关性的模型；VAE 是其神经版本。],
    [*高斯过程*], [函数上的高斯先验，由核函数指定任意两点取值的相关性。],
  )
]

#v(0.6cm)
#outline(title: [目录], indent: 1.4em, depth: 2)

#pagebreak()

#include "parts/c00_foundation.typ"
#include "parts/c01_uq_basics.typ"
#include "parts/c02_survey.typ"
#include "parts/c02b_methods.typ"
#include "parts/c02c_eval.typ"
#include "parts/c03_ibuq.typ"
#include "parts/c03b_ibuq_more.typ"
#include "parts/c04_followups.typ"
#include "parts/c04c_edg_lab.typ"
#include "parts/c04b_more.typ"
#include "parts/c04d_fmrc_theory.typ"
#include "parts/c04e_sonet_lab.typ"
#include "parts/c04f_lvmgp_lab.typ"
#include "parts/c05_related.typ"
#include "parts/c06_appendix.typ"
