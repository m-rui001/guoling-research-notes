#import "../../preamble.typ": *

= 工具：稀疏恢复的两把尺子

== RIP：矩阵级判据

#definition("定义（RIP，有限等距性质）")[
  矩阵 $bold(D) in RR^(M times N)$ 满足阶数 $s$、常数 $delta_s$ 的 RIP，
  如果对所有 $s$-稀疏向量 $bold(v)$，

  $ (1 - delta_s) norm(bold(v))_2^2 <= norm(bold(D) bold(v))_2^2 <= (1 + delta_s) norm(bold(v))_2^2 . $

  即 $bold(D)$ 近似保持所有稀疏向量的长度。$delta_s$ 越小，矩阵越“等距”。
]

RIP 下的恢复定理（Candès 型）如下。

#proposition([定理（$ell_1$ 恢复，RIP 版）])[
  设 $bold(D)$ 的 RIC（受限等距常数）满足 $delta_s < 0.307$。
  对任意 $tilde(bold(c)) in RR^N$，$bold(c)^sharp$ 是基追踪

  $ arg min norm(bold(c))_1 quad "s.t." quad bold(D) bold(c) = bold(D) tilde(bold(c)) $

  的解，则

  $ norm(bold(c)^sharp - tilde(bold(c)))_2 <= C (sigma_(s, 1)(tilde(bold(c))))/(sqrt(s)) , $

  其中 $sigma_(s, 1)(tilde(bold(c)))$ 是 $tilde(bold(c))$ 的最佳 $s$-项逼近误差
  （把最小的 $N - s$ 个系数抹掉后的 $ell_1$ 尾巴）。
  特别地，$tilde(bold(c))$ 恰好 $s$-稀疏时恢复精确。
]

读法：恢复误差由“去掉 $s$ 个最大系数后剩下的尾巴”控制。
目标越稀疏（尾巴越小），恢复越准。整个压缩感知游戏的玩法因此清晰：
设计采样测度，让由它导出的 $bold(D)$ 以高概率满足小 $delta_s$。

样本量代价：随机采样下 RIP 不是白来的。

#proposition("定理（随机采样的 RIP，Rauhut–Ward）")[
  点从“正交化测度”$nu$ 独立采样（含义见下节），
  若样本量满足

  $ M >= C delta^(-2) L s log^3 (s) log (N) , $

  则以概率至少 $1 - N^(-gamma log^3 (s))$，矩阵 $bold(D)\/sqrt(M)$ 的
  $delta_s <= delta$。
]

这个形式的每一块都有具体含义：$s$ 线性项（稀疏度，省不掉）；
$delta^(-2)$（要求越严越贵）；两个对数因子（log 项是这类证明的技术代价）；
以及 $L$，*一个由采样测度与多项式系统决定的常数*。六篇论文的第一篇
就是要把它在各类分布下显式估计出来，$L$ 越小，需要的样本越少。

== MIP：列级判据

第二个判据更直观，作用在矩阵的列上。

#definition("定义（MIC，互不相干常数）")[
  $ mu = mu(bold(Phi)) = max_(k != j) abs(⟨bold(Phi)_k, bold(Phi)_j ⟩) / (norm(bold(Phi)_k)_2 norm(bold(Phi)_j)_2) , $

  即任意两列归一化内积的绝对值的最大值。
  经典结论（Donoho–Huo、Fuchs、Cai–Wang–Xu）：若 $mu < 1\/(2s - 1)$，
  则 $s$-稀疏向量被基追踪精确恢复；带噪情形稳定。
]

RIP 与 MIP 的分工：RIP 给出“接近最优”的整体保证，样本复杂度好；
MIP 条件强（$1\/(2s-1)$ 难达到），但它可以用*逐列的上界*
（不依赖样本随机性，只看基函数与采样测度）直接验证。
梯度增强那篇（第 9 章）走的就是 MIP 路线，因为它的测量矩阵是拼出来的，
列的逐点结构可以手工控制。

两把尺子都指向同一个问题：*采样点从哪个分布抽，矩阵才“好”？*
下一节开始，六篇论文按这个问题的六个答案依次展开。
