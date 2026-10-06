#import "../common.typ": *

= 从拉格朗日量到可运行的算法

理论到这里已经完整。这一节把它翻译成可以反向传播的损失函数，并给出两份算法。

== 统一拉格朗日量

#definition("定义 15（信息论拉格朗日量，论文式 8）")[
  $ max_z quad alpha log I(z_(t-n); z_t) - beta I(z_t; x_t | z_(t-n)) + gamma S( cal(C)/op(tr)(cal(C)) ) + log p(x_t | z_t) , $

  其中 $alpha, beta, gamma$ 是拉格朗日乘子。
]

四项的来历，逐项对照前面三个命题：

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  align: (left, left, left),
  stroke: 0.4pt + luma(190),
  table.header([*项*], [*管什么*], [*来自哪个命题*]),
  [$alpha log I(z_(t-n); z_t)$], [保留时间相干信息（长期预测能力的唯一来源）], [命题 3、2 的下界],
  [$-beta I(z_t; x_t | z_(t-n))$], [压缩快速耗散 + 残差（可压缩的成分）], [命题 3],
  [$gamma S(cal(C)\/op(tr)(cal(C)))$], [维持谱多样性，防模式坍缩], [命题 4、5],
  [$log p(x_t|z_t)$], [重构（解码器）], [基础自监督目标],
)

注意第二项里的一个恒等式，它说明一项惩罚覆盖了两个成分：

$ I(z_t; x_t | z_(t-n)) = I(z_t; x_(t-1)|z_(t-n)) + I(z_t; x_t | x_(t-1)) . $

也就是说，*一个 $beta$ 同时惩罚了快速耗散信息与残差信息*，这两者正是命题 3 里
判定为“可压缩”的两类。

#warn("这一项是“最大化”的形式，容易看错符号")[
  拉格朗日量是 $max$，而第二项带负号，所以它是在*最小化*条件互信息。这符合直觉：
  耗散和残差成分不该被保留。

  但落到代码里时，最小化形式会让符号反过来（见下面的 Algorithm 1）。
  阅读论文时请始终盯住它是在 $max$ 还是在 $min$。
]

=== 为什么不能直接算

$I(z_(t-n); z_t)$ 有个紧凑的闭式（命题 3a），但它在 $d$ 很大时行列式代价高，
而且当潜在是高维时数值不稳定。$I(z_t; x_t|z_(t-n))$ 更是需要条件密度。
$S$ 需要每步做一次特征分解。

所以论文给出一组*可替代的、逐项的*实现。

== 可训练损失

#definition("定义 16（可训练损失，论文式 9）")[
  $ max sum_n [ alpha I(z_n; cal(P)_n) + beta EE_(p_theta(z_n|x_n))[log q_psi(z_n|z_(n-1))] + beta H_(p_theta)(z_n | x_n) + log p_omega(x_n|z_n) ] + gamma S( cal(C)/op(tr)(cal(C)) ) + cal(L)_"ELBO" . $
]

下面逐项解释每一项是什么、为什么是它。

=== 第 1 项：时间相干（InfoNCE）

$ I(z_n; cal(P)_n) $ 把 $z_n$ 与它的*时间邻域* $cal(P)_n = {z_(n plus.minus i) | 1 <= i <= k}$
联系起来。实现上用 InfoNCE：

$ I(z_n; cal(P)_n) approx 1/abs(cal(P)_n) sum_(p in cal(P)_n) log (exp(z_n^T z_p \/ tau))/(sum_(j=1)^B exp(z_n^T z_j \/ tau)) . $

#insight("为什么用邻域而不是只用相邻两步")[
  1. 用 $k > 1$ 的邻域相当于把 $n$ 取大，专门去捕捉 $abs(lambda) approx 1$ 的相干成分，
     同时把只有一步相关性的瞬时噪声压掉。
  2. InfoNCE 是互信息的下界（§2），所以这一项是保守地优化时间相干。
  3. 低维潜在时也可以用命题 3a 的闭式直接算，二者等价但代价不同。
]

=== 第 2、3 项：结构一致性

这是最新颖的一项，也是把“信息论”和“Koopman 线性”扣在一起的那颗螺钉。

#proposition("结构一致性的等价形式（论文式 50–51）")[
  给定 $z_(n-1)$ 时，条件互信息为

  $ I(z_n; x_n | z_(n-1)) = H(z_n | z_(n-1)) - H(z_n | x_n, z_(n-1)) . $

  由于编码器与 $z_(n-1)$ 无关（Fact F.1），$H(z_n|x_n,z_(n-1)) = H_(p_theta)(z_n|x_n)$，于是

  $ I(z_n; x_n | z_(n-1)) = H(z_n|z_(n-1)) - H_(p_theta)(z_n | x_n) . $

  又按定义

  $ I(z_n; x_n | z_(n-1)) = EE_(p(z_n|x_n)) [ log (p(z_n|x_n))/(q(z_n|z_(n-1))) ] = EE_(p(z_n|x_n))[ -log q(z_n|z_(n-1)) ] - H_(p_theta)(z_n|x_n) . $
]

#proof[
  *第一式。* 条件互信息的定义 $I(a;b|c) = H(a|c) - H(a|b,c)$（第 2 章定义 6 之前给过）
  取 $a = z_n$、$b = x_n$、$c = z_(n-1)$，逐符号代入即得。

  *第二式，第一步：写成 KL。* 条件互信息有 KL 形式：
  $I(a;b|c) = EE_(p(b,c)) D_"KL"(p(a|b,c) ‖ q(a|c))$，
  它度量“真实的条件分布 $p(a|b,c)$ 与先验 $q(a|c)$ 差多远”。
  这里 $q$ 是 Koopman 先验 $q_psi(z_n|z_(n-1))$，表达式里本来就不含 $b$；
  又由 Fact F.1，编码器只看 $x_n$，$p(z_n|x_n,z_(n-1)) = p(z_n|x_n)$。两个实参都与 $b$
  无关，对 $b$ 的期望积掉（内层是常数）：

  $ I(z_n; x_n | z_(n-1)) = EE_(p(z_n|x_n)) [ log (p(z_n|x_n))/(q_psi(z_n|z_(n-1))) ] . $

  *第二式，第二步：拆对数。* 把对数拆成
  $log p - log q = (-log q) - (-log p)$，两项分别取期望：第二项是
  $EE[-log p_theta(z_n|x_n)] = H_(p_theta)(z_n|x_n)$（条件熵的定义），于是

  $ I(z_n; x_n | z_(n-1)) = EE_(p(z_n|x_n))[ -log q_psi(z_n|z_(n-1)) ] - H_(p_theta)(z_n|x_n) . quad square $
]

这两步的结论值得单独列出：

#keypoint("结构一致性为什么是一个互信息")[
  $ I(z_n; x_n | z_(n-1)) = EE_(p(z_n|x_n))[ - log q_psi(z_n|z_(n-1)) ] - H_(p_theta)(z_n|x_n) . $

  左端正是拉格朗日量第二项要*最小化*的东西。因此最小化它等价于在代码里*最大化*

  $ EE_(p_theta(z_n|x_n))[ log q_psi(z_n|z_(n-1)) ] + H_(p_theta)(z_n|x_n) . $

  这就是可训练损失里的第 2、3 项。
]

两个效应各自做什么：

- *与 Koopman 动力学对齐*：$EE[-log q(z_n|z_(n-1))]$ 要求从编码器采样的点落在
  Koopman 先验的*高似然区域*。因为先验被参数化为*线性高斯转移*
  $q(z_n|z_(n-1)) = cal(N)(z_n | cal(K) z_(n-1), Sigma)$，最小化这个 KL 就迫使编码器输出
  可以在线性结构下被预测。
- *熵正则*：$H_(p_theta)(z_n|x_n)$ 防止编码器变成确定性的
  （否则 $I(x;z)$ 在连续情形下退化，见 §2 的警告）。

合起来，$p(z_n|x_n) -> q(z_n|z_(n-1))$，也就是

$ p(z_n|x_n) approx q(z_n|z_(n-1)) quad <==> quad z_n "在" cal(K) "下近似线性演化" , $

这就是结构一致性。

=== 第 4 项：重构

$ log p_omega(x_n|z_n) $ 是解码器损失。在 AE 版里它退化成 $L_2$：$norm(x_n - g_omega(z_n))^2$。

=== 第 5 项：von Neumann 熵（预测充分性）

从一个小批量里估计潜在协方差：

$ cal(C) = 1/B sum_(i=1)^B (z_i - bar(z))(z_i - bar(z))^T, quad bar(z) = 1/B sum_i z_i, quad P = cal(C)/op(tr)(cal(C)), quad S(P) = - sum_j lambda_j log lambda_j , $

其中 $lambda_j$ 是 $P$ 的第 $j$ 个特征值。这一项直接对应命题 5，防坍缩。

=== 第 6 项：ELBO

$ cal(L)_"ELBO" = log p_omega(x_(n-1)|z_(n-1)) - D_"KL"( p_theta(z_(n-1)|x_(n-1)) \|\| cal(N)(0, I) ) . $

这一项来自 VAE 的标准目标，作用有两个：训练稳定性，以及额外的重构信号。
它是把 IB/VIB 那套变分机制接到自编码器上的接口。

== 两份算法

=== VAE 版（概率）

#boxed(title: "Algorithm 1（论文版，概率）")[
  输入：数据 $cal(D) = {x_n}_(n=0)^T$；拉格朗日乘子 $(alpha, beta, gamma)$；学习率 $eta$；
  轮数 $K$；批量 $B$；邻域窗口 $k$；温度 $tau$。

  + 初始化编码器 $p_theta(z|x)$、解码器 $p_omega(x|z)$、潜在动力网络 $q_psi(z_n|z_(n-1))$。
  + *for* epoch $= 1$ to $K$:
    + *for* 每个小批量 ${x_1, ..., x_B}$:
      + 采样潜在 $z_i tilde p_theta(z|x_i)$。
      + *时间相干*：对每个 $z_n$，把时间邻居 $cal(P)_n = {z_(n plus.minus i) | 1 <= i <= k}$ 当作正样本，算 InfoNCE。
      + *结构一致性*：算 $EE_(p_theta(z_n|x_n))[log q_psi(z_n|z_(n-1))]$。
      + *预测充分性*：算 $cal(C)$、$P = cal(C)\/op(tr)(cal(C))$、$S(P)$。
      + *ELBO*：算 $cal(L)_"ELBO"$。
      + 总损失：
        $ cal(L) = -[ alpha I(z_n;cal(P)_n) + EE[log q_psi(z_n|z_(n-1))] + H_(p_theta)(z_n|x_n) + log p_omega(x_n|z_n) + gamma S(P) + cal(L)_"ELBO" ] . $
      + 用 Adam 更新 $theta, omega, psi$。
]

#warn("注意论文里两处的写法略有差别")[
  式 (9) 把 $beta$ 显式写在 $EE[log q]$ 与 $H$ 前面；
  Algorithm 1 的总损失里没有显式写出 $beta$（等价于把它吸收进了那两项的定义）。

  这不是错误，只是同一件事的两种记法。但如果你要复现，必须知道：
  *实际起作用的是三个独立的乘子，一个管 InfoNCE（$alpha$），一个管结构一致性（$beta$），
  一个管熵（$gamma$）。*
]

=== AE 版（确定性）

当编码器是确定性的 $z_n = f_theta(x_n)$ 时，概率量退化：

- $EE_(p(z_n|x_n))[log q_psi(z_n|z_(n-1))]$ 退化成 $L_2$ 一致性损失
  $cal(L)_"Koop" = norm(z_(n+1) - cal(K)_psi z_n)^2$；
- $cal(L)_"ELBO"$ 退化成 AE 重构项。

#boxed(title: "Algorithm 2（论文版，确定性）")[
  + 初始化确定性编码器 $z_n = f_theta(x_n)$、解码器 $hat(x)_n = g_omega(z_n)$、Koopman 算子 $cal(K)_psi$。
  + *for* epoch、*for* 小批量：
    + 编码 $z_i = f_theta(x_i)$。
    + 时间相干：同一个 InfoNCE。
    + 结构一致性（确定性）：$cal(L)_"Koop" = norm(z_(n+1) - cal(K)_psi z_n)^2$。
    + 预测充分性：从 $cal(C)$ 算 $S(P)$。
    + 重构：$cal(L)_"rec" = norm(x_n - g_omega(z_n))^2$。
    + 总损失：
      $ cal(L) = cal(L)_"rec" - alpha I(z_n; cal(P)_n) + beta cal(L)_"Koop" - gamma S(P) . $
    + 用 Adam 更新。
]

#insight("AE 版更清楚，适合用来理解三项的作用")[
  $cal(L) = underbrace(cal(L)_"rec", "重构") - underbrace(alpha I(z_n;cal(P)_n), "时间相干（最大化）") + underbrace(beta cal(L)_"Koop", "结构一致性（惩罚）") - underbrace(gamma S(P), "谱多样性（最大化）") . $

  四个符号：重构是“$+$”（要小），时间相干是“$-$”（要*大*，所以最小化时带负号），
  结构一致性是“$+$”（要小），熵是“$-$”（要*大*）。
  这份最简洁的式子把整篇论文讲完了。
]

== 架构细节

潜在空间里的线性算子在三种任务上有三种形态（论文 Table 5）：

#table(
  columns: (auto, auto, auto),
  inset: 6pt,
  align: (left, left, left),
  stroke: 0.4pt + luma(190),
  table.header([*环境*], [*结构*], [*关键特征*]),
  [物理仿真], [AE / VAE], [$z_(t+1) = cal(K) z_t$；重构；Koopman 线性前向；InfoNCE；von Neumann 熵],
  [视觉输入（控制）], [AE], [$z_(t+1) = cal(K) z_t + B a_t + epsilon$（线性高斯）；VAE ELBO；InfoNCE；VNE],
  [图结构动力学], [AE], [$z_(t+1) = cal(K)(A) z_t$（邻接条件化）；重构；InfoNCE；VNE],
)

三个细节：

- *控制任务*里线性算子是 $cal(K) z_t + B a_t$。这是控制理论里的标准线性时不变形式，
  它能让“在潜在空间里做 LQR”变成一件合法的事，而 LQR 只有在潜在动力学线性时才成立。
  这是 Koopman 表示相对一般表示学习的*实际好处*，不只是美学上的。
- *图任务*里 $cal(K)(A) := A times cal(K)$，即由邻接矩阵 $A$ 调制的共享 Koopman 复合
  （沿用 CKO 的设定）。这样参数量与节点数无关，可以外推到更多节点。
- *架构主体*：物理仿真用卷积编码器（$3$ 层，逐步下采样到 $H\/16 times W\/16$）$->$
  展平 $->$ 全连接到 $d_s$ 维潜在 $->$ 单层线性 $cal(K)$ $->$ 对称解码器。
  ERA5 任务用 factorized attention 编码器。视觉输入与图任务沿用 E2C 与 CKO 的原始设置
  以保证公平比较。

超参数（论文 Table 12）：

#table(
  columns: (auto, auto, auto, auto, auto),
  inset: 6pt,
  align: (left, center, center, center, center),
  stroke: 0.4pt + luma(190),
  table.header([*含义*], [*符号*], [*物理仿真*], [*视觉输入*], [*图动力学*]),
  [时间相干权重], [$alpha$], [2.00], [3.00], [2.00],
  [结构一致性权重], [$beta$], [--], [2.00], [--],
  [von Neumann 熵权重], [$gamma$], [0.10], [0.50], [0.10],
  [InfoNCE 邻域], [$k$], [3], [5], [5],
)

#insight("从这张表能读出两件事")[
  1. *物理仿真任务里 $beta$ 是空的。* 因为那里用的是 AE 结构，结构一致性已经是
     一个硬性的 $L_2$ 项（$cal(L)_"Koop"$），不需要额外的概率权重。
     视觉与控制任务用 VAE 结构，才需要 $beta$ 来平衡 KL。

  2. *$gamma$ 的取值差了一个量级*（0.10 vs 0.50）。这符合命题 5 的分析：
     熵正则的合适强度取决于潜在协方差的谱有多倾斜。控制任务的高维图像输入
     更容易坍缩，所以需要更强的熵正则。
]

#boxed(title: "计算开销（论文 Table 7）", fill: rgb("#f2fbf5"), stroke: rgb("#8dc9a3"))[
  加上 InfoNCE 与 von Neumann 熵后，每个 epoch 的额外开销很小：

  - Kármán 涡：总 epoch 时间 $201.23$ s，其中 InfoNCE $13.78$ s、熵 $0.97$ s；
  - Dam Flow：$18.92$ s，其中 $0.76$ s / $0.56$ s；
  - ERA5：$253.24$ s，其中 $15.09$ s / $7.43$ s。

  也就是说，理论换来的收益额外付出的计算代价可以忽略。
]

#pagebreak()
