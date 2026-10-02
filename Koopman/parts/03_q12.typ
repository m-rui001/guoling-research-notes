#import "../common.typ": *

= 论文的论证链（上）：信息会丢吗？丢多少？

到这里两条线索汇合了。Koopman 一侧的问题（表示欠定、模式坍缩、缺乏一般原则）
遇上信息论一侧的工具（IB、变分界、InfoNCE、von Neumann 熵）。

论文的论证链是五个依次递进的提问。这个顺序不是我加的，是论文自己在附录 B 里明确写出的：

#keypoint("论文自己的提问链")[
  1. 信息会丢吗？（命题 1）
  2. 如果会，丢多少？（命题 2）
  3. 丢的是哪一类信息？（命题 3）
  4. 如何优化信息保留？（命题 4）
  5. 如何避免副作用（模式坍缩）？（命题 5）
]

命题块用蓝色边框标出；证明块给出完整推导，不跳步。建议先读命题再看证明。

== 先把信息流的结构写清楚

定义记号：$x_(1:t)$ 是 1 到 $t$ 时刻的状态，$z_(1:t)$ 是自回归生成的潜在变量。
由两张表的对照（论文 Table 1），Koopman 表示诱导出的轨迹分布是

$ p_"KR"(x_(1:t) | x_0) = integral p(z_0|x_0) product_(n=1)^t p(z_n|z_(n-1)) p(x_n|z_n) d z_0 d z_1 dots d z_t, $

三个因子各自对应架构的一个部件，必须认清楚：

- $p(z_0|x_0)$ 是*编码器*，把初始状态映进潜在空间；
- $p(z_n|z_(n-1)) = cal(N)(z_n | cal(K) z_(n-1), Sigma)$ 是*潜在前向*，
  它是 Koopman 算子定义的式 1 的概率版本（$Sigma$ 是噪声方差）；
- $p(x_n|z_n)$ 是*解码器*，通常是高斯的。

#figure(
  align(center, canvas({
    import draw: *
    let c1 = rgb("#1a4f9a")
    let c2 = rgb("#b8860b")
    // 状态行
    content((0,0), text(fill: c1, size: 9.5pt)[$x_(t-n)$])
    content((2.9,0), text(fill: c1, size: 9.5pt)[$x_(t-n+1)$])
    content((5.6,0), text(fill: c1, size: 9.5pt)[$dots.c$])
    content((7.0,0), text(fill: c1, size: 9.5pt)[$x_t$])
    line((0.75,0), (2.45,0), stroke: 1.1pt + c1, mark: (end: ">", fill: c1))
    line((3.4,0), (5.3,0), stroke: 1.1pt + c1, mark: (end: ">", fill: c1))
    line((5.95,0), (6.7,0), stroke: 1.1pt + c1, mark: (end: ">", fill: c1))
    content((1.6,0.35), text(size: 8pt, fill: c1)[$T$])
    content((4.35,0.35), text(size: 8pt, fill: c1)[$T$])
    content((9.9,0), text(size: 8pt, fill: c1)[真实非线性动力学])
    // 潜在行
    content((0,-1.9), text(fill: c2, size: 9.5pt)[$z_(t-n)$])
    content((2.9,-1.9), text(fill: c2, size: 9.5pt)[$z_(t-n+1)$])
    content((5.6,-1.9), text(fill: c2, size: 9.5pt)[$dots.c$])
    content((7.0,-1.9), text(fill: c2, size: 9.5pt)[$z_t$])
    line((0.75,-1.9), (2.45,-1.9), stroke: 1.1pt + c2, mark: (end: ">", fill: c2), dash: "densely-dashed")
    line((3.4,-1.9), (5.3,-1.9), stroke: 1.1pt + c2, mark: (end: ">", fill: c2), dash: "densely-dashed")
    line((5.95,-1.9), (6.7,-1.9), stroke: 1.1pt + c2, mark: (end: ">", fill: c2), dash: "densely-dashed")
    content((1.6,-1.55), text(size: 8pt, fill: c2)[$cal(K)$])
    content((4.35,-1.55), text(size: 8pt, fill: c2)[$cal(K)$])
    content((9.9,-1.9), text(size: 8pt, fill: c2)[线性 Koopman 推进])
    // 编解码
    line((0.0,-0.42), (0.0,-1.48), stroke: 0.9pt + luma(120), mark: (end: ">", fill: luma(120)))
    line((2.9,-0.42), (2.9,-1.48), stroke: 0.9pt + luma(120), mark: (end: ">", fill: luma(120)))
    line((7.0,-1.48), (7.0,-0.42), stroke: 0.9pt + luma(120), mark: (end: ">", fill: luma(120)))
    content((0.55,-0.95), text(size: 7.5pt, fill: luma(100))[(1) 编码])
    content((3.45,-0.95), text(size: 7.5pt, fill: luma(100))[(1) 编码])
    content((7.6,-0.95), text(size: 7.5pt, fill: luma(100))[(3) 解码])
    // 损失标注
    content((4.9,-2.9), cnode(10.4cm, 0.75cm, rgb("#1a7a45"), [(2) 结构一致性：$z$ 必须沿 $cal(K)$ 线性演化 ， 这是 Koopman 表示比普通 VAE 严格之处], size: 8pt))
  }), ),
  caption: [Koopman 表示的自回归信息流结构（对应论文附录 Figure 8）。上排实线是真实的非线性演化 $T$；
  下排虚线是潜在空间里由 Koopman 算子 $cal(K)$ 支配的线性演化。信息可能在三个位置流失：
  (1) 编码 $x -> z$；(2) 线性推进 $z -> z$；(3) 解码 $z -> x$。整篇论文的理论绕这三个流失点展开。],
)

#keypoint("这张表定义的对照关系")[
  #ttable(
    ([目标], [一般表示学习], [Koopman 表示]),
    (
      ([信息流], [$x -> z -> y$], [$x_(n-1) -> z_(n-1) ->^"Koopman" z_n -> x_n$]),
      ([拉格朗日], [$I(x;z) - beta I(z;y)$], [$I(x_(n-1); z_(n-1)) - beta I(z_(n-1); x_n)$]),
    ),
    widths: (auto, auto, auto),
  )

  右列多了一节：$z$ 必须支持到 $z'$ 的*线性*推进。这就是 Koopman 表示比普通表示学习严格的地方。
]

== 命题 1：信息会丢吗？

#proposition("命题 1（潜在演化中的信息丢失）")[
  设 $x_(n-1) -> z_(n-1) ->^cal(K) z_n -> x_n$ 表示 Koopman 表示中的信息传播，则

  $ I(x_(n-1); x_n) >= I(z_(n-1); x_n) >= I(z_(n-1); z_n). $
]

先建立两条基础事实。

#lemma("Fact F.1（编码器只看当前状态）")[
  设 $z_t$ 由只依赖 $x_t$ 的概率编码器 $p(z_t|x_t)$ 给出。于是

  (a) $H(z_t) <= H(x_t)$（数据处理不等式）；

  (b) $I(z_t; xi | x_t) = 0$，即 $z_t perp xi | x_t$（给定 $x_t$，$z_t$ 与系统中任何变量 $xi$ 条件独立）。
]

(b) 是后面所有展开中“修正项归零”的来源。直觉：编码器是 $x_t$ 的函数，所以给定 $x_t$，
$z_t$ 不再携带关于其它任何东西的额外信息。

#proof[
  *第一步：证明 $I(x_(t-1); x_t) >= I(z_(t-1); x_t)$。*

  用互信息的链式法则（引理 1）把 $I(x_(t-1);x_t)$ 展开。注意

  $ I((x_(t-1), z_(t-1)); x_t) = I(x_(t-1); x_t) + I(z_(t-1); x_t | x_(t-1)) . $

  由 Fact F.1(b)，给定 $x_(t-1)$ 时 $z_(t-1)$ 与 $x_t$ 条件独立，故第二项为零。于是

  $ I(x_(t-1); x_t) = I((x_(t-1), z_(t-1)); x_t) . $

  再把右边的联合变量拆开：

  $ I((x_(t-1), z_(t-1)); x_t) = I(z_(t-1); x_t) + I(x_(t-1); x_t | z_(t-1)) . $

  由互信息的非负性，第二项 $>= 0$，于是

  $ I(x_(t-1); x_t) = I(z_(t-1); x_t) + I(x_(t-1); x_t | z_(t-1)) >= I(z_(t-1); x_t) . quad square $

  *第二步：证明 $I(z_(t-1); x_t) >= I(z_(t-1); z_t)$。*

  同样用链式法则：

  $ I(z_(t-1); (x_t, z_t)) = I(z_(t-1); x_t) + I(z_(t-1); z_t | x_t) . $

  由 Fact F.1(b)（$z_t$ 只依赖 $x_t$），第三项为零，故 $I(z_(t-1); x_t) = I(z_(t-1); (x_t,z_t))$。再拆一次：

  $ I(z_(t-1); (x_t, z_t)) = I(z_(t-1); z_t) + I(z_(t-1); x_t | z_t) >= I(z_(t-1); z_t) . quad square $
]

把两步合起来，我们还免费得到信息缺口的*精确*表达式（这是后面要用的）：

$ I(x_(t-1); x_t) - I(z_(t-1); z_t) = underbrace(I(z_(t-1); x_t | z_t), "演化中损失") + underbrace(I(x_(t-1); x_t | z_(t-1)), "编码中损失") . $

两条不等式的物理含义值得逐字读：

#insight("两条不等式的含义")[
  *第一条 $I(x_(t-1);x_t) >= I(z_(t-1);x_t)$*：映射 $x_(t-1) -> z_(t-1)$ 是有损压缩，
  它可能丢掉对未来状态 $x_t$ 有用的信息。这是*编码器的责任*。

  *第二条 $I(z_(t-1);x_t) >= I(z_(t-1);z_t)$*：潜在前向由 Koopman 算子 $cal(K)$ 支配，
  它*本质上限制*了潜在空间能保留的信息。注意这里比较的两个量不对称，
  $x_t$ 比 $z_t$ 携带更多对 $z_(t-1)$ 的依赖，所以中间项通常更大。
  $I(z_(t-1);z_t)$ 就是 $cal(K)$ 给 Koopman 表示设置的*信息上限*。
]

=== 多步版本

单步不够用，自回归预测是几十步展开的。多步版本可以递归地得到。

#proof[
  先从 $I(x_(t-n); x_t)$ 出发，完全重复单步的过程：

  $ I(x_(t-n); x_t) = I((x_(t-n), z_(t-n)); x_t) = I(z_(t-n); x_t) + I(x_(t-n); x_t | z_(t-n)) = I(z_(t-n); z_t) + I(z_(t-n); x_t|z_t) + I(x_(t-n); x_t | z_(t-n)) . $

  接下来处理最后一项 $I(x_(t-n); x_t | z_(t-n))$。引入 $z_(t-n+1)$：

  $ I(x_(t-n); x_t | z_(t-n)) = I(x_(t-n); z_(t-n+1) | z_(t-n)) + I(x_(t-n); x_t | z_(t-n), z_(t-n+1)) . $

  这里之所以没有修正项，是因为被积掉的修正项 $I(x_(t-n); z_(t-n+1) | x_t, z_(t-n))$
  为零：由 Fact F.1(b)，$z_(t-n+1)$ 只依赖 $x_t$，给定 $x_t$ 后它与 $x_(t-n)$ 无关。
  再引入 $z_(t-n+2)$，重复这一操作。记 $z_(t-n:i) := (z_(t-n), ..., z_i)$，对
  $i = t-n+1, ..., t$ 逐次递归，得到

  $ I(x_(t-n); x_t | z_(t-n)) = sum_(i=t-n+1)^t I(x_(t-n); z_i | z_(t-n:i-1)) + I(x_(t-n); x_t | z_(t-n:t)) . $

  合并，多步信息缺口为

  $ I(x_(t-n); x_t) - I(z_(t-n); z_t) = sum_(i=t-n+1)^t underbrace(I(x_(t-n); z_i | z_(t-n:i-1)), "快速耗散") + underbrace(I(x_(t-n); x_t | z_(t-n:t)), "残差") + underbrace(I(z_(t-n); x_t | z_t), "演化损失") . quad square $
]

这三项的含义在命题 3 里会被重新捡起来，先记住它们：

- $I(x_(t-n); z_i | z_(t-n:i-1))$：在已知此前所有潜在状态的条件下，新的潜在变量 $z_i$
  还能揭示多少关于过去状态 $x_(t-n)$ 的信息。这是*快速耗散*的部分，每多走一步它就迅速衰减。
- $I(x_(t-n); x_t | z_(t-n:t))$：由于压缩表示而无法被潜在序列完全表达的残差依赖。
- $I(z_(t-n); x_t | z_t)$：$z_(t-n)$ 中关于 $x_t$ 的信息量减去 $z_t$ 保留的部分。
  这一项为正说明*潜在演化过程中信息丢失*，对应 $abs(lambda) < 1$ 的 Koopman 模式。

== 命题 2：丢多少？

命题 1 是抽象的，它只说“信息会丢”，没给出能算的量。要把它变成能上界的量，
需要三样工具，这里先把它们各自的作用说清楚：

- *一把衡量“两个分布差多远”的尺子*。因为我们要比较的是“真实轨迹分布”
  与“Koopman 模型产生的轨迹分布”。下面用 KL 散度与全变差距离。
- *Pinsker 不等式*：它把“KL 散度”换算成“全变差距离”，也就是把抽象的分布差异
  换成更接近直觉的量，两个分布有多少概率质量对不上。
- *率失真理论*（定义见 §2）：把“信息量”换算成“误差的量级”。
  这一步才让结论从“信息会丢”升级成“误差有多大”。

先建立一个变分分解。

#lemma("引理 3（变分分解，论文式 18）")[
  设 $p$ 是真实轨迹分布，$p_"KR"$ 是式 (2) 中的理想 Koopman 模型分布，$q_"KR"$ 是它的变分近似。则

  $ D_"KL"(p(x_(1:t)|x_0) ‖ q_"KR"(x_(1:t)|x_0)) <= D_"KL"(p(x_(1:t)|x_0) ‖ p_"KR"(x_(1:t)|x_0)) + cal(E)_"enc" + cal(E)_"tra" + cal(E)_"rec" . $
]

#proof[
  先拆出理想模型与真实分布的差：

  $ D_"KL"(p ‖ q_"KR") = EE log (p(x_(1:t)|x_0))/(p_"KR"(x_(1:t)|x_0)) + EE log (p_"KR"(x_(1:t)|x_0))/(q_"KR"(x_(1:t)|x_0)) . $

  第一项就是 $D_"KL"(p ‖ p_"KR")$。第二项把 $p_"KR", q_"KR"$ 都按式 (2) 展开成
  编码器、转移、解码器的乘积之积。由“边缘化不会增加 KL 散度”（边缘化指把联合分布对某些变量积分掉、得到边缘分布；这一步丢弃信息，因此两个分布之间的差异只会变小。原式 $D_"KL"(p_X ‖ q_X) <= D_"KL"(p_(X Y) ‖ q_(X Y))$），

  $ EE log (p_"KR"(z_(0:t), x_(1:t)|x_0))/(q_"KR"(z_(0:t), x_(1:t)|x_0)) = underbrace(D_"KL"(p(z_0|x_0) ‖ q_"KR"(z_0|x_0)), cal(E)_"enc") + sum_(n=1)^t underbrace(D_"KL"(p(z_n|z_(n-1)) ‖ q_"KR"(z_n|z_(n-1))), cal(E)_"tra") + sum_(n=1)^t underbrace(D_"KL"(p(x_n|z_n) ‖ q_"KR"(x_n|z_n)), cal(E)_"rec") . quad square $
]

三个误差项分别对应三个架构部件：$cal(E)_"enc"$（潜在表示误差）、$cal(E)_"tra"$（Koopman 算子误差）、
$cal(E)_"rec"$（重构误差）。下面记 $cal(E) := cal(E)_"enc" + cal(E)_"tra" + cal(E)_"rec"$。

#lemma("引理 4（Pinsker 不等式）")[
  对同一空间上的两个分布 $p, q$，全变差距离

  $ norm(p - q)_"TV" := sup_X abs(p(X) - q(X)) = 1/2 integral abs(p(x) - q(x)) d x, $

  且它被 KL 散度控制：

  $ norm(p - q)_"TV" <= sqrt(1/2 D_"KL"(p ‖ q)) . $
]

#proposition("命题 2（Koopman 表示的自回归误差上界）")[
  真实轨迹分布与 Koopman 诱导轨迹分布之间的分布差异被信息缺口控制：

  $ norm(p(x_(1:t) | x_0) - q_"KR"(x_(1:t) | x_0))_"TV" <= 1/2 [ D_"KL"(p(x_(1:t)|x_0) ‖ p_"KR"(x_(1:t)|x_0)) + cal(E) ] <= 1/2 [ sum_(n=1)^t ( I(x_(n-1); x_n) - I(z_(n-1); z_n) ) + cal(E) ] . $

  并且误差上界为

  $ norm(EE_(q_"KR")[x_(1:t) | x_0] - EE_p[x_(1:t) | x_0])_2 <= C_0 sqrt( 2 sum_(n=1)^t (I(x_(n-1);x_n) - I(z_(n-1);z_n)) + cal(E) ) , $

  其中 $C_0$ 是正常数。
]

#proof[
  *第一步：Pinsker 加变分分解。* 直接由引理 4 与引理 3：

  $ norm(p(x_(1:t)|x_0) - q_"KR"(x_(1:t)|x_0))_"TV" <= sqrt(1/2 (D_"KL"(p ‖ p_"KR") + cal(E)_"enc" + cal(E)_"tra" + cal(E)_"rec")) . $

  *第二步：把 KL 与信息缺口联系起来。* 这是关键一步。先看单步：

  $ I(x_(t-1);x_t) - I(z_(t-1);z_t) = EE[log (p(x_t|x_(t-1)))/(p(x_t))] - EE[log (p(z_t|z_(t-1)))/(p(z_t))] . $

  合并两个对数：

  $ = EE[log (p(x_t|x_(t-1)) p(z_t))/(p(x_t) p(z_t|z_(t-1)))] . $

  由贝叶斯法则 $p(x_t, z_t) = p(x_t) p(z_t|x_t) = p(z_t) p(x_t|z_t)$，故
  $p(z_t) \/ p(x_t) = p(z_t|x_t) \/ p(x_t|z_t)$。代入：

  $ I(x_(t-1);x_t) - I(z_(t-1);z_t) = EE[log ((p(x_t|x_(t-1)))/(p(z_t|z_(t-1))) dot (p(z_t|x_t))/(p(x_t|z_t)))] . $

  对 $n = 1, ..., t$ 求和，左右两侧的因子逐项配对并 telescopes：

  $ sum_(n=1)^t ( I(x_(n-1);x_n) - I(z_(n-1);z_n) ) = EE[log (p(z_(0:t), x_(1:t)))/(p_"KR"(z_(0:t), x_(1:t)))] = D_"KL"(p(x_(1:t)|x_0) ‖ p_"KR"(x_(1:t)|x_0)) . $

  代回第一步即得全变差上界。

  *第三步：从全变差到 $L^2$。* 由三角形不等式与全变差的定义，

  $ norm(EE_(q_"KR")[x_(1:t)|x_0] - EE_p[x_(1:t)|x_0])_2 = norm(integral x_(1:t) (q_"KR" - p) d x_(1:t))_2 <= integral norm(x_(1:t)) abs(q_"KR" - p) d x_(1:t) <= 2 C_0 norm(p - q_"KR")_"TV" . $

  其中用到状态位于紧空间 $M$ 上（紧 = 有界且闭，因此状态变量取不到无穷大），故 $norm(x_(1:t)) <= C_0 < infinity$。
  再代入 Pinsker 与第二步的等式，得到式 (6)。$square$
]

命题 2 的读法：

#insight("信息缺口与误差的关系")[
  不等式右侧是 $sum_n [I(x_(n-1);x_n) - I(z_(n-1);z_n)]$。

  - $I(x_(n-1);x_n)$ 度量*原始系统*中非线性映射 $T$ 的内在动力学耦合强度；
  - $I(z_(n-1);z_n)$ 度量这种耦合在 Koopman 表示下还剩多少；
  - 两者之差就是“用线性 Koopman 表示逼近非线性动力学时丢掉的信息”。

  预测误差被这个逐项信息缺口控制。*想减小误差，就要缩小这个缺口*。
  而 $I(z_(n-1);z_n)$ 由命题 1 是信息上限，所以缺口的缩小只能靠提升 $I(z_(n-1);z_n)$，
  这恰好就是拉格朗日量第一项的由来。
]

$abs(lambda) = 1$ 特殊情形（遍历系统）值得一提，因为它把命题 2 变成一句干净的话。

#lemma("Remark F.4（遍历系统）")[
  对遍历系统，条件在 $x_0$ 上取长时间平均：

  $ lim_(t->infinity) 1/t D_"KL"(p(x_(1:t)|x_0) ‖ q_"KR"(x_(1:t)|x_0)) <= lim_(t->infinity) 1/t sum_(n=1)^t ( I(x_(n-1); x_n) - I(z_(n-1); z_n) ) . $
]

左端是相对熵率，右端是每步信息差的平均。也就是说，*在遍历系统上，动态偏差由每步信息差控制*。

=== 下界：误差不能比什么更小

上界告诉我们“做得好有多好”，下界告诉我们“做不做得更好”。

#proposition("命题 2 的对偶（率失真下界）")[
  在 $L^2$ 失真下经理想 Koopman 模型 $p_"KR"$，最小可达失真 $D$ 满足

  $ D >= (n t)/(2 e) exp( 2/(n t) H(x_(1:t)) ) exp( - 2/(n t) sum_(n=1)^t I(z_(n-1); z_n) ) . $

  于是给定 $x_0$ 后 $t$ 步的累积均方误差有下界

  $ EE_p[ norm(x_(1:t) - EE_(q_"KR")[x_(1:t)|x_0])^2 | x_0 ] >= C exp( - ( 2/(n t) sum_(n=1)^t I(z_(n-1); z_n) + cal(E) ) ) , $

  其中常数 $C = (n t)/(2e) exp( (2)/(n t) sum_(n=1)^t I(x_(n-1); x_n) )$ 吸收了轨迹的边缘熵。
]

#proof[
  经典率失真定理指出：$x_(1:t) in RR^(n times t)$ 在 $L^2$ 失真下的最小可达失真 $D$
  被 Shannon 下界控制。把轨迹视为 $n t$ 个分量、总码率由潜在信道提供，
  这个信道每步可传的信息量正是 $I(z_(n-1); z_n)$，于是

  $ D >= (n t)/(2 e) exp( (2)/(n t) H(x_(1:t)) ) exp( - (2)/(n t) sum_(n=1)^t I(z_(n-1); z_n) ) . $

  又因为 $H(x_(1:t))$ 完全可以由原始动力学 $T$ 的互信息和 $sum_(n=1)^t I(x_(n-1); x_n)$ 度量，
  代入并合并常数即得第二式。$square$
]

#insight("下界的意义：信息与精度的兑换率")[
  $E >= C exp(-(2)/(n t) sum_n I(z_(n-1); z_n))$，预测误差随保留的信息*指数*衰减。

  反过来读：*要让误差减少一个因子 $e^2$，潜在空间必须多保留 1 nat 的互信息*。
  这个“汇率”是后面所有实验解读的基础：任何预测精度的改进都可以翻译成
  “多保留了 $Delta R$ nat 关于未来的信息”。

  同时它也解释了为什么“用掉 $d$ 维但只用了 1 维”是真正的损失，
  多出来的维度本可以提供更多的 $I(z_(n-1);z_n)$，但被浪费了。
]

现在两个问题回答完了：信息一定会丢（命题 1），丢多少有上下界（命题 2）。
下一个问题是：*丢的是哪一类信息？* 这一问题把信息论接到 Koopman 的谱理论上。

#pagebreak()
