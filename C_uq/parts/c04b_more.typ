#import "../../preamble.typ": *

= 后继方法补充：推导、对照与画廊

== EDG 补充一：为什么 subVP 的转移核解析可采样

解码过程要求 $p_D (z_t | z_0)$ 可直接采样（否则每步要数值积分
Fokker–Planck 方程，simulation-free 就不成立）。subVP 型 SDE
$ d z_t = -1/2 beta(t) z_t d t + sqrt(beta(t)) d W_t $
是线性 SDE，条件转移是高斯：

$ z_t | z_0 tilde cal(N)(e^(-1/2 integral_0^t beta) z_0, (1 - e^(-integral_0^t beta)) bold(I)) . $

推导完整写出。记 $Lambda(t) = integral_0^t beta(s) d s$（累积系数，
下面反复出现）。线性 SDE 的解有显式公式：对两边从 $0$ 到 $t$ 积分，
套用一维线性 ODE 的常数变易法并补上噪声项，

$ z_t = e^(-Lambda(t)\/2) z_0 + integral_0^t e^(-(Lambda(t) - Lambda(s))\/2) sqrt(beta(s)) d W_s , $

读法：第一项是确定性部分的衰减，第二项把每个时刻的噪声
乘上"从它发生到 $t$ 的衰减因子"后累加。
这个公式直接给出后面要用的两条性质：随机积分的期望为零
（被积函数与 $W$ 适应且平方可积，随机积分的定义性质）；
随机积分的二阶矩公式 $E[(integral h d W)^2] = integral norm(h)^2 d s$。

*第一步：均值。* 对解取期望，第二项期望为零：

$ E[z_t] = e^(-Lambda(t)\/2) z_0 , $

读法：均值按指数衰减，速度由 $beta$ 的累积量决定。

*第二步：方差。* $op("Var")[z_t] = E[z_t^2] - (E[z_t])^2$。
$z_0$ 与随机积分独立（后者只含 $t > 0$ 的噪声增量），
平方展开后的交叉项期望为零；对第二项用二阶矩公式
（这里 $norm(h)^2 = e^(-(Lambda(t) - Lambda(s))) beta(s)$）：

$ op("Var")[z_t] = integral_0^t e^(-(Lambda(t) - Lambda(s))) beta(s) d s . $

*第三步：算出这个积分。* 链式法则给出
$d\/d s[-(Lambda(t) - Lambda(s))] = beta(s)$，所以被积函数可以收进一个导数：
$e^(-(Lambda(t) - Lambda(s))) beta(s) = d\/d s e^(-(Lambda(t) - Lambda(s)))$，
积分值等于原函数在两端的差：
$e^0 - e^(-Lambda(t)) = 1 - e^(-Lambda(t))$。于是

$ z_t | z_0 tilde cal(N)(e^(-1/2 integral_0^t beta) z_0, (1 - e^(-integral_0^t beta)) bold(I)) . $

“方差保持”名字的来历在公式里直接可见：
$t -> infinity$ 时 $e^(-integral_0^t beta) -> 0$，方差趋于 1，
与先验 $cal(N)(0, bold(I))$ 一致。
两个矩都解析，采样直接按公式来。$beta_min = 0.1$、
$beta_max = 20$ 是标准配置。

#insight("用日常语言读一遍：这三步合起来干了什么")[
  解码器每走一步都需要"从 $z_0$ 出发 $t$ 时刻后在哪里"的答案。
  这三步把答案从"解一个随机微分方程"降为"查一个高斯公式"：
  线性结构让解可以整条写出，期望和方差各自只剩一个能精确积分的式子。
  于是训练与采样全程不碰数值积分，这是 subVP 型 SDE 被选中的原因。
]

== EDG 补充二：概率流 ODE 与密度计算

重加权需要算 $p_E (z_0 | x)$：编码过程的密度。
score 网络除了定义逆时 SDE，还定义了一个*概率流 ODE*
（同一密度的确定性表示）：

$ d z_t = (f(z_t, t) - 1/2 g(t)^2 s(z_t, x, t; theta)) d t . $

#proposition("命题（概率流 ODE 的密度公式）")[
  沿 ODE 的流映射 $z_T |-> z_0$ 改变密度的规则是连续性方程：

  $ log p_E (z_0 | x) = log p_E (z_T | x) - integral_0^T nabla dot F_theta (z_t, t) d t , $

  $F_theta = f - 1/2 g^2 s$ 是 ODE 的速度场。散度用 Hutchinson
  估计器估计（随机向量 $epsilon$，
  $nabla dot F approx epsilon^T (partial F\/partial z) epsilon$）。
]

#proof[
  *第一步：连续性方程从质量守恒来。* 取任意区域 $Omega$，
  它里面的概率质量是 $integral_Omega p(z, t) d z$。质量随时间变化的唯一途径
  是穿过边界流动：流出边界 $partial Omega$ 的总流量是
  $integral_(partial Omega) p F dot bold(n) d A$（$bold(n)$ 是单位外法向，
  $F dot bold(n)$ 是垂直于边界的速度分量）。
  于是 $d\/d t integral_Omega p d z = -integral_(partial Omega) p F dot bold(n) d A$。
  右端由散度定理换成体积分 $-integral_Omega nabla dot (p F) d z$，移项：

  $ d/d t integral_Omega p d z + integral_Omega nabla dot (p F) d z = 0 . $

  这个等式对*每个*区域 $Omega$ 成立，被积函数不逐点为零就能造出反例区域，
  所以 $partial_t p + nabla dot (p F) = 0$。这就是连续性方程，
  它的全部内容是"概率质量既不凭空产生也不凭空消失"。

  *第二步：展开并除以 $p$。* 乘积法则给出
  $nabla dot (p F) = F dot nabla p + p nabla dot F$，方程变成
  $partial_t p + F dot nabla p + p nabla dot F = 0$。
  在密度为正处两边除以 $p$：
  $partial_t log p + F dot nabla log p = -nabla dot F$。

  *第三步：左端是沿轨迹的全导数。* ODE 的轨迹满足
  $d z_t \/ d t = F(z_t, t)$，多元链式法则给出

  $ d/d t log p(z_t, t) = partial_t log p + (d z_t \/ d t) dot nabla log p
    = partial_t log p + F dot nabla log p , $

  正是第二步方程的左端。所以沿轨迹
  $d log p \/ d t = -nabla dot F_theta (z_t, t)$：
  密度的对数变化率等于速度场散度的负值。

  *第四步：沿轨迹积分。* 从 $T$ 积到 $0$（时间倒行）：

  $ log p_E (z_0 | x) - log p_E (z_T | x) = integral_T^0 (-nabla dot F) d t
    = -integral_0^T nabla dot F_theta (z_t, t) d t . $

  端点 $p_E (z_T)$ 已知：条件 (b) 说它与 $p_D (z_T)$ 重合，
  而 $p_D (z_T)$ 由转移核解析可得。$square$
]

散度项没法逐维求偏导再求和（维数高时代价是每步 $d$ 次反向传播），
实际用 Hutchinson 估计器。它为什么无偏，两步算完。
记 $A = partial F \/ partial z$，$epsilon$ 的分量独立地以等概率取 $plus.minus 1$
（Rademacher 随机向量）：

$ E_epsilon[epsilon^T A epsilon] = E_epsilon["tr"(A dot epsilon epsilon^T)] = "tr"(A dot E_epsilon[epsilon epsilon^T]) = "tr"(A) = nabla dot F . $

第一步：标量等于自己的迹，再用迹的循环性质把 $epsilon$ 挪到一起；
第二步：期望与迹都是线性运算，可以交换；
$E[epsilon epsilon^T]$ 的第 $i, j$ 元是 $E[epsilon_i epsilon_j] = bold(1){i = j}$，
即单位矩阵。所以期望恰是 $"tr"(A) = nabla dot F$，估计无偏，
代价从 $d$ 次偏导数降为一次雅可比向量积。

这套“采样走 SDE、密度走 ODE”的双表示是扩散模型的标准结构，
EDG 把它用在潜空间而不是数据空间，于是数据空间的
高维性（LJ55 的 165 维）不直接冲击散度估计的成本。

#figure(
  align(center, image("../fig_c/x_c3_ising.png", width: 86%)),
  caption: [Ising 模型的连续松弛采样（论文 2401.02080）：固定潜变量
  $z_0$，从温度 $T = 2.0$ 生成到 $T = 2.7$ 的状态序列。
  相变附近（磁化翻转）被自然捕捉；离散构型按
  $pi(s | x) = product_i (1 + e^(-2 s_i x_i))^(-1)$ 从连续样本抽取。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_c/x_c3_rew.png"), image("../fig_c/x_c3_cost.png"))),
  caption: [左：重加权的效果（MoG2(i)）：能量直方图从未加权（偏差）
  修正到贴合目标，加权 MMD 从 0.50 降到 0.09。
  右：训练时间与 MMD 的对比：EDG 用更少时间达到更低 MMD。],
)

== FMRC 补充：误差链的完整叙述

主定理把 $cal(L)^"FMRC"$ 与算子误差连起来，中间的链条值得完整走一遍：

$ "速度场误差" epsilon^2 = integral_0^1 E[norm(v_i^* - v_(i, theta))^2] d s quad => quad W_2("联合分布") <= epsilon e^(integral L_s d s) quad => quad norm(cal(K)_tau - hat(cal(K))_L)_H <= sqrt(beta) epsilon e^(integral L_s) . $

链条有两个箭头。第一个箭头在笔记内拆成三步走完；
第二个箭头的两半各引一条已给出的引理（下一章砖一与砖四）。

*第一个箭头：速度场误差到分布误差。* 输入两条：训练收敛使
速度场误差的平方平均不超过 $epsilon^2$；速度场对状态 $L_s$-Lipschitz。

+ *耦合两条轨迹。* 让全速度场与瓶颈速度场从同一个起点出发；
  把 $X$ 分量的速度补一个零块，$X$ 在两条轨迹里走成同一条线，
  偏差全部集中在 $X_1$ 分量。记 $e(s) = norm(Y^s - hat(Y)^s)$，
  求导、用三角不等式、再对速度差用 Lipschitz 条件：

  $ d e \/ d s &<= norm(v(s, Y^s) - v_theta(s, hat(Y)^s)) \
    &= norm(v(s, Y^s) - v_theta(s, Y^s)) + norm(v_theta(s, Y^s) - v_theta(s, hat(Y)^s)) \
    &<= Delta_s + L_s e(s) , $

  第一行是三角不等式拆成"同一点的速度差"与"同一场里的状态差"两段；
  第二行第一段记作 $Delta_s$（两个速度场在同一点的差，
  平方沿时间的平均被 $epsilon^2$ 控制，这是训练误差假设），
  第二段用 $v_theta$ 的 Lipschitz 条件。
+ *解微分不等式（Gronwall）。* $e' <= L_s e + Delta_s$ 两边乘积分因子
  $e^(-integral_0^s L_r d r)$，左边变成一个函数的导数；从 $0$ 积到 $1$，
  再对起点取期望并沿 $Delta$ 用 Cauchy–Schwarz：

  $ (E[norm((X, X_1) - (X, hat(X)_1))^2])^(1/2) <= epsilon dot e^(integral_0^1 L_s d s) . $

  读法：指数因子是 Lipschitz 常数沿路径的累积放大，
  流走得越远、对状态越敏感，同样的速度场误差被放得越大。
+ *认出 $W_2$。* 左端是一条具体耦合的 $L^2$ 代价；
  $W_2$ 的定义是对一切耦合取下确界（下一章砖一给出定义式），
  下确界不超过任何一条具体的值，于是
  $W_2((X, X_1), (X, hat(X)_1)) <= epsilon e^(integral L_s d s)$。

*第二个箭头：$W_2$ 到算子范数。* 分两半。前一半是
$W_2$ 与 $dot(H)^(-1)$ 范数的比较（Benamou–Brenier 线性化加 Loeper 型引理，
即下一章砖一：密度有 $beta$ 倍体积支配时
$norm(mu - nu)_(dot(H)^(-1)) <= sqrt(beta) W_2(mu, nu)$；
$beta$ 在这里进场，多峰分布峰间的低密度沟谷正是它可能失效的位置）。
后一半是"弱形式算子范数被联合密度的 $dot(H)^(-1)$ 误差控制"，
依据是算子范数的上确界只遍历 $g(x) h(y)$ 形状的测试函数，
而 $dot(H)^(-1)$ 范数遍历全部测试函数；逐步展开见下一章砖四。
两个箭头接起来就是主定理的右端。

#insight("用日常语言读一遍：这条链保证了什么")[
  它把"网络训练得好不好"翻译成"降维有没有丢信息"：
  速度场误差先被流的敏感性指数放大成两个联合分布的搬运距离，
  再被密度支配条件折算成算子的弱范数误差。
  每一环的方向都朝保守一侧，所以损失收敛就能保证算子误差收敛；
  代价全在常数里（指数因子与 $sqrt(beta)$），见下面的两个折扣。
]

#warn("读这条链时打两个折扣")[
  一，$e^(integral L_s)$ 里的 Lipschitz 常数没有实用上界，
  理论意义是“形状正确”（误差确实随速度场误差消失）多于
  “数值精确”。二，第二段的密度上下界 $beta$ 在多峰分布处
  可能不成立（峰间密度低）。论文正稿只保留 $H$ 范数结果，
  谱误差分析在源码里被注释掉，原因是论证不完备，
  这本身是值得学习的诚实处理。
]

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_c/c4_swiss.png"), image("../fig_c/x_c4_meta.png"))),
  caption: [FMRC 的数据构造与聚类检查（论文 2408.15981）。
  左：Swiss roll 映射把二维七井势的轨迹卷进三维
  （观测空间是非线性流形）。右：PCCA+ 聚类结果与各簇的
  RC 值分布：簇内直方图窄、簇间间隙清晰。],
)

== UQ-SONet 与 LVM-GP 补充

#figure(
  align(center, image("../fig_c/x_c5_ens.png", width: 86%)),
  caption: [VIDON + 深度集成的对照实验（论文 2509.25646 附录）：
  均值预测准确，但不确定度带系统性偏窄：多次独立训练的
  分歧没有覆盖真实的条件波动。这是“集成低估条件不确定度”
  在算子学习里的又一次出现。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_c/x_c5_2dns.png"), image("../fig_c/x_c5_sde.png"))),
  caption: [左：二维 Navier–Stokes（涡量-速度形式，$nu = 0.001$）：
  UQ-SONet 的均值误差 1.89 到 2.69 对 VIDON 的 5.98 到 8.74
  （$times 10^(-2)$），置换不变嵌入在高维问题上的优势最明显。
  右：随机椭圆方程（算子本身随机）在传感器数不同时的
  条件预测：不确定度不随传感器加密消失：算子内在随机性
  被潜变量保留。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_c/tmp_c6_enc.png"), image("../fig_c/x_c6_lc.png"))),
  caption: [LVM-GP 的编码器结构与先验自适应（论文 2507.22493）。
  左：确定性特征 $overline(z)$、置信门 $m$ 与 GP 先验 $z_0$ 的插值结构
  （与 IB-UQ 同构，先验换成带空间相关的 GP）。
  右：GP 核长度尺度的训练曲线：从初始 1.0 自适应收敛到 0.58，
  先验自己学到了数据的相关长度。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_c/x_c6_inv.png"), image("../fig_c/x_c6_ens.png"))),
  caption: [左：一维非线性 Poisson 反问题（识别 $lambda = 0.7$）：
  深度集成的解出现明显伪振荡，LVM-GP 与 HMC 平滑且贴近真值。
  右：多孔介质流动问题三种方法的对照：LVM-GP 与 HMC 的
  不确定度带接近，深度集成在高噪声下失效。],
)
