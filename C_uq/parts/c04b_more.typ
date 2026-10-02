#import "../../preamble.typ": *

= 后继方法补充：推导、对照与画廊

== EDG 补充一：为什么 subVP 的转移核解析可采样

解码过程要求 $p_D (z_t | z_0)$ 可直接采样（否则每步要数值积分
Fokker–Planck 方程，simulation-free 就不成立）。subVP 型 SDE
$ d z_t = -1/2 beta(t) z_t d t + sqrt(beta(t)) d W_t $
是线性 SDE，条件转移是高斯：

$ z_t | z_0 tilde cal(N)(e^(-1/2 integral_0^t beta) z_0, (1 - e^(-integral_0^t beta)) bold(I)) . $

推导：线性 SDE 的均值满足 $d E[z_t]\/d t = -1/2 beta(t) E[z_t]$
（直接解出指数衰减），方差满足 $d op("Var") \/ d t = -beta op("Var") + beta$
（解出 $op("Var") = 1 - e^(-integral_0^t beta)$，即“方差保持”名字的来历：
$t -> infinity$ 时方差趋于 1，与先验一致）。
两个矩都解析，采样直接按公式来。$beta_min = 0.1$、
$beta_max = 20$ 是标准配置。

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
  连续性方程 $partial_t p + nabla dot (p F) = 0$（概率质量守恒）
  沿流线积分：密度的对数变化率等于速度场的负散度
  $d log p \/ d t = -nabla dot F$。从 $T$ 积到 $0$ 即得
  （端点 $p_E (z_T)$ 已知：条件 (b) 说它与 $p_D (z_T)$ 重合，
  而 $p_D (z_T)$ 由转移核解析可得）。$square$
]

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

第一段：flow matching 的标准稳定性（速度场误差经 Lipschitz
流放大成分布误差，Gronwall 型论证）。
第二段：$W_2$ 与 $dot(H)^(-1)$ 范数的比较（Benamou–Brenier 线性化：
$W_2$ 控制负一阶齐次 Sobolev 范数，密度有上下界 $beta$ 时等价）；
而弱形式算子范数恰好是 $dot(H)^(-1)$ 范数限制在
乘积测试函数 $g(x) h(y)$ 上的部分，所以被 $W_2$ 控制。

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
