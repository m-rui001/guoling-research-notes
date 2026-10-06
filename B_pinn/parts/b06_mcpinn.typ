#import "../../preamble.typ": *

= 第四篇：分数阶算子的随机估计（MC-PINN，2022）

== 新困难：分数阶导数无法自动微分

回到开头的例三：带分数阶导数的方程。分数阶导数描述记忆效应与长程相互作用，
在多孔介质溶质输运、粘弹性力学里是标准模型。它有两个数值上的坏性质：

- *非局部*：$(partial^gamma u)/(partial t^gamma)$ 依赖 $u$ 的整段历史，
  每算一个点都要扫一遍过去；
- *核奇异*：积分核 $(t - tau)^(-gamma)$ 在 $tau -> t$ 处爆炸。

fPINN（Pang–Lu–Karniadakis 2019）的处理方式是：整数阶导数交给自动微分，
分数阶导数用 Grünwald–Letnikov 离散格式 + 求积点近似。
代价是每个训练点需要大量*辅助点*（评估 $u$ 在附近网格点的值）：
二阶 GL 格式在 $d$ 维需要 $q^d times L$ 次网络评估，
$q >= 8$ 是角向求积点数，$L approx 100$ 是辅助点数。三维已经勉强，十维直接不可能。

Guo, Wu, Yu, Zhou（2022）的替换方案：把分数阶算子*改写成随机量的期望*，
然后用蒙特卡洛采样估计。每个训练点的评估次数从 $q^d times L$ 降到 $8 m + 1$
（$m approx 20$ 个随机样本），与维数几乎无关。

== 两个分数阶算子

时间方向的 Caputo 导数（$0 < gamma < 1$）：

$ (partial^gamma u(x, t))/(partial t^gamma) = 1/(Gamma(1 - gamma)) integral_0^t (t - tau)^(-gamma) (partial u(x, tau))/(partial tau) d tau . $

空间方向的分数阶 Laplacian（$0 < alpha < 2$），由一个超奇异积分定义：

$ (-Delta)^(alpha\/2) u(x) = C_(d, alpha) "P.V." integral_(RR^d) (u(x) - u(y))/(norm(x - y)_2^(d + alpha)) d y , quad C_(d, alpha) = (2^alpha Gamma((alpha + d)\/2))/(pi^(d\/2) |Gamma(-alpha\/2)|) . $

逐项读：被积函数是“点 $x$ 的值与周围点 $y$ 的值之差，按距离的 $d + alpha$ 次幂衰减加权”。
距离近的邻居权重大（核奇异），远处的贡献衰减慢（非局部）。
$"P.V."$ 是柯西主值：积分在 $y = x$ 处发散，取对称极限后收敛；
$u(x) - u(y)$ 在 $y$ 接近 $x$ 时与 $norm(x - y)$ 同阶，恰好抵消奇性。

== 分数阶 Laplacian 的 MC 逼近

核心动作：把积分改写成“对随机方向 $xi$ 和随机半径 $r$ 的期望”。
以 $x$ 为球心、$r_0$ 为半径把积分分内外两段。

*内段*（$r in [0, r_0]$，被积函数奇性所在）：利用球坐标与对称化
（$y = x + r xi$ 与 $y = x - r xi$ 两项相加消掉奇次项），

$ integral_(y in B_(r_0)(x)) (u_(n n)(x) - u_(n n)(y))/(norm(x - y)_2^(d + alpha)) d y = (|S^(d-1)| r_0^(2 - alpha))/(2(2 - alpha)) E_(xi, r tilde f_I (r)) [ (2 u_(n n)(x) - u_(n n)(x - r xi) - u_(n n)(x + r xi))/(r^2) ] , $

其中 $xi$ 均匀分布于单位球面 $S^(d-1)$，半径的采样密度是

$ f_I (r) = (2 - alpha)/(r_0^(2 - alpha)) r^(1 - alpha) bold(1)_(r in [0, r_0]) , quad r\/r_0 tilde "Beta"(2 - alpha, 1) . $

*外段*（$r in [r_0, infinity)$，远处的贡献）：

$ integral_(y in.not B_(r_0)(x)) (u_(n n)(x) - u_(n n)(y))/(norm(x - y)_2^(d + alpha)) d y = (|S^(d-1)| r_0^(-alpha))/(2 alpha) E_(xi, r tilde f_O (r)) [ 2 u_(n n)(x) - u_(n n)(x - r xi) - u_(n n)(x + r xi) ] , $

$ f_O (r) = alpha r_0^alpha r^(-1 - alpha) bold(1)_(r in [r_0, infinity)) , quad r_0\/r tilde "Beta"(alpha, 1) . $

两个采样密度都是幂律，用 Beta 分布采样（$r \/ r_0 tilde op("Beta")(2 - alpha, 1)$
即半径 $r = r_0 beta^(1\/(2 - alpha))$，$beta tilde "Beta"$），一次采样一次评估。
权重常数恰好吸收了被积函数的 $r$ 依赖：这是重要性采样的标准构造，
采样密度正比于积分核的径向部分，让被估的量变得有界。

Caputo 导数同样处理：分部变换后

$ (partial^gamma u_(n n))/(partial t^gamma) = gamma/(1 - gamma) t^(1 - gamma) E_(tau tilde f_(I, t)) [ (u_(n n)(x, t) - u_(n n)(x, t - tau_epsilon t))/(tau_epsilon t) ] + (u_(n n)(x, t) - u_(n n)(x, 0))/(t^gamma) , $

$ f_(I, t)(tau) = (1 - gamma) tau^(-gamma) bold(1)_(tau in [0, 1]) , quad tau tilde "Beta"(1 - gamma, 1) . $

#warn("两个工程细节")[
  *舍入保护*：内段公式里 $r -> 0$ 时差商 $(2 u(x) - u(x - r xi) - u(x + r xi))\/r^2$
  趋于径向二阶导，数值上两个几乎相等的数相减会吃掉有效数字。
  实际用 $r_epsilon = max{epsilon, r}$ 设置下限。
  *硬边界*：边界条件可以编码进网络结构
  （$u_(n n) = "ReLU"(1 - norm(x)^2) tilde(u)_(n n)$ 自动满足球域上的零边界），
  省掉边界损失。
]

== 无偏损失（本篇的理论核心）

物理损失是残差的平方：$L_"equ" (theta) = norm(L[u_(n n)] - f)^2$。
但 $L[u_(n n)]$ 现在是随机估计，平方的期望不等于期望的平方：
直接对同一个随机估计平方取平均，得到的是 $E[hat(L)^2] = L^2 + op("Var")(hat(L))$，
多出的一块方差项会让梯度有系统性偏移。

解法是*独立双采样*：对同一个训练点做两次独立估计 $hat(L), hat(L)'$（各自独立采
$tau, xi, r$），损失用它们的乘积：

$ hat(L)_"equ" (theta) = 1/(m N_u) sum_(i, j) hat(L)[u_(n n)(x_i, t_i; theta)] dot hat(L)'[u_(n n)(x_i, t_i; theta)] . $

#proposition("无偏性")[
  设 $hat(L)$ 与 $hat(L)'$ 独立同分布，且都是 $L[u_(n n)] - f$ 的无偏估计
  （$epsilon = epsilon_t = 0$、舍入误差可忽略时成立），则

  $ E[hat(L)_"equ" (theta)] = E[hat(L)] E[hat(L)'] = (L[u_(n n)] - f)^2 = L_"equ" (theta) . $
]

#proof[
  两步。*第一步（独立性拆开期望）*：把两次估计写成
  $hat(L) = mu + delta$、$hat(L)' = mu + delta'$，其中
  $mu = L[u_(n n)] - f$ 是公共期望，$delta, delta'$ 是零均值的波动。
  乘积展开后逐项取期望：

  $ E[hat(L) hat(L)'] = mu^2 + mu E[delta] + mu E[delta'] + E[delta delta'] = mu^2 , $

  前两项由波动零均值消失；最后一项 $E[delta delta'] = E[delta] E[delta'] = 0$
  用独立性（独立随机变量期望的乘法性质，概率论标准结论）。
  方差项就在这里被排除：如果损失直接用 $hat(L)^2$，
  对应的展开会多出 $E[delta^2] = op("Var")(hat(L))$，梯度被系统性抬高。
  *第二步（无偏性）*：由第 10 章各估计公式的构造，
  每个 $hat(L)$ 的期望恰是它要估计的积分，所以
  $E[hat(L)] = E[hat(L)'] = mu$，相乘即得结论。$square$
]

证明本身短，构造才是重点：
*它把“随机化的算子”变成了“可以放心放进损失函数的对象”*。
同样的技巧在 NFF 的变分损失里已经出现过一次（第 8 章），
两篇论文独立使用它，说明这是随机化物理损失的通用模式。

#figure(
  align(center, image("../fig_b/tmp_mcsketch.png", width: 92%)),
  caption: [MC-PINN 总框架（论文 2203.08501）。对每个训练点：随机采样方向 $xi$
  与半径 $r$（内段与外段各一组分布），在 $x plus.minus r xi$ 处评估网络，
  构造分数阶算子的无偏估计；两次独立估计相乘构成物理损失。
  与 fPINN 的 GL 求积相比，评估次数从 $q^d times L$ 降到 $8 m + 1$，与维数几乎无关。],
)

== 实验

*正向*：分数阶 Laplacian 方程 $(-Delta)^(alpha\/2) u = f$ 在 2 维、3 维、10 维求解，
精确解 $u = (1 - norm(x)_2^2)^(1 + alpha\/2)$（Dyda 2012 的经典构造）。
误差在 $10^(-3)$ 附近饱和，且对 MC 样本数 $m$ 与截断半径 $r_0$ 都不敏感，
10 维也照常工作（GL 路线在此维数下需要约 $8^10 times 100$ 次评估，完全不可行）。

*反向*：时空分数阶平流扩散方程，同时识别四个参数
（阶数 $alpha, gamma$、速度 $c, v$）。一维观测下识别结果与真值
$(1.5, 0.5, 0.1, 1.0)$ 的偏差在千分位，维数升到 5 误差仍在 $3 times 10^(-3)$ 量级。

*不确定性*：随机参数分数阶扩散方程（$alpha, mu$ 均匀随机），
MC-PINN 充当廉价代理，再用 #term("ABC", "近似贝叶斯计算，按“模拟数据与观测的距离不超过阈值”筛选先验样本，避开写似然函数") 从 $10^5$ 个先验样本中
按传感器失配筛选，得到参数的后验分布，真值落在后验峰附近。

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_b/tmp_mcsample.png"), image("../fig_b/tmp_ade1d.png"))),
  caption: [左：相对 $L^2$ 误差随 MC 样本数 $m$ 的变化（2D/3D/10D，论文 2203.08501）。
  误差在约 $10^(-3)$ 处饱和，对 $m$ 不敏感，10 维与 2 维行为一致。
  右：反问题中四个参数 $(alpha, gamma, c, v)$ 的识别随训练的收敛，
  真值以虚线标出。],
)

== 遗留问题

MC 估计把算子随机化，损失面因此带噪声，训练只能到 $10^(-3)$ 量级；
分数阶参数 $alpha$ 出现在采样密度 $f_I$ 里，反向问题中它是可训练参数而采样分布
依赖它，需要小心处理梯度。另外，随机 + 时间依赖 + 分数阶三个难点
到目前还是分开解决的：模态空间管随机与时间（第 6 章），MC 管分数阶（本章）。
把它们拼起来是下一篇的事。
