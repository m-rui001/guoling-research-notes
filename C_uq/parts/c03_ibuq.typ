#import "../../preamble.typ": *

= 方法核心：信息瓶颈做 UQ（IB-UQ，2023）

== 它要解决的问题

综述的框架给出预测分布的构造方式（后验样本 + 混合），但构造方式本身
没回答一个关键问题：*网络在 OOD 输入上的不确定度可信吗？*
贝叶斯方法与集成方法的参数分歧在训练数据附近有意义，
在外推区未必变大（综述自己的实验显示深度集成的 OOD 不确定度经常不涨）。

Guo–Wu–Zhou–Wang–Zhou（CMAME 2023）从信息瓶颈（IB）出发构造 UQ，
出发点是一个观察：IB 的目标天然包含“表示要简单”这一项，
而“简单表示对没见过的输入无话可说”，这恰好是 OOD 不确定度想要的形状。

先回顾 IB 的目标。经典 IB（Tishby 1999）寻找一个表示 $Z$，最大化

$ I(Z; Y) - beta I(Z; X) , $

即“保留关于目标 $Y$ 的信息，压缩关于输入 $X$ 的信息”，$beta$ 是权衡系数。
Alemi 等（2017）的 VIB 用变分下界把它变成可训练的损失。
IB 谱系与本文的关系逐条对照：

- Tishby–Pereira–Bialek（1999）：IB 原始定义；Tishby–Zaslavsky（2015）
  用 IB 分析深度网络的逐层表示；
- Alemi 等（2017）：VIB 变分近似，本文的下界推导直接沿用其技术；
- Alemi 等（2018）：指出 VIB 天然蕴含不确定性、可改善校准与 OOD 检测，
  这是本文的动机文献；
- 本文的改造有两处，见下文（压缩项换成宽分布版本；编码器换成置信度混合）。

== 改造一：压缩项换成宽分布版本

原始 IB 压缩 $I(Z; X)$：惩罚表示对输入的依赖。但“对输入的依赖”
在 OOD 处应该是*消失的*，惩罚总量做不到这一点。本文的数据增强版目标：

$ cal(L)_"IB" = I(Z; Y) - beta I(tilde(Z); tilde(X)) , $

其中 $tilde(X) tilde tilde(p)(x)$，而 $tilde(p)(x) prop p(x)^(1\/tau)$（$tau > 1$）
是比训练分布*更宽*的分布（加温：$tau$ 越大越平）。

#insight("为什么宽分布压缩就是 OOD 保险")[
  在 $tilde(p)$ 下做压缩，等于要求“在训练分布及其周边（含外推区）都不能
  依赖输入的具体取值”。最优解会在低密度区自动放弃编码，
  于是 OOD 输入的表示退化为无信息的先验样本，预测分布退化为先验，
  不确定度自动增大。本节末尾的命题把这个论证变成不等式。
]

“加温分布”如何实现？用 *GIN*（General Incompressible-flow Networks，
保体积的 normalizing flow）：它类似 RealNVP 但雅可比行列式恒为 1，
于是变换不改变密度值，只改变密度的*位置*。把可训练参数 $W$ 缩放为
$sqrt(tau) W$，就得到密度 $p_V (v)^(1\/tau)$ 的精确采样器。
“$1\/tau$ 次幂”从哪来，三步说清：

+ *密度换算没有雅可比因子。* 设 $v = g(x)$ 是 GIN 的变换。
  一般的 flow 换密度要乘 $|det partial g \/ partial x|$；
  GIN 的行列式恒为 1，所以 $p(x) = p_V (g(x))$：密度值沿着变换原样搬运
  （第 0 章 §0.9 的换算公式取雅可比因子为 1 的特例）。
+ *取幂等于除指数。* 若基分布是高斯，$p_V (u) prop exp(-norm(u)^2 \/ 2)$，
  则 $p_V (u)^(1\/tau) prop exp(-norm(u)^2 \/ (2 tau))$：
  密度取 $1\/tau$ 次幂，就是把指数上的能量除以 $tau$。
+ *除以 $tau$ 等于放大方差。* $exp(-norm(u)^2 \/ (2 tau))$
  作为 $u$ 的函数正比于 $cal(N)(0, tau bold(I))$ 的密度。
  所以采样程序是：抽 $v tilde cal(N)(0, tau bold(I))$
  （等价于抽标准高斯再乘 $sqrt(tau)$，第 0 章 §0.5 的线性变换封闭性），
  令 $x = g^(-1)(v)$。换回 $x$ 的密度时本应再乘雅可比因子，
  但它是 1，于是 $x$ 的密度恰为 $exp(-norm(g(x))^2 \/ (2 tau)) prop p(x)^(1\/tau)$。

第三步里"再乘雅可比因子"这一处正是普通 flow 做不到的地方：
行列式不是 1 时会多出一个无法吸收的因子，加温就不再精确。

== 改造二：置信度混合编码器

编码器显式地把“离训练数据多远”编码成一个逐维的门控：

$ z = "diag"(m(x)) overline(z)(x) + "diag"(bold(1) - m(x)) z_0, quad z_0 tilde cal(N)(0, I) , $

$m(x), overline(z)(x)$ 都是神经网络，$m(x) in [0, 1]^d$ 是逐维置信度。
这两个量为什么能当"门"用，看一下 $z | x$ 的分布：固定 $x$ 后
$m(x)$ 与 $overline(z)(x)$ 是确定的数与向量，$z$ 是先验噪声 $z_0$
的仿射变换，由第 0 章 §0.5 的线性变换封闭性，$z | x$ 仍是高斯：
均值是 $m(x) overline(z)(x)$（噪声项均值为 0），
协方差是 $"diag"(bold(1) - m(x)) dot "diag"(bold(1) - m(x)) = "diag"((1 - m(x))^2)$
（对角矩阵相乘就是对角元相乘）。于是
$op("Var")(z_i | x) = (1 - m_i (x))^2$：
门控的每个分量直接就是该维潜变量的方差参数。
两个极端：$m = bold(1)$（训练数据附近）时 $z$ 是确定性特征 $overline(z)(x)$；
$m = bold(0)$（OOD）时 $z$ 退化为纯噪声 $z_0$，编码器不再携带任何输入信息。
解码器是条件高斯 $q_D (y | z) = cal(N)(mu_D (z), Sigma_D (z))$，
均值与对角协方差都由网络输出。推断时的预测分布
$hat(p)(y | x) = integral q_E (z | x) q_D (y | z) d z$ 用 MC 采样计算。

#figure(
  align(center, image("../fig_c/tmp_c2_ab.png", width: 90%)),
  caption: [IB-UQ 回归框架总图（论文 2302.03271）。GIN 从训练分布生成
  加温样本 $tilde(x)$；编码器输出确定性特征 $overline(z)$ 与置信门 $m$，
  与先验噪声 $z_0$ 混合成 $z$；高斯解码器输出预测均值与方差；
  变分损失 $cal(L)_"VIB"$ 由似然项与 KL 型压缩项组成。],
)

== 变分下界（推导完整保留）

IB 目标里的两个互信息都算不出，各用一个变分量替换。记号：
$q_E (z | x)$ 是编码器，$q_D (y | z)$ 是解码器，$e(z)$ 是任取的边际密度代理
（实现里用一个 RealNVP flow）。

*第一项：$I(Z; Y)$ 的下界。* 记号先说清：联合分布由模型与数据分布给出，
$p(y, z) = integral p(y | x) q_E (z | x) p(x) d x$
（从数据分布抽 $x$、编码得 $z$、解码出 $y$，三步合成的联合），
推导只用到它的两个边缘 $p(y)$ 与 $p(z)$，以及条件分布 $p(y | z)$。
从互信息的定义 $I(Z; Y) = E_(p(y, z))[log(p(y, z) \/ (p(y) p(z)))]$ 出发
（第 0 章 §0.4 的定义式按连续情形积分改写）：

*第一步：拆成两个条件对数。* 分子分母同除 $p(z)$，
$p(y, z) \/ (p(y) p(z)) = p(y | z) \/ p(y)$，于是

$ I(Z; Y) = E_(p(y, z))[log p(y | z)] - E_(p(y, z))[log p(y)] . $

*第二步：把第二个期望里的 $z$ 积掉。* $log p(y)$ 不含 $z$，
对联合分布取期望可以先对 $z$ 积分：
$E_(p(y, z))[log p(y)] = E_(p(y))[log p(y)]$，
按熵的定义（第 0 章 §0.4）它等于 $-cal(H)(Y)$。

*第三步：把 $log p(y | z)$ 拆成解码器加比值。*
$log p(y | z) = log q_D (y | z) + log(p(y | z) \/ q_D (y | z))$
（右端展开时 $log q_D$ 一正一负抵消）。取期望时第二项在固定 $z$ 下
对 $y$ 的平均恰是 $op("KL")(p(y | z) ‖ q_D (y | z))$（KL 的定义），
再对 $z$ 取平均：

$ I(Z; Y) &= E_(p(y, z))[log q_D (y | z)] + E_(p(z))[op("KL")(p(y | z) ‖ q_D (y | z))] + cal(H)(Y) . $

*第四步：取下界。* 中间一项是 KL 的期望，KL 非负（第 0 章 §0.4 的性质二），
扔掉它不等号朝"低估相关性"的方向走：

$ I(Z; Y) - cal(H)(Y) >= hat(I(Z; Y)) = E_(p(y, z))[log q_D (y | z)] , $

等号当且仅当每个 $z$ 处 $q_D (y | z) = p(y | z)$
（KL 为零当且仅当两个分布相同）。*最大化解码似然就是在最大化
相关性的下界*。

*第二项：$I(tilde(Z); tilde(X))$ 的上界。*
这里的联合分布由编码器自己定义：
$tilde(p)(tilde(x), tilde(z)) = tilde(p)(tilde(x)) q_E (tilde(z) | tilde(x))$
（抽一个加温输入，编码器给它一个随机表示），
所以条件分布 $tilde(p)(tilde(z) | tilde(x))$ 就是 $q_E$ 本身。互信息为

$ I(tilde(Z); tilde(X)) = E_(tilde(p)(tilde(x), tilde(z)))[log q_E (tilde(z) | tilde(x))] - E_(tilde(p)(tilde(z)))[log tilde(p)(tilde(z))] , $

第一项与上一项的第一步同形；第二项里的边缘 $tilde(p)(tilde(z))$
是编码器输出在加温分布上的真实边缘，要对 $tilde(x)$ 积分，算不出。

*第一步：插入代理边缘。* 把 $log tilde(p)(tilde(z))$ 拆成
$log e(tilde(z)) + log(tilde(p)(tilde(z)) \/ e(tilde(z)))$
（展开时 $log e$ 一正一负抵消），代入：

$ I(tilde(Z); tilde(X)) &= E[log q_E] - E_(tilde(p)(tilde(z)))[log e(tilde(z))] \
  &quad - E_(tilde(p)(tilde(z)))[log(tilde(p)(tilde(z)) \/ e(tilde(z)))] . $

*第二步：认出 KL 并取上界。* 第三个期望按定义是
$op("KL")(tilde(p)(tilde(z)) ‖ e(tilde(z)))$，KL 非负，
扔掉它不等号朝"高估压缩项"的方向走：

$ I(tilde(Z); tilde(X)) <= hat(I(tilde(Z); tilde(X))) = E_(tilde(p))[log q_E (tilde(z) | tilde(x)) - log e(tilde(z))] , $

等号当且仅当 $e(tilde(z)) = tilde(p)(tilde(z))$（KL 为零）。*最小化“编码器对宽分布样本的
编码长度减去代理边际的编码长度”就是在最小化压缩项的上界*。

*合并。* 两个界都朝安全方向偏（相关性项低估、压缩项高估），得到可训练的
VIB 目标：

$ cal(L)_"VIB" = E_(x, y, z)[log q_D (y | z)] - beta E_(tilde(x), tilde(z))[log q_E (tilde(z) | tilde(x))\/e(tilde(z))] <= cal(L)_"IB" - cal(H)(Y) . $

为什么 $cal(H)(Y)$ 不再出现：它只由数据的边缘分布决定，
不含 $q_D, q_E, e$ 的任何参数，对优化而言是常数，
从目标里去掉不改变极小值点。

训练算法（每步）：抽 mini-batch；可选 Mixup 平滑（$lambda tilde "Beta"(alpha, alpha)$，
$alpha = 0.005$ 取得极小，只起平滑经验分布的作用）；
GIN 抽 $tilde(x)$；编码得 $z, tilde(z)$；用样本平均估计上式；
梯度上升更新 $q_D, q_E, e$ 的全部参数。

#insight("用日常语言读一遍：这两条界合起来干了什么")[
  互信息的两个项都要求知道真实的条件分布与真实的边缘分布，都算不出。
  推导把每一项拆成"一个能算的期望"加"一个非负的 KL"：
  第一项扔掉 KL 得到相关性的下界，第二项扔掉 KL 得到压缩的上界，
  两个"扔掉"的方向都让目标变得保守，所以最大化这个替代目标
  不会把模型引向高估自己的方向。训练时唯一要算的是
  三个网络的输出对数，样本平均即可估计，这就是整个框架可以端到端训练的原因。
]

== 数据增强的理论保证

#proposition("命题（OOD 处的编码退化）")[
  记 $cal(M) = (q_E, q_D, e)$，最优模型 $cal(M)^*$。若存在常数 $C > 1$ 使
  $p(y)^(-1) p(y | x) <= C$ 对所有 $x, y$，则对满足 $p(x) <= epsilon tilde(p)(x)$
  的 $x$ 与 $epsilon$：

  $ op("KL")(q_E^*(dot | x) ‖ e^*(dot)) <= epsilon beta^(-1) log C . $
]

读法：$op("KL")(q_E^* ‖ e^*)$ 衡量“编码器输出与代理边际的偏离”，
也就是编码器还在多大程度上“依赖这个输入”。不等式说：
在训练密度比加温密度低 $epsilon$ 倍的位置（外推区），
最优编码器的输出与输入几乎独立，$z$ 落回先验 $z_0$，
预测分布退化为先验，不确定度由解码器在纯噪声输入下的行为决定
（自动变大）。这个命题是“宽分布压缩”设计的数据理论保证。

#figure(
  align(center, image("../fig_c/tmp_c2_beta.png", width: 88%)),
  caption: [信息平面分析（论文 2302.03271）：不同 $beta$ 下相关性下界
  $hat(I(Z;Y))$ 与压缩上界 $hat(I(tilde(Z);tilde(X)))$ 的轨迹。
  $beta$ 小时压缩不足（点右上），$beta$ 大时过度压缩（点左下），
  实验里的最优区间是 $beta in [0.2, 0.3]$。],
)

== DeepONet 版（算子学习）

算子学习的 UQ 把瓶颈放在 branch 网络（处理传感器测量的一侧）。
编码器同构：$z(o) = "diag"(m(o)) overline(z)(o) + "diag"(bold(1) - m(o)) z_0$，
$o$ 是传感器读数向量。解码器用确定性 DeepONet 同时输出均值与对数标准差：

$ cal(F)_theta (u)(y) = sum_(i=1)^n b_i (o) t_i (y) , quad (mu_D (y, z), log sigma_D (y, z))^T = sum_(i=1)^n b_i (z) t_i (y) , $

IB 目标变为 $cal(L)_"IBONet" = I((Y, Z); s(Y)) - beta I(tilde(Z); tilde(O))$，
变分下界与回归版逐字平行（把 $x$ 换成 $o$、$y$ 换成采样点处的 $s(y)$）。
两个 remark：只做单点 $y$ 的 UQ（多点联合预测是开放问题）；
trunk 网络也可以加瓶颈，本文选择只在 branch。

== 实验检验

四个场景，每个都把 IB-UQ 与 GP、HMC、深度集成对照：

- *不连续函数*（32 个含噪样本）：IB-UQ 与 GP/HMC 一样在外推区
  不确定度显著增大；深度集成的 OOD 不确定度不明显增大，
  印证了综述的“集成盲区”；
- *算子学习*（扩散反应方程与对流方程，OOD 输入取不同的相关长度）：
  误差被 2 倍标准差覆盖，输入越粗糙误差越大且不确定度越大；
- *California housing*：用 #term("LOF", "Local Outlier Factor，局部离群因子：按每个点的邻居密度给它打异常分") 划分 ID/OOD，
  误差随 LOF 降低（越异常）增大，不确定度同步覆盖；
- *气候模型*（温度场到气压场，删掉部分季度数据做 OOD）：
  IB-UQ 的不确定度在缺失区域增大，而普通 DeepONet 不输出任何置信度。

#warn("IB-UQ 的边界")[
  三点。第一，GIN 加温与 Mixup 都引入新超参数（$tau, alpha$），
  虽然实验显示结果对它们不敏感，但理论没有给出选择准则。
  第二，$e(z)$ 的表达力决定压缩上界的紧度，表达不足时上界松。
  第三，OOD 命题的常数条件（$p(y)^(-1) p(y|x) <= C$）在重尾目标下
  可能不成立。这些都是后续工作（第 9、10 章）继续改造的入口。
]
