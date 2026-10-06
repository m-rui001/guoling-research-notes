#import "../../preamble.typ": *

= 概念准备：不确定性从哪里来

== 为什么神经网络需要 UQ

神经网络正在进入计算科学的内核：替代数值求解器解 PDE（物理信息网络）、
学算子（#term("DeepONet", "神经算子：两个网络的内积，branch 网络吃输入函数的离散表示，trunk 网络吃查询点坐标")）、参数化随机过程。随之而来的问题是：
这些网络的预测误差怎么量化？传统数值分析有成熟的误差理论
（离散误差、收敛阶），神经网络没有对应物，它的误差来源更多：

- *数据*：观测有噪声、有缺失；
- *物理模型*：方程本身是对现实的近似（模型误设）；
- *网络结构*：超参数选择、过参数化；
- *优化*：非凸目标只找到局部极小，采样有随机性。

SciML 的 UQ 框架（Psaros 等，SIAM Review 2023）把不确定性分成两大类，
这个区分贯穿全部后续方法：

- *aleatoric（偶然不确定性）*：数据噪声带来的，*不可约*，数据再多也在；
- *epistemic（认知不确定性）*：有限数据与模型选择带来的，*可约*，
  数据越多越小。

两者的区分决定下游动作：传感器读数的抖动（aleatoric）
不该激励你去补实验；网络在外推区的大幅波动（epistemic）应该。

#figure(
  align(center, image("../fig_c/c1_pie.png", width: 72%)),
  caption: [总不确定性的来源分解（综述 2201.07766 的饼图）：数据（噪声、缺失）、
  物理模型（误设、随机性）、网络（结构、超参数、过参数化）、后验推断误差。
  其中模型误设一项是该框架没有显式建模的，也是后续多篇工作
  （如模型修正方向）要补的缺口。],
)

== 统一的问题设定

综述给出一个覆盖四类问题的统一格式：

$ cal(F)_(lambda, xi)[u(x; xi)] = f(x; xi), quad x in Omega, quad xi in Xi, quad cal(B)_(lambda, xi)[u(x; xi)] = b(x; xi), quad x in Gamma . $

四类场景由“谁未知”区分：

1. *正向确定性*：算子与参数已知，求 $u$；
2. *混合确定性*：$lambda(x)$ 与 $u(x)$ 都部分未知，共同重构；
3. *混合随机*：算子随机，求 $u$ 与 $lambda$ 的统计量；
4. *算子学习*：$cal(F), cal(B)$ 本身未知，先预训练再处理有限含噪推断数据。

== 贝叶斯模型平均：参数不确定性怎么变成预测不确定性

设网络权重为 $theta$，数据为 $cal(D)$。贝叶斯框架下权重有后验分布
$p(theta | cal(D)) prop p(cal(D) | theta) p(theta)$：
先验 $p(theta)$ 是看数据之前对权重的假设，似然 $p(cal(D) | theta)$
是给定权重时数据出现的机会，后验是两者按贝叶斯定理合成、
看数据之后的更新（硬币版完整数值例见第 0 章 §0.3）。
预测的分布由 #term("BMA", "贝叶斯模型平均") 给出：

$ p(u | x, cal(D)) = integral p(u | x, theta) p(theta | cal(D)) d theta , quad p(cal(D)) = integral p(cal(D) | theta) p(theta) d theta . $

这个积分算不出来（$theta$ 维数是百万级），实践用
#term("MC", "Monte Carlo，蒙特卡洛：从分布抽样本、用样本平均代替积分，第 0 章 §0.6") 近似：
从后验抽 $M$ 个样本 $hat(theta)_1, ..., hat(theta)_M$，预测分布近似为混合

$ p(u | x, cal(D)) approx overline(p)(u | x) = 1/M sum_(j=1)^M p(u | x, hat(theta)_j) . $

预测均值与总方差有显式公式（方差的分解定律）：

$ hat(u)(x) approx overline(mu)(x) = 1/M sum_(j=1)^M u_(hat(theta)_j)(x) , $

$ op("Var")(u | x, cal(D)) approx overline(sigma)^2(x) = underbrace(Sigma_u^2, overline(sigma)_a^2(x)) + underbrace(1/M sum_(j=1)^M (u_(hat(theta)_j)(x) - overline(mu)(x))^2, overline(sigma)_e^2(x)) . $

逐项读：总方差 = aleatoric 项（噪声水平 $Sigma_u^2$，各成员共享）加
epistemic 项（各成员预测在均值附近的散度）。*两类不确定性从此有了
可计算的化身*：网络组内分歧就是 epistemic 的度量。异方差
（heteroscedastic：噪声水平随输入 $x$ 变化）版本把
$Sigma_(u, theta)^2(x)$ 也作为网络的额外输出（每点一个噪声水平），
$overline(sigma)_a^2$ 改为对各成员噪声输出的平均。

#figure(
  align(center, image("../fig_c/c1_total.png", width: 88%)),
  caption: [从 BMA 到预测分布的四步（综述 2201.07766）：参数后验 → MC 采样
  → 各成员预测叠加噪声 → 混合出总预测分布。图中两个预测成员的散度
  给出 epistemic 项，噪声宽度给出 aleatoric 项。],
)

各种 UQ 方法的差别，只在于“怎么得到后验样本”这一步：
HMC 与朗之万做 #term("MCMC", "马尔可夫链蒙特卡洛：构造平稳分布等于后验的随机链，链上位置当样本，第 0 章 §0.6")；
变分推断（MFVI、MC dropout）用可优化分布近似；
深度集成直接训练多个网络取它们的极小值点；拉普拉斯近似在训练终点
局部展开。综述把这些统一进一个“模型 + 物理 + 后验推断”的框架
（第 2 章的图），下一节展开。

== 校准：被忽略的第二个维度

有了不确定度还不够，它必须*诚实*。#term("校准", "calibration") 的定义：

#definition("什么叫“校准”")[
  理想的概率校准满足：真实值落在预测分布 $p$ 分位数以内的频率等于 $p$。
  实用判据：对每个分位水平 $p in [0, 1]$，

  $ P(U_X <= overline(P)_X^(-1)(p)) = p , $

  校准误差在分位水平上取均方根（RMSCE）：

  $ "RMSCE" = sqrt(E_(p in [0,1]) (p - hat(p)(p))^2) approx sqrt(1/N_p sum_(j=1)^N_p [p_j - 1/N_"test" sum_(i=1)^N_"test" bold(1){u_i <= hat(u)(x_i)_(p_j)}]^2) . $
]

#figure(
  align(center, image("../fig_c/c1_misscal.png", width: 86%)),
  caption: [过置信与欠置信（综述 2201.07766）。左：过置信模型的预测区间
  比数据生成分布窄，校准曲线在低位就达到 1；右：欠置信相反。
  中间的对角线是理想校准。深度网络的常态是过置信，尤其在
  OOD（out-of-distribution，分布外：输入落在训练数据没有覆盖的区域）区域。],
)

事后校准有三条标准路线，都只需要一小块留出数据：
方差缩放（Levi 等：找一个最优常数 $s$，用 $cal(N)(overline(mu)(x), s^2 overline(sigma)^2(x))$
替换预测分布，只改方差）；
CDF 校准（Kuleshov 等：拟合失准函数 $Q(p) = P(U_X <= overline(P)_X^(-1)(p))$ 并把它
作用到所有预测 CDF 上，均值方差都改）；
CRUDE（Zelikman 等：对标准化残差 $(u - overline(mu))\/overline(sigma)$ 拟合经验分布并回代）。
这三条在第 2 章的比较实验里都会出场。
