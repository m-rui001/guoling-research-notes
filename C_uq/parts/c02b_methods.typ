#import "../../preamble.typ": *

= 经典 UQ 方法逐个精讲

综述把方法分成两族各五六种，主章节只给了一句话版。
这一章逐个展开：每种方法的机制、代价、失效模式。
它们是后续所有方法的对照基线，也是读综述比较实验的解码表。

== HMC：金标准及其代价

HMC（混合蒙特卡洛）把后验采样想象成一个物理系统：
权重 $theta$ 是粒子位置，引入动量 $p$，在“负对数后验”构成的
势能面上滚动。每步先按哈密顿动力学用 leapfrog 积分走
$L$ 步（可逆、保体积），再按 Metropolis 准则
（按新旧位置的密度比决定接受概率）接受或拒绝。

它解决随机游走采样的大问题：哈密顿轨迹一次可以横穿后验的
大范围，混合效率在高维远超逐坐标采样。
代价同样清楚：每个 leapfrog 步要算一次完整梯度
（对百万级参数的网络，一次前向 + 反向），
$L$ 步 × 5000 样本 × burn-in（预热：丢弃链初期还没进入平稳段的样本），
综述算例里 131.5 秒
（函数逼近这种小网络），PDE 场景要乘几个量级。
另一个实践坑：步长与 $L$ 要调（NUTS 是自动调步长与 $L$ 的 HMC 变体），
调不好要么拒绝率爆炸要么轨迹折返。

== 朗之万动力学：一步版

朗之万（LD）可以读成“带各向同性噪声的梯度下降”：
每步 $theta_(t+1) = theta_t + eta nabla log p(theta_t | cal(D)) + sqrt(2 eta) bold(z)_t$。
步长 $eta$ 趋零时马尔可夫链的不变分布（链跑久之后位置服从的
平稳分布，第 0 章 §0.6）恰好是后验。为什么是后验：
后验的对数梯度 $nabla log p$ 指向密度高处，噪声项按扩散方程
把概率从密度高处向低处搬运；对精确的连续时间过程，
这两股流在密度恰好等于 $p(theta | cal(D))$ 时互相抵消，
分布不再变化，所以后验是平稳解。
它比 HMC 便宜（每样本一次梯度），混合慢；
综述里 44.4 秒。常数步长版本不精确（不变分布有偏），
有偏性在 $eta$ 小时可以忽略但永不严格为零，
这也是它常被当成“便宜的近似”而不是“正确的方法”的原因。

“常数步长有偏”这句话可以用一维高斯目标算到底，
看清偏差从哪来、有多大。取目标 $cal(N)(0, 1)$：
此时 $log p(theta) = -theta^2 \/ 2$，梯度是 $-theta$，
一步更新为 $theta' = (1 - eta) theta + sqrt(2 eta) dot z$，
$z tilde cal(N)(0, 1)$。

#proof[
  设链已进入平稳状态，$theta tilde cal(N)(0, sigma^2)$，求 $sigma^2$。

  *第一步：算一步之后 $theta'$ 的分布。* $theta'$ 等于 $theta$ 的
  线性变换加一项独立高斯噪声。由第 0 章 §0.5 的线性变换封闭性，
  $(1 - eta) theta tilde cal(N)(0, (1 - eta)^2 sigma^2)$；
  两个独立高斯之和仍是高斯，方差相加（协方差为零，交叉项在期望中消失），
  所以 $theta' tilde cal(N)(0, (1 - eta)^2 sigma^2 + 2 eta)$。

  *第二步：解不变条件。* 平稳意味着 $theta'$ 与 $theta$ 同分布：
  $sigma^2 = (1 - eta)^2 sigma^2 + 2 eta$。把含 $sigma^2$ 的项移到左边，
  提出 $sigma^2$：
  $sigma^2 [1 - (1 - eta)^2] = 2 eta$。
  括号展开 $1 - (1 - 2 eta + eta^2) = eta (2 - eta)$，
  两边除以 $eta (2 - eta)$：

  $ sigma^2 = 2/(2 - eta) = 1/(1 - eta \/ 2) approx 1 + eta \/ 2 . $

  读法：链最终停留的分布方差比目标方差 1 大。
  $eta = 0.1$ 时 $sigma^2 approx 1.053$；步长减半，偏差也约减半。
]

#insight("用日常语言读一遍：偏差从哪来")[
  每一步都拿常数近似连续变化的漂移与噪声，这个线性化误差
  每步一点点积累，最后凝成不变分布里那份约 $eta \/ 2$ 的方差膨胀。
  高维情形结论同形：偏差是 $eta$ 的一阶小量，所以实践里靠调小步长压制，
  或者像 MALA 那样在每步后面加一次 Metropolis 接受拒绝
  （按新旧位置的密度比决定接受概率），把偏差精确修掉，代价是每步多算一次密度。
]


== MFVI：把采样变成优化

平均场变分推断用对角高斯 $q_omega (theta) = cal(N)(mu, "diag" sigma^2)$
近似后验，目标是最大化证据下界（ELBO；证据是边际似然
$p(cal(D))$ 的别名）。推导值得完整写出，
因为它是后面 VIB（IB-UQ）下界的骨架：

#proof[
  对任意分布 $q$，逐步展开 $log p(cal(D))$。

  *第一步：把常数写成积分。* $log p(cal(D))$ 不含 $theta$，是常数；
  任何分布满足 $integral q(theta) d theta = 1$。于是
  $log p(cal(D)) = log p(cal(D)) dot integral q(theta) d theta
    = integral q(theta) log p(cal(D)) d theta$：
  常数乘 1 不变，而常数可以自由移进积分号（对 $theta$ 的积分不改变它）。

  *第二步：插入贝叶斯恒等式。* 由乘法公式
  $p(cal(D), theta) = p(theta | cal(D)) p(cal(D))$（第 0 章 §0.3）
  得 $log p(cal(D)) = log p(cal(D), theta) - log p(theta | cal(D))$。
  再同时加减 $log q(theta)$：

  $ log p(cal(D)) = log(p(cal(D), theta) \/ q(theta)) + log(q(theta) \/ p(theta | cal(D))) . $

  读法：右端两项相加时 $log q(theta)$ 一负一正抵消，
  剩下的正是上式的 $log p(cal(D), theta) - log p(theta | cal(D))$。

  *第三步：两边对 $q$ 取期望。* 第一步的积分形式两边替换成第二步的右端：

  $ log p(cal(D)) &= integral q(theta) log(p(cal(D), theta) \/ q(theta)) d theta \
    &quad + integral q(theta) log(q(theta) \/ p(theta | cal(D))) d theta \
    &= "ELBO"(q) + op("KL")(q ‖ p(theta | cal(D))) . $

  第二个积分与第 0 章 §0.4 的 KL 定义逐字相同，它是 $q$ 与真后验的 KL；
  第一个积分按定义记作 ELBO。KL 非负（第 0 章 §0.4 的性质二），
  所以 $log p(cal(D)) >= "ELBO"(q)$。

  *第四步：从界到优化问题。* $log p(cal(D))$ 由数据决定，是固定常数
  （算不出它的值，但它不随 $q$ 变）。于是
  $max_q "ELBO"(q)$ 与 $min_q op("KL")(q ‖ p(theta | cal(D)))$
  是同一个问题：两者之和恒等于这个常数，ELBO 每升一分，KL 就降一分。

  *第五步：把 ELBO 拆成重构项与先验惩罚项。* 对第一个积分用
  $p(cal(D), theta) = p(cal(D) | theta) p(theta)$ 拆开对数，
  第三个积分按熵的定义（第 0 章 §0.4，$cal(H)(q) = -integral q log q$）改写：

  $ "ELBO"(q) &= integral q(theta) log p(cal(D) | theta) d theta
      + integral q(theta) log p(theta) d theta
      - integral q(theta) log q(theta) d theta \
    &= underbrace(E_(q(theta))[log p(cal(D) | theta)], "重构项：采样出的参数解释数据多好")
      - underbrace(op("KL")(q(theta) ‖ p(theta)), "变分分布偏离先验的惩罚") . $

  前一行到后一行的中间代数：后两个积分合并为
  $integral q(theta) log(p(theta) \/ q(theta)) d theta$，
  它按定义就是 $op("KL")(q ‖ p(theta))$。读法：
  优化 ELBO 就是在"拟合数据"与"不偏离先验"之间取加权。
  第 5 章的 VIB 目标与这一行同构。

  *第六步：为什么整个目标能反向传播。* 目标里有采样
  $theta tilde q_omega$，而"采样"这个操作对参数 $mu, sigma$ 没有导数。
  重参数化 $theta = mu + sigma dot epsilon$、$epsilon tilde cal(N)(0, 1)$
  把随机性转移到与参数无关的 $epsilon$ 上（两者同分布由第 0 章 §0.5 的
  线性变换封闭性保证）。于是目标写成 $E_(epsilon)[g(mu, sigma, epsilon)]$，
  是可导函数的期望；期望与求导交换需要的条件
  （被积函数对参数连续可微、且被可积函数控制）在高斯族与常见网络下成立，
  梯度就能用 $epsilon$ 的样本平均来估计。$square$
]

#insight("用日常语言读一遍：这六步合起来干了什么")[
  它把算不出的证据 $log p(cal(D))$ 换成了一个能算的下界：
  前四步证明 $log p(cal(D))$ 永远压着 ELBO，两者的差距恰好是
  $q$ 与真后验的距离；第五步把 ELBO 改写成"数据拟合减去先验惩罚"，
  这才是可以直接写进代码的形式；第六步解决"目标里有采样、怎么求梯度"。
  代价在平均场假设里：$q$ 的族太窄时，最优的 $q$ 只能贴住后验的一个峰。
]

失效模式来自平均场假设：对角高斯表达不了权重之间的后验相关，
不确定性往往被低估；深度网络的 ELBO 优化本身非凸，
$q$ 可能塌到后验的某一个峰（多峰后验的模式坍缩）。
综述里 19.5 秒，全表最便宜。

== MC dropout：不重新训练的不确定度

标准训练（带 dropout），推断时也开着 dropout 采样 $T$ 次，
统计预测均值与方差。Gal–Ghahramani 的等价性定理说这近似于
对网络权重做变分推断（Bernoulli 掩码的变分分布）。
它的吸引力是零额外训练成本；代价是“不确定度的质量”
依赖 dropout 率这个超参数，且 OOD 区域经常过置信
（综述实验的直接观察）。41.7 秒。

== 拉普拉斯近似：训练后的一次性计算

在训练终点 $hat(theta)$ 对负对数后验做二阶泰勒展开：
$-log p(theta|cal(D)) approx "常数" + 1/2 (theta - hat(theta))^T bold(H) (theta - hat(theta))$，
后验近似为 $cal(N)(hat(theta), bold(H)^(-1))$。
Hessian（损失对参数的二阶偏导数矩阵）太大
（百万量级参数的平方规模），实际用对角近似、
Kronecker 分解（K-FAC）或 Gauss-Newton 结构。
60.3 秒，精度常与 MFVI 相当。失效模式：单模态假设
（多峰后验只取一个峰）。

== 深度集成与其省钱变体

深度集成（DEns）训练 $M$ 个不同初始化的网络。
它的“贝叶斯解读”至今有争议：多初始化捕到的是损失面的
多个局部极小，与后验多峰有对应关系但没有等价定理。
实证上它是最强的便宜基线之一（359.8 秒，即 $M$ 个训练串行跑完）。
两个省钱变体：快照集成 SEns（31.9 秒，单轨迹取多个极小）
与 SWAG（36.6 秒，快照再拟合高斯采样）。

#warn("综述比较实验的三个可迁移结论")[
  一，*便宜方法 + 学好先验* 可以打平昂贵方法：
  用留出数据调噪声水平与先验强度后，MFVI/LA 的校准与精度
  接近 HMC。二，*深度集成的不确定度在 OOD 处经常不涨*：
  各成员在训练分布外的行为趋于一致（都外推得差不多差），
  散度不再标记误差。三，*单指标排名不可信*：
  dropout 在某些指标上大幅领先 HMC、在 OOD 校准上大幅落后，
  论文结论依赖指标选择。
]

== 非贝叶斯三替代（综述列出的未来方向）

*Conformal prediction*：不需要分布假设。拿一个留出集，
计算“真实值相对预测分位的偏差”的经验分位数 $q$，
新预测的区间宽度乘上 $q$。有限样本覆盖保证
（$P(u in "区间") >= 1 - alpha$）是分布自由的，
代价是区间只有宽度没有形状、且要求数据可交换
（把样本顺序打乱，联合分布不变；比独立同分布略弱的假设）。

*Evidential deep learning*：网络直接输出一个证据分布的参数
（回归用 Normal-Inverse-Gamma：一次前向同时给出
均值、aleatoric、epistemic 与总方差），推断最便宜；
训练需要精心设计的损失（对错配的证据重罚），校准难做。

*神经过程（NP）*：函数空间的潜变量模型：给一组
“输入-输出”对，推断潜变量 $z$ 的后验，再预测新点的条件分布。
UQ-SONet（第 9 章）的 cVAE 解码器与它同族，差别在后者
把 PDE 物理损失接了进来。

== 指标体系补全

主章节给了 RL2E / MPL / RMSCE 三个一级指标。补两个常被引用的：

*CRPS*（连续排序概率得分）：预测分布 $F$ 对真值 $u$ 的得分
$ "CRPS"(F, u) = integral (F(x) - bold(1){x >= u})^2 d x $，
它度量“预测 CDF 与以真值为中心的阶跃函数”的距离，
有解析式（对高斯预测），是气象学几十年的标准得分，
综合了校准与锐度。

*覆盖率与区间宽度（PIW）*：名义 95% 区间的实际覆盖率
与平均宽度。覆盖率达标只是底线（宽到荒谬也能达标），
覆盖率与宽度的联合才是质量，RMSCE 与 PIW
要成对报告的原因就在这里。

#figure(
  align(center, image("../fig_c/c1_metrics.png", width: 90%)),
  caption: [评价指标全景（综述 2201.07766）：按数据需求分三类。
  蓝色需要测试集（RL2E、MPL、RMSCE），橙色需要金标准参考解
  （NIP-G、KL-G），绿色两者都不需要（区间宽度、离散度）。
  实践建议是三类各取至少一个。],
)
