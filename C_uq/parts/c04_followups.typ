#import "../../preamble.typ": *

= 后继一：Boltzmann 分布的扩散采样器（EDG，2024）

== 问题：没有样本的生成建模

统计力学里的标准对象是 Boltzmann 分布

$ pi(x) = 1/Z exp(-U(x)) , quad Z = integral e^(-U(x)) d x , $

配分函数 $Z$ 是高维积分，通常算不出。难点与普通生成模型不同：
*没有训练数据*，拿到 $pi$ 的样本正是任务本身。于是 JSD、MMD、Wasserstein 这些“真样本对假样本”的分布差异度量（分别是 Jensen–Shannon 散度、最大均值差异、最优传输距离）全部不可用，只剩 KL 可选；
MCMC（HMC、朗之万）渐近无偏但混合慢；变分路线（Boltzmann Generator）
用 normalizing flow 最小化 $D_"KL" (p_D ‖ pi)$，但 flow 的双射性限制表达力；
扩散采样（PIS 等）训练中要反复数值积分 SDE，昂贵。

Wang–Guo–Wu–Zhou（2024）的 EDG（energy based diffusion generator）：
*VAE 框架装一个潜空间扩散编码器*。解码过程不需要迭代采样，
编码过程（score 网络）训练时不需要数值积分 SDE（simulation-free），
两头都省。

== 框架与主定理

解码过程（生成方向）：$z_0 tilde p_D (z_0) = cal(N)(0, I)$，
$x | z_0 tilde p_D (x | z_0; phi)$（可训练解码器），
然后潜变量沿一个 SDE 扩散：

$ d z_t = f(z_t, t) d t + g(t) d W_t, quad t in [0, T] . $

两个设计条件：(a) $p_D (z_t | z_0)$ 可解析采样（选 subVP 型 SDE，即分数方差保持型，score-based 模型的标准选型，
好处是转移核解析可采样）；(b) $T$ 足够大使得 $z_T$ 近似与 $z_0$ 独立（信息耗尽）。

编码过程（推断方向）：条件 score 网络
$s(z_t, x, t; theta) approx nabla_(z_t) log p_D (z_t | x)$，
沿逆时 SDE 从 $(x, z_T)$ 走回 $(x, z_0)$。

#proposition("定理（VAE 型 KL 上界，主定理）")[
  $ D_"KL" (p_D (x) ‖ pi(x)) &<= D_"KL" (p_D (x, z_[(dot)]) ‖ p_E (x, z_[(dot)])) \
    &= cal(L)(theta, phi) + integral_0^T g(t)^2/2 E_(p_D)[norm(nabla_(z_t) log p_D (z_t))^2] d t + log Z , $

  $ cal(L)(theta, phi) = E_(p_D)[log p_D (x | z_0; phi) + U(x)] + integral_0^T g(t)^2/2 E_(p_D)[norm(s(z_t, x, t; theta))^2 + 2 nabla_(z_t) dot s(z_t, x, t; theta)] d t , $

  且等号成立当且仅当三条同时满足：解码过程中 $z_0$ 与 $z_T$ 独立，
  编码过程在 $T$ 时刻的边缘与解码过程重合（$p_E (z_T) = p_D (z_T)$），
  且 $s$ 恰为条件 score $nabla_(z_t) log p_D (z_t | x)$。
]

证明分两段：第一段是 VAE 型上界加独立性（附录证明 $x perp z_T$）；
第二段用 #term("Girsanov 定理", "随机分析中处理 SDE 漂移项变化的测度变换定理，给出两个漂移对应的分布之间的密度比") 处理逆时 SDE 的测度变换，
配合分部积分恒等式 $E_p[s^T nabla_z log p] = -E_p["tr"](partial_z s)$
把 score 的期望项化成散度项。最后两块与参数无关，
所以训练损失就是 $cal(L)(theta, phi)$。

训练用 #term("Hutchinson 估计器", "用随机向量的期望估计矩阵的迹或算子的散度，代价远低于显式求和") 估计散度项（Rademacher 随机向量 $epsilon$）：

$ cal(L)_t = (T g(t)^2)/2 (norm(s)^2 + 2 (partial[epsilon^T s])/(partial z_t) epsilon) , $

时间 $t$ 均匀采样，整个训练过程只从 $p_D (z_t | z_0)$ 精确采样，
不需要数值积分任何 SDE。*与 PIS 的差别就在这一句*。

== 样本重加权与配分函数

训练分布 $p_D$ 与目标 $pi$ 有偏离，用重要性采样去偏：以 $p_D (x, z_0)$
为提案，未归一化权重

$ w(x, z_0) = (exp(-U(x)) p_E (z_0 | x))/(p_D (z_0) p_D (x | z_0)) , $

$p_E (z_0 | x)$ 用概率流 ODE 与散度估计计算。
副产品是一个配分函数下界 $log Z >= E_(p_D)[log w]$，
它同时充当训练质量的指示器（越大越好）。

网络设计有两个针对物理问题的细节：score 模型把 $t = 0, T$ 两端的
已知边界 score 显式嵌入参数化（两端误差为零）；
解码器用*广义哈密顿动力学*（GHD）：从 $z_0 = (zeta, v_1, ..., v_K)$ 出发做
可训练步长的 leapfrog 式更新，末端加一步布朗动力学输出 $x$，
解码密度可解析求值并与 score 网络联合训练。

#figure(
  align(center, image("../fig_c/tmp_c3_edg.png", width: 90%)),
  caption: [EDG 的解码/编码结构（论文 2401.02080）。实线是解码依赖链
  （$z_0$ → GHD 解码 → $x$ → 潜空间 SDE 扩散），
  虚线是编码依赖（score 网络逆时回到 $z_0$）。
  联合损失 $cal(L)(theta, phi)$ 同时训练两支。],
)

实验：2D 能量函数的六个分布（MoG2、MoG2(i)、MoG6、MoG9、Ring、Ring5）
上，EDG 的加权 MMD 在五列并列最低，Ring5 一列以 $0.02$ 对
L2HMC 的 $0.01$ 排第二；
贝叶斯逻辑回归四个数据集精度与 AUC 全部第一；
Lennard-Jones 原子簇（LJ13 维数 39、LJ55 维数 165）的原子间距分布
最贴合参考数据，重加权有效样本比 #term("rESS", [$ (sum_i w_i)^2 / (N sum_i w_i^2) $，衡量重要性权重的均匀程度，1 为完美、趋 0 为退化]) 比 Boltzmann Generator 与 PIS 高一个量级；Ising 模型（连续松弛）的 $log Z$ 估计在多数温度下最优。

#figure(
  align(center, image("../fig_c/c3_lj.png", width: 86%)),
  caption: [LJ13 与 LJ55 的原子间距分布（论文 2401.02080）：
  Boltzmann Generator（左列）与 EDG（右列）对比参考数据（灰）。
  维数 165 时 EDG 仍贴合双峰结构。],
)

= 后继二：反应坐标的最优降维（FMRC，2024）

== 问题：降维有没有“最优”，能不能端到端学

高维动力系统的长期行为由少数 #term("反应坐标", "reaction coordinate") 刻画
（分子动力学里构象变化的主轴）。基于转移算子谱分解的方法（#term("TICA 与 VAMPnets", "分子动力学里用转移算子的特征函数做降维的方法族，前者线性、后者用两个网络")）给出线性降阶模型，转移密度位于非线性低维流形上时失效；
深度方法（time-lagged 自编码器、信息瓶颈式反应坐标、normalizing flow）
缺少“损失函数与降维引起的算子误差”之间的定量联系。
Bittracher 等（2023）给了两个评估准则（见下），但数值实现
需要估计与积分转移密度，难以扩展。

Zhang–Guo–Wu（2024）把两个准则*等价改写为两个条件概率模型的可计算损失*，
用 flow matching 端到端优化，并证明损失值直接控制算子误差。

== 两个准则

数据是 Markov 过程 $X_t$ 的转移对，lag 为 $tau$，转移密度 $p_tau (x, y)$，
经验联合密度 $rho(x, y) = rho_0(x) p_tau (x, y)$（*不假设平衡*）。
两个算子：Perron–Frobenius 算子 $cal(T)_tau u_t (y) = E[u_t (X_t) | X_(t+tau) = y]$
（密度前传），Koopman 算子 $cal(K)_tau f(x) = E[f(X_(t+tau)) | X_t = x]$
（观测量前传）。

#definition("两个准则")[
  *Lumpability*（可归并性）：$p_tau (x, y) = p_tau^L (r(x), y)$ 对所有 $x, y$。
  即“$x$ 的转移密度只依赖 $r(x)$”：反应坐标是预测未来的充分统计量。

  *Decomposability*（可分解性）：$p_tau (x, y) = p_tau^D (x, r(y)) p_"local" (r(y), y)$。
  即“两步采样”：先由 $p_tau^D$ 定未来落在哪个水平集，
  再由与 $x$ 无关的 $p_"local"$ 在集合内定位。反向等价刻画
  （第 13 章的引理）：decomposable 当且仅当反向转移密度
  $p_(-tau)(y, x) = rho_0(x) p_tau (x, y) \/ rho_1(y)$ 满足
  $p_(-tau)(y, x) = p_(-tau)^D (r(y), x)$。
]

直觉：lumpability 管“过去的哪些信息对未来必要”（正向充分性），
decomposability 管“未来的状态能不能由反应坐标加局部噪声解释”（反向充分性）。
两者合起来才是“这个反应坐标无损”。

== Flow matching 损失与主定理

Flow matching 的元素：路径 $Y^s = (1 - s) Y^0 + s Y^1$，
目标速度场 $v^*(s, y) = E[Y^1 - Y^0 | Y^s = y]$，训练损失
$integral_0^1 E[norm(v(s, Y^s) - (Y^1 - Y^0))^2] d s$。
FMRC 同时恢复 $p_tau$ 与 $p_(-tau)$（两个方向的速度场 $v_0, v_1$），
然后在速度场的第三个参数处塞进瓶颈：

$ min_(theta) cal(L)^"FMRC" = cal(L)_0[v_0^"RC" (dot, *, r(dot))] + cal(L)_1[v_1^"RC" (dot, *, r(dot))] , $

即速度场只能看到 $r(X_t)$ 或 $r(X_(t+tau))$，看不到完整状态。

#proposition("定理（主定理）")[
  对给定反应坐标 $r$：

  $ min_(v_0^"RC", v_1^"RC") cal(L)^"FMRC" >= min_(v_0, v_1) cal(L) , $

  且*等号成立当且仅当 lumpability 与 decomposability 同时满足*。
]

证明要点：等号方向，最优瓶颈速度场等于全场，其诱导的流恰好给出
lumpability 的 $p_tau^L$ 与反向的 $p_(-tau)^D$；反方向，在两个准则下
条件期望允许把 $X_t$ 换成 $r(X_t)$，瓶颈不损失信息。

于是 $cal(L)^"FMRC"$ 的*数值*就是反应坐标质量的度量，
联合最小化 $min_theta cal(L)^"FMRC"[v_(0, theta)^"RC", v_(1, theta)^"RC", r_theta]$
端到端学最优 RC。误差的定量版（定理 5.7）：速度场误差
$epsilon$ 经 Lipschitz 常数指数放大成联合分布的 $W_2$ 误差，
再用 Benamou–Brenier 线性化（$W_2$ 与 $dot(H)^(-1)$ 范数的比较）
转成算子的弱范数误差：

$ norm(cal(K)_tau - hat(cal(K))_L)_H <= sqrt(beta) epsilon e^(integral_0^1 L_s d s) . $

实验（七井圆形势，再经 #term("Swiss roll 映射", "把二维平面卷成三维流形的经典非线性测试变换") 扭曲的数据构造）：
一维 RC 清晰区分七个亚稳井，#term("PCCA+", "把转移算子的特征向量聚类成亚稳态的谱聚类方法") 聚类的簇内 RC 值一致、
簇间有清晰间隙，与理论刻画（lumpable ⇒ 簇内转移同质）一致。

#figure(
  align(center, image("../fig_c/c4_rc.png", width: 84%)),
  caption: [学到的反应坐标（论文 2408.15981）：$r(x)$ 投影在二维
  Swiss roll 数据上，颜色区分七个亚稳井。同一井内的点 RC 值连续，
  井与井之间有清晰分界，这正是 lumpability 的几何形态。],
)

= 后继三：置换不变的算子学习 UQ（UQ-SONet，2025）

== 问题：变长传感器 + 不确定度

DeepONet 的两个限制：传感器数量与位置必须训练测试一致；
没有内置 UQ。VIDON（Prasthofer 等）用 Deep Sets + set transformer
把 branch 换成对变长集合的置换不变编码，但它是确定性映射，
不处理稀疏、噪声或算子本身的随机性。

Ma–Guo–Wu–Zhou（2025）的 UQ-SONet 组合两样东西：#term("set transformer", "处理集合输入的网络：逐元素编码后用注意力做置换不变的汇合，源自 Deep Sets") 嵌入与条件 VAE。

== 结构

*嵌入*：每个传感器读数先嵌入 $Lambda_i = Lambda_x (x_i) + Lambda_kappa (kappa(x_i))$，
然后注意力池化成固定维表示（$H$ 个头）：

$ h^(l)(cal(O)) = sum_(i=1)^m (exp(w_l (Lambda_i) \/ sqrt(d_"emb")))/(sum_(k=1)^m exp(w_l (Lambda_k) \/ sqrt(d_"emb"))) dot v_l (Lambda_i) , $

$h(cal(O))$ 是各头输出的拼接。*置换不变性*来自求和：传感器换顺序，
注意力权重的集合不变，和不变。

*条件 VAE*：潜变量 $z tilde cal(N)(0, I)$ 承载“观测不完备 + 算子固有随机性”
两层信息；解码器是 DeepONet 型 $hat(cal(G))(h(cal(O)), z)(y) = sum_(n=1)^p b^n (h, z) t^n (y)$，
输出各点条件独立的乘积高斯
$p_D (overline(u) | cal(O), z) = product_i cal(N)(u(y_i) | hat(cal(G))(y_i), M sigma_u^2)$。
*人工噪声方差乘 $M$ 是刻意的*：它保证 $M -> infinity$ 时离散损失收敛到
*函数空间 cVAE* 的损失（基于 Seidman 等的 functional VAE 理论，
附录严格推导）。编码器 $q_E (z | cal(O), overline(u))$ 是对角高斯。
损失是负 ELBO：

$ cal(L) = E[ D_"KL" (q_E (z | overline(u), cal(O)) ‖ p(z)) - E_(q_E)[log p_D (overline(u) | cal(O), z)] ] , $

展开后 KL 项有解析式（对角高斯对标准高斯），似然项用重参数化采样。
推断时多次采 $z$、前向、统计条件分布。

#figure(
  align(center, image("../fig_c/tmp_c5_struct.png", width: 92%)),
  caption: [UQ-SONet 总结构（论文 2509.25646）。上半：set transformer
  把变长传感器集合编码成置换不变表示 $h(cal(O))$；
  下半：条件 VAE 的编码器（训练时）与先验采样（推断时），
  DeepONet 型解码器输出均值场。不确定度来自潜变量 $z$ 的多次采样。],
)

== 实验：不确定度随传感器数的两种行为

一维扩散方程（$10^4$ 训练对）：传感器数 $m$ 从 1 增到 10 时
*不确定度带随 $m$ 收缩*，均值误差全面优于 VIDON（如 $m = 1$：
4.47 对 5.35，$times 10^(-2)$）；乘性噪声下误差几乎不变。
二维 Navier–Stokes（涡量-速度形式，$nu = 0.001$）：均值误差 1.89 到 2.69 对
VIDON 的 5.98 到 8.74（$times 10^(-2)$），优势最明显。

第四个实验是*随机椭圆方程*（系数与源项都由 GP 抽样，算子本身随机）：
即使传感器很密，不确定度也不消失。这里要说清模型的立场：
UQ-SONet 没有把 aleatoric 与 epistemic 写成两项分解，两类来源
都塞进同一个潜变量 $z$；能分开的是*观测行为*，即随 $m$ 增大而收缩的那部分
来自观测不完备，收缩不到零的剩余部分对应算子的内在随机性。
比较的基准是条件 GP 后验（每套观测采 1000 个后验样本），
指标是均值与标准差的相对 $L_2$ 误差，不是覆盖率类指标。

#figure(
  align(center, image("../fig_c/tmp_c5_var.png", width: 90%)),
  caption: [不同传感器数下的条件预测（论文 2509.25646）：
  $m = 2, 4, 7, 10$ 的代表样本，参考解（灰）、预测均值（蓝）
  与不确定度带。带宽度随 $m$ 增大系统性收缩，收缩行为本身
  就是“网络知道数据在变多”的证据。],
)

附录还有一个对照：给 VIDON 加深度集成，均值准，
但不确定度带*系统性偏窄*（低估条件不确定性），
再次呼应综述与 IB-UQ 里“集成低估 OOD/条件不确定度”的发现。

= 后继四：物理约束的潜变量模型（LVM-GP，2025）

== 问题：潜变量要不要空间相关

LVM-GP（Feng–Guo–Wan–Wu–Zhou–Zhou，2025）处理含噪数据下的
PDE 正反问题。动机来自一个一维观察：$lambda u_(x x) = f$ 中
$u$ 是 $f$ 的二重积分，因此

$ op("Var")(u(x)) tilde integral integral_(Omega_x) op("Var")(f(xi)) d xi , $

*u 的不确定性是源项不确定性的全域能量聚合*，带空间相关结构。
常规潜变量模型（VAE 类）用低维标准高斯潜变量，
表达不了函数级的空间相关性；B-PINN-HMC 精确但极贵；
深度集成在高噪与数据稀疏下不可靠（产生伪振荡）。

== 方法：把 IB-UQ 的编码器接上 GP 先验

编码器与 IB-UQ 的置信度混合*同构*，差别只在一处：先验 $z_0$
从标准高斯换成*高斯过程*：

$ bold(z)(bold(x); theta_E) = "diag"(m(bold(x))) overline(bold(z))(bold(x)) + "diag"(bold(1) - m(bold(x))) bold(z)_0, quad bold(z)_0 tilde op("GP")(0, K(bold(x), bold(x)')) , $

$K = sigma_K^2 exp(-norm(x - x')^2 \/ (2 l^2)) bold(I)_(d_z)$。
$m -> 1$ 时 $z$ 退化为确定性特征；$m -> 0$（远离数据）时 $z$ 落入 GP 先验，
它带着空间相关性。对固定 $x$，$op("Var")(z(x)) = (1 - m(x))^2$，
置信函数因此是“非平稳核 GP 的另一种参数化”。

解码器是神经算子型的深层积分算子（位置变换核），
条件高斯输出 $u(x) = z_L (x) + sigma_u (x) omega_D$。
物理约束体现在三处：$f$ 与 $b$ 的均值由 PDE 诱导
（$mu_f = cal(N)_x[mu_u]$）；GP 的导数有解析处理（KL 展开逐项求导）；
反问题里 $lambda$ 与 $u$ 共享同一个编码器（一致潜空间）。

损失 = 数据负对数似然 + IB 式正则。正则有两个版本：
基本版是 $E D_"KL" (q_E ‖ e_1)$（$e_1$ 为标准高斯），防止编码器把
$m$ 全推到 1（潜 GP 塌缩成确定性映射）；
改进版对一批随机点取联合 GP 边际 $e_2 = cal(N)(0, bold(K))$，
KL 有解析式，同时优化 $m$ 与 GP 核超参数（实验中长度尺度从 1.0
自适应收敛到 0.58，先验自己学到了数据的相关长度）。

#figure(
  align(center, image("../fig_c/tmp_c6_fwd.png", width: 92%)),
  caption: [LVM-GP 正问题全流程（论文 2507.22493）：置信度编码器把数据
  映到潜随机场（GP 先验兜底），深层积分算子解码出 $u$ 与 $f$，
  两者的均值由 PDE 耦合，损失由三份数据似然加 IB 正则组成。],
)

实验（对照 B-PINN-HMC 与 20 个深度集成）：一维 Poisson 噪声
$sigma = 0.1$ 时深度集成的 $f$ 不确定度失效，LVM-GP 比 HMC 均值更准、
带更紧；一维非线性 Poisson 反问题识别 $lambda = 0.7$，
深度集成的解出现明显伪振荡而 LVM-GP 平滑；二维扩散反应
$sigma = 0.1$ 时 LVM-GP 的标准差 $5.70 times 10^(-3)$ 对 HMC 的
$4.98 times 10^(-2)$，*预测标准差紧一个量级且均值更准*；
六维污染源反演定位三个高斯源，三个定位均值都接近真值，
但 LVM-GP 报告的标准差比 HMC 大一个量级（$10^(-2)$ 对 $10^(-3)$），
带宽的取舍在下一节按偏差与标准差之比逐项看。
训练策略有一条实用经验：先 5000 步只训均值、再加 std 微调，
分阶段训练是稳定性的关键。
