#import "../common.typ": *

= 附录

== A. 记号表

#table(
  columns: (auto, auto),
  inset: 6pt,
  align: (left, left),
  stroke: 0.4pt + luma(190),
  table.header([*记号*], [*含义*]),
  [$x in M subset.eq RR^n$], [动力学系统的状态；$n$ 是状态维数],
  [$T: M -> M$], [离散时间非线性映射，$x_t = T(x_(t-1))$],
  [$H = L^2(M, mu)$], [可观测函数构成的希尔伯特空间],
  [$phi, psi_i$], [观测函数 / 特征（observable）],
  [$z in cal(Z)$], [潜在变量；$cal(Z)$ 是由 ${phi_1, ..., phi_d}$ 张成的潜在空间],
  [$d$], [潜在维数],
  [$cal(K)$], [Koopman 算子；矩阵形式记作 $cal(K)$ 或 $K$],
  [$cal(C)$], [潜在协方差矩阵 $op("Cov")(z)$],
  [$M_n$], [$n$ 步线性前向协方差 $sum_(i=0)^(n-1) cal(K)^i Sigma (cal(K)^i)^T$],
  [$I(a; b)$], [互信息],
  [$D_"KL"(p ‖ q)$], [KL 散度],
  [$H(x)$], [（微分）熵],
  [$S(rho)$], [von Neumann 熵],
  [$rho$], [密度矩阵（半正定、迹为 1）],
  [$lambda$], [特征值（Koopman 算子或协方差矩阵的，视上下文）],
  [$gamma_i$（在注水问题里记作 $g_i$）], [谱增益],
  [$p_i$], [注水问题里的谱权],
  [$alpha, beta, gamma$], [拉格朗日乘子（时间相干 / 结构一致性 / 谱熵）],
  [$theta, omega, psi$], [编码器 / 解码器 / 潜在动力学的网络参数],
  [$t$], [时间步；$n$ 为滞后步数],
  [$cal(N)$], [高斯分布],
  [$op("tr")$], [矩阵的迹],
)

== B. 关键公式速查

#boxed(title: "定义与架构")[
  - *Koopman 算子*：$(cal(K)phi)(x) = phi(T(x))$
  - *不变子空间条件*：$cal(K) Psi(x) = K^T Psi(x)$
  - *EDMD 估计*：$K^T = hat(Y)^T hat(X)(hat(X)^T hat(X))^dagger$
  - *Koopman 自编码器*：$x ->^f z ->^(cal(K)) z' ->^g hat(x)$
  - *LRA 多步一致性*：$cal(L) = EE[(...)"多步重构" + beta (...)"潜在一致性"] + Omega$
  - *VAMP-2 分数*：$hat(R)_2 = norm(cal(C)_(0 0)^(-1\/2) cal(C)_(0 1) cal(C)_(1 1)^(-1\/2))_F^2$
  - *CKO 组合算子*：$cal(K) = sigma times hat(K)$
]

#boxed(title: "信息论工具", fill: rgb("#f2fbf5"), stroke: rgb("#8dc9a3"))[
  - *互信息*：$I(x;y) = EE[log p(x,y)/(p(x)p(y))] = H(x) - H(x\|y)$
  - *链式法则*：$I(x;y) = I(x;(y,z)) - I(x;z\|y)$
  - *DPI*：$X -> Y -> Z ==> I(X;Z) <= I(X;Y)$
  - *IB 拉格朗日*：$cal(L) = I(x;z) - beta I(z;y)$，等价于 $I(x;z) + beta I(x;y\|z) + "const"$
  - *VIB 码率上界*：$I(x;z) <= EE D_"KL"(q(z|x) ‖ r(z))$
  - *VIB 相关性下界*：$I(z;y) >= H(y) + EE log q(y|z)$
  - *InfoNCE*：$hat(I) = log K - cal(L)_"NCE" <= I(x;c)$
  - *von Neumann 熵*：$S(rho) = -op("tr")(rho log rho) = H(lambda(rho))$
  - *有效维度*：$d_"eff" = exp(S(rho))$
]

#boxed(title: "五个命题", fill: rgb("#fff9ec"), stroke: rgb("#ecc76a"))[
  - *命题 1*：$I(x_(n-1);x_n) >= I(z_(n-1);x_n) >= I(z_(n-1);z_n)$
  - *命题 2*：$norm(p - q_"KR")_"TV" <= 1/2 [ sum_(n=1)^t (I(x_(n-1);x_n) - I(z_(n-1);z_n)) + cal(E) ]$，
    且 $EE norm(x_(1:t) - EE_(q_"KR")[x_(1:t)|x_0])^2 >= C exp( - (2)/(n t) sum_n I(z_(n-1);z_n) )$
  - *命题 3*：$I(z_t;x_t) = I(z_(t-n);z_t) + I(z_t;x_(t-1)\|z_(t-n)) + I(z_t;x_t\|x_(t-1))$
  - *命题 4*：注水解 $p_i = max(0, 1/(2 mu) - 1/g_i)$ $->$ 谱倾斜时低秩 $->$ 模式坍缩
  - *命题 5*：加熵正则后所有 $p_i > 0$ $->$ 有效维度被证明提升
]

#boxed(title: "拉格朗日量与损失")[
  - *拉格朗日量*：$max_z alpha log I(z_(t-n); z_t) - beta I(z_t;x_t\|z_(t-n)) + gamma S(cal(C)\/op("tr")(cal(C))) + log p(x_t\|z_t)$
  - *结构一致性等价*：$I(z_n;x_n\|z_(n-1)) = EE[log q_psi(z_n\|z_(n-1))] - H_(p_theta)(z_n\|x_n)$
  - *AE 总损失*：$cal(L) = cal(L)_"rec" - alpha I(z_n;cal(P)_n) + beta norm(z_(n+1) - cal(K) z_n)^2 - gamma S(P)$
]

== C. 补充：命题 2 的率失真下界

命题 2 的第四个证明步骤值得单独展开，因为它把“信息缺口”换算成了“误差的量级”。

#proof[
  *出发点：Shannon 下界。* 对任意源 $X$ 与 $L^2$ 失真，率失真函数满足

  $ R(D) >= h(X) - 1/2 log(2 pi e D) . $

  移项（取自然对数、以 $D$ 为主元）：

  $ D >= 1/(2 pi e) e^(2 h(X) - 2 R) . $

  这就是“误差随码率指数衰减”的一般形式。

  *应用到 Koopman 轨迹。* 把 $x_(1:t) in RR^(n times t)$ 看成 $n t$ 维源（$n$ 是状态维数）。
  码率 $R$ 由潜在信道提供，而潜在信道每一步可传的信息量正是 $I(z_(n-1); z_n)$，
  所以单位维度的码率近似为

  $ R approx 1/(n t) sum_(n=1)^t I(z_(n-1); z_n) . $

  代入 Shannon 下界：

  $ D >= (n t)/(2 e) exp( 2/(n t) h(x_(1:t)) ) exp( - 2/(n t) sum_(n=1)^t I(z_(n-1);z_n) ) . $

  *吸收常数。* 由 $sum_(n=1)^t I(x_(n-1);x_n)$ 决定了 $h(x_(1:t))$ 中与动力学相关的部分
  （剩余部分是加性常数），把这两项一起收进常数

  $ C := (n t)/(2 e) exp( 2/(n t) sum_(n=1)^t I(x_(n-1); x_n) ) , $

  得到

  $ EE_p[ norm(x_(1:t) - EE_(q_"KR")[x_(1:t)|x_0])^2 | x_0 ] >= C exp( - ( 2/(n t) sum_(n=1)^t I(z_(n-1);z_n) + cal(E) ) ) . quad square $
]

#insight("这个下界的三种读法")[
  1. *误差随保留的信息指数下降*。多保留 1 nat，误差降 $e^2$ 倍（约 $7.4$ 倍）。
  2. *常数 $C$ 由原系统的熵决定*。这意味着一个“更难”的系统（熵更高、动力学耦合更强）
     有更高的误差地板，这与 Koopman 表示无关，它是问题本身的难度。
  3. *加上了 $cal(E)$*。近似误差的三个部分（编码、转移、解码）也进入了指数，
     这提醒我们“完美的信息保留”并不能弥补糟糕的参数化。
]

== D. 补充：结构一致性的两个效应

§6 给出了 $I(z_n;x_n|z_(n-1)) = EE[log q_psi] - H_(p_theta)$ 的等价性。
它的两个效应值得分开看。

*(1) 与 Koopman 动力学对齐。* 期望项 $EE_(p_theta(z_n|x_n))[-log q_psi(z_n|z_(n-1))]$
要求从编码器采样的点落在 Koopman 先验的高似然区域。因为先验被参数化为线性高斯转移，

$ q_psi(z_n|z_(n-1)) = cal(N)(z_n | cal(K) z_(n-1), Sigma) , $

最小化这项 KL 迫使编码器输出可以在*线性结构*下被预测。

换一种说法：如果你把编码器输出的点画出来，它们必须“看起来像是从一个线性系统采出来的”。
这就是结构一致性的几何含义。

*(2) 熵正则。* $H_(p_theta)(z_n|x_n)$ 鼓励编码器保持随机性。

#warn("这一项为什么必要：确定性编码器会让互信息退化")[
  如果编码器是确定性的（$z_n = f_theta(x_n)$），则 $H_(p_theta)(z_n|x_n) = -infinity$，
  于是 $I(z_n;x_n|z_(n-1))$ 无定义，整个损失函数崩溃。

  在代码层面这就是“后验坍缩”的反面问题：编码器方差趋零。
  加上 $+H$ 项（在最大化形式里）就是在防止这件事。
]

两者合起来：

$ p(z_n | x_n) approx q(z_n | z_(n-1)) quad <==> quad z_n "在" cal(K) "下近似线性演化" , $

这正是 Koopman 表示的定义要求，被写成了一个可微的损失项。

== E. 参考文献与阅读顺序

=== 核心论文

- X. Cheng, W. Yuan, Y. Yang, Y. Zhang, S. Cheng, Y. He, Z. Sun.
  *Information Shapes Koopman Representation.* ICLR 2026. arXiv:2510.13025.
  代码：`github.com/Wenxuan52/InformationKoopman`

=== Koopman 理论基础

- B. O. Koopman. *Hamiltonian systems and transformation in Hilbert space.* PNAS, 1931.
  Koopman 算子的原始定义。
- B. O. Koopman, J. von Neumann. *Dynamical systems of continuous spectra.* PNAS, 1932.
  连续谱现象。理解“为什么可能根本没有特征值”的关键。
- M. O. Williams, I. G. Kevrekidis, C. W. Rowley. *A data-driven approximation of the
  Koopman operator: extending dynamic mode decomposition.* J. Nonlinear Science, 2015.
  EDMD 的奠基工作。
- S. L. Brunton, B. W. Brunton, J. L. Proctor, J. N. Kutz. *Koopman invariant subspaces and
  finite linear representations of nonlinear dynamical systems for control.* PLoS ONE, 2016.
  证明了多吸引子系统不存在包含状态的有限维线性表示（§1 引用的 Brunton 定理）。
- N. Takeishi, Y. Kawahara, T. Yairi. *Learning Koopman invariant subspaces for DMD.* NeurIPS 2017.
  第一次用网络学观测函数。
- S. L. Brunton, M. Budisic, E. Kaiser, J. N. Kutz. *Modern Koopman theory for dynamical
  systems.* SIAM Review, 2022. 综述，适合作为坐标系。

=== 深度 Koopman 学习

- S. E. Otto, C. W. Rowley. *Linearly recurrent autoencoder networks for learning dynamics.*
  SIAM J. Appl. Dyn. Syst., 2019. LRA，本文架构的祖先。
- B. Lusch, J. N. Kutz, S. L. Brunton. *Deep learning for universal linear embeddings of
  nonlinear dynamics.* Nature Communications, 2018. 同一时期的平行工作，同样值得读。
- S. Pan, K. Duraisamy. *Physics-informed probabilistic learning of linear embeddings of
  nonlinear dynamics with guaranteed stability.* SIAM J. Appl. Dyn. Syst., 2020.
  结构参数化分支。
- O. Azencot, N. B. Erichson, V. Lin, M. Mahoney. *Forecasting sequential data using
  consistent Koopman autoencoders.* ICML 2020. 一致性分支。
- Y. Li, H. He, J. Wu, D. Katabi, A. Torralba. *Learning compositional Koopman operators
  for model-based control.* ICLR 2020. CKO，图/组合结构分支。
- S. Pan, E. Kaiser, B. M. de Silva, J. N. Kutz, S. L. Brunton. *PyKoopman: A Python package
  for data-driven approximation of the Koopman operator.* 2023. KAE 基线。
- Y. Liu, C. Li, J. Wang, M. Long. *Koopa: Learning non-stationary time series dynamics with
  Koopman predictors.* NeurIPS 2023. 时间序列分支。
- X. Wu, X. Qiu, H. Gao, J. Hu, B. Yang, C. Guo. *$K^2$VAE: A Koopman-Kalman enhanced
  variational autoencoder for probabilistic time series forecasting.* ICML 2025.
  Koopman + Kalman + VAE。
- Y. Xu, et al. *ResKoopNet: Learning Koopman representations for complex dynamics with
  spectral residuals.* ICML 2025. 谱残差分支。
- X. Cheng, Y. He, Y. Yang, X. Xue, S. Cheng, D. Giles, X. Tang, Y. Hu. *Learning chaos
  in a linear way.* ICLR 2025. PFNN，混沌结构分支。

=== 变分、核与统计理论

- A. Mardt, L. Pasquali, H. Wu, F. Noé. *VAMPnets for deep learning of molecular kinetics.*
  Nature Communications, 2018. 变分谱分数，*必读*。
- H. Wu, A. Nüske, F. Paul, S. Klus, P. Koltai, F. Noé. *Variational Koopman models.*
  J. Chem. Phys., 2017. 非平衡数据下的变分 Koopman。
- V. Kostic, P. Novelli, A. Maurer, C. Ciliberto, L. Rosasco, M. Pontil. *Learning dynamical
  systems via Koopman operator regression in RKHS.* NeurIPS 2022. KORR。
- V. Kostic, K. Lounici, P. Novelli, M. Pontil. *Sharp spectral rates for Koopman operator
  learning.* NeurIPS 2023.
- P. Bevanda, M. Beier, A. Lederer, S. Sosnowski, E. Hüllermeier, S. Hirche. *Koopman kernel
  regression.* NeurIPS 2023. KKR，*必读*。
- P. Bevanda, et al. *Koopman-equivariant Gaussian processes.* 2025.
- S. Das, D. Giannakis. *Koopman spectra in reproducing kernel Hilbert spaces.* ACHA, 2020.

=== 信息论与表示学习

- N. Tishby, F. C. Pereira, W. Bialek. *The information bottleneck method.* 1999.
  IB 的原始论文。
- N. Tishby, N. Zaslavsky. *Deep learning and the information bottleneck principle.* ITW 2015.
- A. A. Alemi, I. Fischer, J. V. Dillon, K. Murphy. *Deep variational information bottleneck.*
  ICLR 2017. VIB，两块变分界的来源。
- A. van den Oord, Y. Li, O. Vinyals. *Representation learning with contrastive predictive
  coding.* 2018. InfoNCE 的来源。
- M. Wu, C. Zhuang, M. Mosse, D. Yamins, N. Goodman. *On mutual information in contrastive
  learning for visual representations.* 2020. 本文引用的 InfoNCE 参考。
- M. Federici, P. Forré, R. Tomioka, B. S. Veeling. *Latent representation and simulation of
  Markov processes via time-lagged information bottleneck.* ICLR 2024.
  T-IB，*最接近的替代品，必读*。
- O. Roy, M. Vetterli. *The effective rank: A measure of effective dimensionality.* EUSIPCO 2007.
  有效维度。
- F. Bach. *Information theory with kernel methods.* IEEE Trans. IT, 2022.
  论文点名的未来方向。
- E. Witten. *A mini-introduction to information theory.* 2020. 密度矩阵的入门读物。

=== 神经算子（平行线索）

- L. Lu, P. Jin, G. Pang, Z. Zhang, G. E. Karniadakis. *Learning nonlinear operators via
  DeepONet.* Nature Machine Intelligence, 2021.
- Z. Li, et al. *Fourier neural operator for parametric partial differential equations.* ICLR 2021.
- W. Xiong, X. Huang, Z. Zhang, R. Deng, P. Sun, Y. Tian. *Koopman neural operator as a
  mesh-free solver of non-linear PDEs.* J. Comput. Phys., 2024. KNO。

=== 建议的阅读顺序

#insight("如果你只有两天")[
  *第一天（理解问题）*：
  §0 本文（只读引言与 §3），然后 Brunton et al. 2016 的摘要与结论（不变子空间的根本困难），
  然后 Otto & Rowley 2019 的架构图与损失函数。

  *第二天（理解答案）*：
  本文 §3.1–3.3 加附录 F，然后 Tishby & Zaslavsky 2015（IB 的直觉），
  Alemi et al. 2017（VIB 的界），最后 Federici et al. 2024（T-IB，最接近的替代品）。

  读完这两天，你就有能力判断这篇论文的位置了。
]
