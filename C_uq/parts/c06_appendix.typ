#import "../../preamble.typ": *

= 附录

== A. 记号表

#table(
  columns: (auto, auto),
  inset: 6pt,
  align: (left, left),
  stroke: 0.4pt + luma(190),
  table.header([*记号*], [*含义*]),
  [$theta$], [网络权重],
  [$cal(D)$], [训练数据],
  [$p(theta | cal(D))$], [权重后验],
  [$hat(theta)_j$], [后验的第 $j$ 个样本（集成成员 / MCMC 样本）],
  [$overline(mu)(x), overline(sigma)^2(x)$], [预测均值与总方差],
  [$Sigma_u^2$], [噪声方差（aleatoric）],
  [$cal(F), cal(B)$], [PDE 算子与边界算子],
  [$Z, X, Y$], [IB 的表示变量、输入、目标],
  [$tilde(X), tilde(p)$], [加温输入与其分布（$prop p^(1\/tau)$）],
  [$q_E, q_D, e$], [编码器、解码器、边际密度代理（flow）],
  [$m(x), overline(z)(x)$], [置信度门控与确定性特征],
  [$beta, tau$], [IB 权衡系数与加温温度],
  [$U(x), Z$], [Boltzmann 势能与配分函数（第 7 章语境）],
  [$p_tau, p_(-tau)$], [正向与反向转移密度],
  [$cal(K)_tau, cal(T)_tau$], [Koopman 算子与 Perron–Frobenius 算子],
  [$r(x)$], [反应坐标],
  [$h(cal(O))$], [传感器集合的置换不变表示],
  [$z$], [cVAE / LVM 的潜变量],
  [$K(x, x')$], [GP 核函数],
)

== B. 关键公式速查

#boxed(title: "基础")[
  - *BMA*：$p(u|x, cal(D)) = integral p(u|x, theta) p(theta|cal(D)) d theta$
  - *MC 近似*：$overline(p)(u|x) = 1/M sum_j p(u|x, hat(theta)_j)$
  - *方差分解*：$overline(sigma)^2 = overline(sigma)_a^2 + overline(sigma)_e^2$
    （噪声项 + 成员散度项）
  - *校准*：$P(U_X <= overline(P)_X^(-1)(p)) = p$ 对所有 $p$；
    误差 RMSCE $= sqrt(E_p (p - hat(p)(p))^2)$
]

#boxed(title: "IB-UQ", fill: rgb("#f2fbf5"), stroke: rgb("#8dc9a3"))[
  - *目标*：$cal(L)_"IB" = I(Z;Y) - beta I(tilde(Z); tilde(X))$，
    $tilde(p) prop p^(1\/tau)$
  - *编码器*：$z = "diag"(m(x)) overline(z)(x) + "diag"(1 - m(x)) z_0$
  - *下界*：$I(Z;Y) - cal(H)(Y) >= E[log q_D (y|z)]$
  - *上界*：$I(tilde(Z);tilde(X)) <= E[log q_E (tilde(z)|tilde(x)) - log e(tilde(z))]$
  - *OOD 命题*：$p(x) <= epsilon tilde(p)(x) ==>
    op("KL")(q_E^*(dot|x) ‖ e^*) <= epsilon beta^(-1) log C$
  - *算子版*：瓶颈在 branch，
    $cal(L)_"IBONet" = I((Y,Z); s(Y)) - beta I(tilde(Z); tilde(O))$
]

#boxed(title: "四个后继的核心量", fill: rgb("#fff9ec"), stroke: rgb("#ecc76a"))[
  - *EDG*：$cal(L) = E[log p_D (x|z_0) + U(x)] + integral g^2\/2 E[norm(s)^2 + 2 div s]$，
    训练不需要积分 SDE；权重 $w = e^(-U) p_E (z_0|x) \/ (p_D (z_0) p_D (x|z_0))$
  - *FMRC*：lumpability $p_tau (x,y) = p_tau^L (r(x), y)$；
    decomposability $p_tau (x,y) = p_tau^D (x, r(y)) p_"local" (r(y), y)$；
    $cal(L)^"FMRC"$ 取等 ⟺ 两准则成立
  - *UQ-SONet*：注意力池化 $h^(l) = sum_i "softmax"(w_l (Lambda_i)\/sqrt(d)) v_l (Lambda_i)$；
    噪声方差乘 $M$ 使损失收敛到函数空间 cVAE
  - *LVM-GP*：$z = "diag"(m) overline(z) + "diag"(1 - m) z_0$，
    $z_0 tilde op("GP")(0, K)$；$op("Var")(z(x)) = (1 - m(x))^2$
]

== C. 六篇论文信息

- K. Psaros, X. Meng, Z. Zou, L. Guo, G. E. Karniadakis.
  *Uncertainty quantification in scientific machine learning: methods,
  metrics, and comparisons.* SIAM Review, 65(2), 2023. arXiv:2201.07766.
- L. Guo, H. Wu, W. Zhou, Y. Wang, T. Zhou. *Information bottleneck based
  uncertainty quantification for neural function regression and operator
  learning.* J. Comput. Phys., 510:113089, 2024. arXiv:2302.03271.
- Y. Wang, L. Guo, H. Wu, T. Zhou. *Energy based diffusion generator for
  efficient sampling of Boltzmann distributions.* J. Comput. Phys., 2024.
  arXiv:2401.02080.
- Z. Zhang, L. Guo, H. Wu. *Optimal low-dimensional approximation of
  transfer operators via flow matching: computation and error analysis.*
  2024. arXiv:2408.15981.
- L. Ma, L. Guo, H. Wu, T. Zhou. *Deep set based operator learning with
  uncertainty quantification.* 2025. arXiv:2509.25646.
- X. Feng, L. Guo, X. Wan, H. Wu, T. Zhou, W. Zhou. *LVM-GP:
  Uncertainty-aware PDE solver via coupling latent variable model and
  Gaussian process.* 2025. arXiv:2507.22493.

#boxed(title: "本笔记读到的版本")[
  笔记里的公式、定理编号与数值一律按 arXiv TeX 源码抄，源码与期刊正稿有几处差别，
  涉及结论时以源码为准并在此说明：

  - 2201.07766 的源码是投 Elsevier 的 preprint 版（elsarticle 模板）。
    它定义的指标是 RL2E、MPL、RMSCE、PIW、SDCV 与两个金标准指标 NIP-G、KL-G；
    正式发表版里常见的 CRPS、PICP、ACE 在这一版源码中没有出现，
    笔记提到它们时会说明出处不在这里。
  - 2401.02080 的源码是 2025 年改投 Neural Networks 的重投稿，管线图的日期戳是 2025-09-09。
    该版只有一个编号定理（KL 上界），正文里没有 MMD 收敛速率界；
    JCP 正稿中的相应内容在这份源码中找不到，笔记只引用这份源码里存在的结果。
  - 2408.15981 的源码是 SIAM 会议短文格式，只有一个实验、没有数值表格，
    引用的补充材料文件不在源码包内；$L^2_(mu)$ 误差与谱误差的一节被整段注释掉。
  - 2509.25646 与 2507.22493 的源码都是 elsarticle preprint，
    两者都没有编号定理（前者只有 Jensen 型 ELBO 与函数空间损失的推导，
    后者只定义了 theorem 环境却未使用）。
]

=== 背景文献

- N. Tishby, F. C. Pereira, W. Bialek. *The information bottleneck method.*
  1999. IB 的原始论文。
- A. Alemi, I. Fischer, J. Dillon, K. Murphy. *Deep variational information
  bottleneck.* ICLR, 2017. VIB 变分下界的来源。
- A. Alemi, I. Fischer, J. Gilmer. *Deep variational information bottleneck
  and uncertainty.* 2018. “VIB 蕴含不确定性”的动机文献。
- M. Raissi, P. Perdikaris, G. E. Karniadakis. *Physics-informed neural
  networks.* J. Comput. Phys., 2019. PINN。
- L. Lu, P. Jin, G. E. Karniadakis. *DeepONet: Learning nonlinear operators.*
  2021. 算子学习基座。
- L. Dinh, J. Sohl-Dickstein, S. Bengio. *Density estimation using RealNVP.*
  arXiv:1605.08803, 2016；与 P. Sorrenson, C. Rother, U. Köthe.
  *Disentanglement by nonlinear ICA with general incompressible-flow
  networks (GIN).* ICLR, 2020（保体积 flow，IB-UQ 用它给输入分布“加温”）。
- J. Sohl-Dickstein, E. Weiss, N. Maheswaranathan, S. Ganguli.
  *Deep unsupervised learning using nonequilibrium thermodynamics.* 2015；
  与 J. Ho, A. Jain, P. Abbeel. *Denoising diffusion probabilistic models.*
  2020. 扩散模型基础。
- M. R. Prasthofer, T. De Ryck, S. Mishra. *Variable-input DeepONet
  (VIDON).* 2022. 置换不变算子学习的前身。
- J. H. Seidman, G. Kissas, G. J. Pappas, P. Perdikaris. *Variational
  autoencoding neural operators.* ICML, PMLR 202 (2023). 函数空间 VAE 理论依据
  （UQ-SONet 附录用它的定理 4.1 推离散无关的泛函损失）。
- F. Noé, S. Olsson, C. Köhler, H. Wu. *Boltzmann generators.* Science,
  2019. 采样对照基线。

== D. 自检清单

#boxed(title: "读完后应该能回答")[
  1. aleatoric 与 epistemic 的区分标准是什么？方差分解定律怎么把
     两者算出来？（第 1 章）
  2. “过置信”的定义？为什么单指标选型会误导？（第 1、2 章）
  3. BMA 的 MC 近似里，成员散度对应哪一类不确定性？（第 1 章）
  4. IB-UQ 把经典 IB 的压缩项换成了什么？为什么宽分布压缩能带来
     OOD 鲁棒？（第 5 章）
  5. 置信度混合编码器的两个极端（$m = bold(1)$ 与 $m = bold(0)$）
     各对应什么行为？（第 5 章）
  6. OOD 退化命题的条件与结论分别是什么？（第 5 章）
  7. EDG 为什么不能沿用普通生成模型的损失？simulation-free 训练
     省掉了什么？（第 7 章）
  8. lumpability 与 decomposability 各自刻画反应坐标的什么性质？
     FMRC 的损失为什么能同时检验两者？（第 8 章）
  9. UQ-SONet 里“乘 $M$ 的噪声方差”解决什么问题？（第 9 章）
  10. LVM-GP 相对 IB-UQ 把先验换成了什么？动机来自哪个一维观察？
      （第 10 章）
]
