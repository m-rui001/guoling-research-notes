#import "../../preamble.typ": *

= 第五篇：分布未知时的数据驱动基（2018）

== 换数据假设：认知不确定性

前四篇都假设随机变量的分布*已知*（至少知道它的密度类型）。
但工程里常见的情形是*认知不确定性*（epistemic）：只有一批样本
（或只有若干阶矩），密度本身未知。这击中了前面所有方法的要害：
$ell_1$ 的采样测度需要密度来定义；正交基（gPC 的 Wiener–Askey 对应）
需要分布来构造。

Guo–Liu–Zhou（2018）的两步法把两件事解耦：

1. *构基*：只从样本的矩出发构造正交多项式（data-driven aPC），
   与任何采样假设无关，可以离线做；
2. *采样与拟合*：Christoffel 加权最小二乘 + 平衡测度采样，
   采样测度只依赖“支撑是否有界”这类粗信息，不依赖具体密度。

== 构基：Hankel 矩阵与 Cholesky

样本矩 $mu_k = 1/M sum_(m=1)^M eta_m^k$。正交多项式的三项递推系数
由矩匹配方程决定，把它整理成 *Hankel 矩阵*

$ bold(H) = mat(mu_0 & mu_1 & dots.c & mu_k; mu_1 & mu_2 & dots.c & mu_(k+1); dots.v & & & dots.v; mu_k & mu_(k+1) & dots.c & mu_(2k)) , $

要求它在 Hamburger 意义下正定（$det bold(H) > 0$，即矩序列对应一个真实测度）。
Cholesky 分解 $bold(H) = bold(R)^T bold(R)$ 后，三项递推系数直接从 $R$ 读出
（Mysovskih 定理 + Golub–Welsch 公式）：

$ a_j = r_(j, j+1)\/r_(j, j) - r_(j-1, j)\/r_(j-1, j-1) , quad b_j = r_(j+1, j+1)\/r_(j, j) . $

#insight("为什么这条路稳")[
  直接解矩匹配方程（经典 aPC 的做法）等价于求一个高次 Vandermonde 矩阵的逆，
  次数一高就病态。Hankel + Cholesky 把同样的问题变成对*条件良好的对称正定
  矩阵*做分解：数值上稳定，而且“矩是否对应合法测度”的判定
  （Hamburger 正定性）恰好就是分解能否完成的判定，算法与理论合一。
]

== 采样与拟合：平衡测度 + Christoffel 加权

拟合用加权最小二乘：

$ f_N = arg min_(p in P_N) 1/M sum_(m=1)^M w_m (p(z_m) - f(z_m))^2 , quad w_m = N / sum_(j=1)^N Phi_j^2(z_m) , $

权重正是 Christoffel 函数。为什么这么加权？来自 Cohen 等的稳定性定理：

#proposition("定理（Cohen–LS 稳定性）")[
  记 $kappa(N) = max_z sum_(j=1)^N Phi_j^2(z)$。若

  $ kappa(N) <= delta M / log M , quad delta = (1 - log 2)/(2 - 2 r) , $

  则无权最小二乘以概率 $1 - 2 M^(-r)$ 稳定。
]

Legendre 基的 $kappa(N) tilde N^2$，要求 $M >> N^2$，不划算。
改用加权基 $hat(Phi)_j = Phi_j \/ sqrt(K(z))$ 后 $hat(kappa)(N) = N$
（*最优*）。严格无偏要求按 $kappa(z) rho(z) = rho(z)\/K(z)$ 采样
（采样偏向能量大的点），而位势论给出
$kappa rho ->$ 平衡测度（$N -> infinity$），于是用平衡测度采样：
有界域上是张量积 Chebyshev 测度（$z_m = cos(u_m)$，*对任意容许输入密度普适*）；
无界域上高斯型用密度 $C (2 - norm(z)^2)^(d\/2)$（球内，Beta 分布采样半径）、
指数型用 $C sqrt((4 - sum_i zeta_i)^d \/ product_i zeta_i)$（Dirichlet 向量法）。

#warn("这篇的边界条件")[
  三个诚实的局限（结论节自述）：离散分布情形的稠密性定理仍是开放问题；
  样本矩本身带估计误差（density error）未量化；无界域需要假设密度类型
  （高斯型还是指数型），否则只能截断。*采样测度与构基解耦*换来的稳定性
  是“准线性”（$M = C N log N$），没有 $ell_1$ 的欠定能力，目标真的极稀疏、样本极缺时仍需要下一章的合流版本。
]

实验覆盖了四组混合分布（二项、Poisson、均匀、正态的两两组合，含离散）、
直方图原始数据、5 维情形与电阻网络（$d = 2 p$ 个不确定电阻的节点电压），
对数线性采样率 $M = N log N$ 下条件数基本有界。

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a5_cond.png"), image("../fig_a/tmp_a5_resistor.png"))),
  caption: [左：2 维混合分布（二项 × Poisson）下的条件数，对数线性采样率
  （log-linear）曲线平稳、线性采样率（$1.5N, 2N$）随 $N$ 缓慢爬升。
  右：电阻网络（$d = 4$，各电阻服从不同指数分布），节点电压的逼近误差
  随采样率收敛。分布可以任意混合（包括离散），这是数据驱动基相对分布假设方法的主要优势。],
)

= 第六篇：合流（induced sampling，2020）

== 问题：两条路线各有半截

到这里有半截 $ell_1$ 路线（样本少、欠定、但采样测度要密度）和半截
数据驱动路线（密度无关、但只到超定最小二乘）。
Guo–Narayan–Liu–Zhou（2020）把它们接上：*data-driven 基 + 预条件 $ell_1$ +
从 induced 测度采样*。这是六篇里集成度最高的一篇。

== 诱导测度：只依赖样本的采样分布

构基与第五篇相同（加权矩 → Hankel → Cholesky → 三项递推；
多元按坐标边际化逐维构基再张量化，注意所得基只尊重边际测度）。
关键是采样。定义*诱导测度（induced measure）*：在经验样本集
$Xi = {z^((1)), ..., z^((Q))}$ 上

$ mu(z) = sum_(j=1)^Q tilde(kappa)_j delta_(z^((j))) , quad tilde(kappa)_j = kappa(z^((j))) / sum_(q=1)^Q kappa(z^((q))) , $

其中 $kappa(z; Lambda_K^"TD") = 1/N sum_(j=1)^N Phi_j^2(z)$ 是归一化
Christoffel 函数。读法有两层，方向容易搞反。这篇论文把 $kappa$ 定义为
*能量*的平均（$1\/N sum Phi_j^2$），所以 $kappa$ 大的点（多项式能量高）
被采得更多；随后每个样本按 $1\/kappa$ 加权，两个因子相乘恰好还原经验分布，
只差一个全局常数，被回归问题吸收。这就是原文所说“用 Radon–Nikodym
导数 $1\/kappa$ 保持无偏”的准确含义。偏向高能量点采样的作用：让
$sqrt(W) A$ 的最大元有一致界：对任何点都有
$abs(phi_k (z))^2 <= sum_j Phi_j^2(z) = N kappa(z)$，于是加权后
$abs(phi_k (z))\/sqrt(kappa(z)) <= sqrt(N)$，与采样位置无关。
支撑就是样本集本身，有限点上的逆变换采样，实现只有几行。

预条件权取 $W_j = 1\/kappa(z^((j)))$，基追踪格式
$min norm(c)_1$ s.t. $sqrt(W) A c = sqrt(W) b$。

#insight("这个组合为什么自洽")[
  压缩感知的恢复误差理论罚项正比于 $sqrt(W) A$ 的最大元，
  $1\/kappa$ 加权把最大元压平。加权等效于测度偏置，要保持估计无偏
  就必须从“偏置后的测度”采样，而偏置后的测度恰好是 $omega\/kappa$ 的
  离散版本，即 induced 测度。采样、加权、无偏三者咬合：
  *从 $mu$ 采样、按 $1\/kappa$ 加权，等于从原分布做了一次重要性采样*。

  与平衡测度的关系需要说清：$omega$ 已知时 induced 测度有唯一极限
  （加权多位势平衡测度；有界超立方上就是乘积 Chebyshev 测度）。
  有限 $K$ 时 induced 不是平衡测度，但它*只依赖样本*，
  而且样本保证落在数据支撑内，平衡测度或 Chebyshev 采样会把点
  落在低概率甚至数据集之外的区域。
]

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a6_cas.png"), image("../fig_a/tmp_a6_ind.png"))),
  caption: [两种采样测度的散点对比（论文 2008.10121，数据为
  二项 × Poisson 模拟）。左：CSA 采样（平衡测度的经验版），点撒满了
  Chebyshev 网格，大量点落在数据支撑之外；右：induced 采样，
  点全部落在数据样本上，采样概率按 $kappa$（多项式能量）偏置，再由 $1\/kappa$ 权重补偿，两者配对保持估计无偏。],
)

实验：稀疏恢复（真实密度是均匀、截断正态、截断对数正态的等权混合，
$Q = 10^5$ 个经验样本）：低维高次（$d = 2, K = 20$）induced 显著优于 MC 与 CSA；
高维低次（$d = 10, K = 3$）induced 与 MC 相当
（此时 $mu approx omega$，两种采样没有差别，作者如实说明了这一点）。
Kirchhoff 板弯曲问题（随机杨氏模量 KL 展开，Morley 非协调有限元 1369 个未知数）
上 induced 在 $d = 2$ 显著省样本。

#figure(
  align(center, image("../fig_a/tmp_a6_pde.png", width: 76%)),
  caption: [Kirchhoff 板问题的 $ell_2$ 误差（$d = 2$）：induced 采样（红色）
  在全部样本量上都低于 MC 与 CSA，达到同精度所需样本最少。],
)

== 六篇合起来：Christoffel 函数的六副面孔

#ttable(
  ([论文], [Christoffel 函数的身份], [解决的核心困难]),
  (
    ([第一篇 (2017)], [求积权 + 逐点一致界（$L$ 因子）], [任意分布下的 RIP 样本量]),
    ([第二篇 (2017)], [点列加权（$1\/sqrt(K)$ 进基）], [超定 LS 的确定性选点]),
    ([第三篇 (2018)], [列归一（$P$ 矩阵）], [导数行的各向同性]),
    ([第四篇 (2018)], [（不直接出现）], [罚函数与采样正交改进]),
    ([第五篇 (2018)], [LS 权重 + 平衡测度], [分布未知时的稳定拟合]),
    ([第六篇 (2020)], [采样概率 + 预条件权], [数据支撑上的稀疏恢复]),
  ),
  widths: (auto, 1fr, 1fr),
)

一条线收束：*凡是“多项式在某点太活跃”导致的问题（采样噪声、病态、
权重失衡），答案都是按 Christoffel 函数的倒数加权或偏置*。
