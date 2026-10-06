#import "../../preamble.typ": *

= 第一篇：Gauss 求积网格上的随机子采样（SISC 2017）

== 要解决的问题

压缩感知用到随机 PDE 时，前人的采样策略有两类：直接从原密度 $rho$ 独立采样
（Doostan–Owhadi；Yan–Guo–Xiu），或从 Chebyshev 测度采样配预条件
（Rauhut–Ward 的 Legendre 情形；Tang–Zhou 把“张量积 Gauss 网格子采样”
用在均匀分布 + Legendre 上）。Guo–Narayan–Zhou–Chen（SISC 2017）
把 Tang–Zhou 的结构化采样推广到*几乎全部经典一元分布*：
任意 Beta 权（有界域）与单、双侧指数权（无界域，含 Hermite 与 Laguerre）。

网格子采样比自由采样多两样东西。
第一，点被限制在求积网格上，采样是离散抽取，实现简单且可复现；
第二，求积网格自带一套权重（Christoffel 函数），用它做预条件后
测量矩阵*恰好是正交矩阵的行子集*，分析可以进行到底。

== 构造

对每个维度取 $n_i$ 点 Gauss 求积法则：节点 $z^i_k$，权重 $w^i_k$，
以及 Christoffel 函数

$ lambda^i_n (x^i) = 1 / sum_(k=0)^(n-1) (phi^i_k (x^i))^2 . $

构造离散经验测度 $nu_(bold(n))^i = 1/n sum_(k=1)^n delta_(z^i_k)$，张量化为
$nu_(bold(n)) = times_i nu_(bold(n))^i$。采样即从这批有限网格点里无放回抽 $M$ 个。

关键的预条件加权：定义

$ psi^i_(k, n)(x^i) = sqrt(n lambda^i_n (x^i)) phi^i_k (x^i) , $

即“多项式乘以 Christoffel 权重的平方根”。

#lemma("引理（离散正交系统）")[
  $psi^i_(k, n)$ 在离散测度 $nu^i_n$ 下正交归一：

  $ 1/n sum_(p=1)^n psi^i_(j, n)(z^i_p) psi^i_(k, n)(z^i_p) = delta_(j, k) . $

  多元版本 $psi_(bold(k), bold(n))(z) = product_i psi^i_(k_i, n_i)(z^i)$ 在 $nu_(bold(n))$ 下正交归一。
]

#proof[
  由 Gauss 求积的精确性（$n$ 点法则对次数 $<= 2n - 1$ 的多项式精确）：

  $ 1/n sum_(p=1)^n n lambda^i_n (z^i_p) phi^i_j (z^i_p) phi^i_k (z^i_p) = sum_(p=1)^n w^i_p phi^i_j (z^i_p) phi^i_k (z^i_p) = integral phi^i_j phi^i_k rho^i = delta_(j, k) . $

  第一步把 $n lambda_n$（Christoffel 权重）搬进取权和后恰好等于求积权；
  第二步用求积精确性换成连续积分；第三步用基的正交性。多元版由乘积结构直接张开。$square$
]

这个引理是全篇的引擎：它说明“带权多项式在离散测度下构成一个正交系统”，
于是张量化的加权测量矩阵是*离散正交矩阵（DOM）*，而正交矩阵的行子采样
恰好是 RIP 理论最成熟的对象（无放回子采样仍保持 RIP，Rauhut）。

== 主定理：样本量准则

#proposition("定理（主定理，样本量）")[
  设 $Lambda$ 是大小为 $N$ 的有限指标集，$bold(n)$ 是使
  $Lambda subset.eq Lambda^P_(bold(n) - bold(1))$ 的最小多重指标。
  从测度 $nu_(bold(n))$ 无放回抽 $M$ 个样本，若

  $ M >= L(bold(n)) C_1 s log^3 (s) log (N) , $

  $L(bold(n)) = product_(i=1)^d L_i (n_i)$，则加权 $ell_1$ 优化问题的解 $bold(c)^sharp$ 满足

  $ Pr[ norm(bold(c) - bold(c)^sharp)_2 <= C_2 (sigma_(s, 1)(bold(c)))/(sqrt(s)) ] >= 1 - N^(-gamma log^3 (s)) . $
]

读法：需要的样本数 = 稀疏度 × 分布常数 × 对数因子。
分布常数 $L$ 是唯一随分布与维数变化的部分，它的大小决定这条路在高维下
是否还便宜。以下引理把 $L_i$ 在三类分布下显式估出。

#lemma([引理（三类分布的 $L_i$ 上界）])[
  - *Jacobi 权*（$rho^i prop (1 - x)^gamma (1 + x)^delta$，$gamma, delta >= -1\/2$，
    含 Legendre 与 Chebyshev）：
    $L_i (n) <= C(gamma, delta)$，*与 $n$ 无关*（一致有界系统）；
  - *双侧指数权*（$rho^i prop exp(-abs(x)^alpha)$，$alpha > 3\/2$，含 Hermite）：
    $L_i (n) <= C n^(2\/3)$；
  - *单侧指数权*（$x >= 0$，$alpha > 3\/4$，含 Laguerre）：$L_i (n) <= C n^(2\/3)$。
]

#proof[
  证明链分三步（附录）。
  *第一步，定位节点区间。* 用 Mhaskar–Rahmanov–Saff 数
  $a^W_n = (n sqrt(pi) Gamma(alpha\/2) \/ Gamma(alpha\/2 + 1\/2))^(1\/alpha)$
  确定时间尺度，Levin–Lubinsky 定理保证 Gauss 节点落在
  $[-hat(a)_n^W, hat(a)_n^W]$ 内，$hat(a)_n^W = a_n (1 + c n^(-2\/3))$。
  *第二步，逐点界。* 在该区间上引用加权多项式的一致上界
  $sup abs(psi_(k, n))^2 <= C(alpha) n^(2\/3)$（Jacobi 情形右端是常数：
  有界区间上 Christoffel 权重恰好抵消端点的多项式增长）。
  *第三步，矩阵级。* 满 DOM 矩阵 $bold(A)$（$(A)_(p, q) = sqrt(w_(bold(l)(p))) phi_(bold(k)(q))(bold(z)_(bold(l)(p)))$）
  是正交矩阵；逐项估计

  $ (product_(i=1)^d sqrt(n_i)) abs(A_(p, q)) = product_(i=1)^d abs(psi_(dot, n_i)) <= product_(i=1)^d sqrt(L_i (n_i)) = sqrt(L(bold(n))) , $

  再引用“正交矩阵无放回行子采样仍保持 RIP”的结果完成。$square$
]

#figure(
  align(center, image("../fig_a/tmp_a1_bounds.png", width: 90%)),
  caption: [上界 $L_i$ 的形状（论文 1602.00995）。左：Jacobi 情形的界随参数变化，
  与 $n$ 无关；右：Hermite/Laguerre 情形的 $n^(2\/3)$ 界。有界分布
  （Legendre/Chebyshev）是这条路的最爱：$L$ 是常数，样本量只比最优多
  $C^d$ 的维数代价（数值上接近 $2^d s log^3 s log N$）。],
)

#warn([无界分布的 $n^(2\/3)$ 是本质的])[
  Hermite（高斯输入）情形 $L_i tilde n^(2\/3)$ 不是技术瑕疵：子采样分析策略下
  它已经是这种分析策略下的下确界（sharp），不是技术瑕疵。直觉：无界域上 Gauss 节点越跑越远，
  Christoffel 权重在尾部的增长赶不上多项式的增长。数值实验也确认：
  $d = 10$ 时按原密度随机采样（$L$ 型分析不适用）反而更好。
  这提示“采样测度与多项式系统要配套”，也是后来 induced sampling（第 14 章）的动机之一。
]

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a1_rec.png"), image("../fig_a/tmp_a1_pde.png"))),
  caption: [左：$M = 85$ 时稀疏 Legendre 展开的恢复概率（$d = 2$ 与 $d = 10$），
  Gauss 网格子采样（Gaussian）与其他采样策略的对比，低维高维都不劣或更优。
  右：$d = 3$ 随机椭圆 PDE 的均值 $ell_2$ 误差，Gauss 子采样与高精度
  稀疏网格参考解一致。],
)

= 第二篇：加权近似 Fekete 点（确定性采样）

== 换一条路：超定最小二乘

第一篇处理 $M < N$（欠定 + $ell_1$）。如果样本量 $M$ 能略超 $N$
（$M = 1.05 N$ 这类轻度过采样），可以用普通最小二乘：
求解稳定、代价低、不需要 $ell_1$ 优化器。但最小二乘对*点的位置*敏感：
等距的 Newton–Cotes 点即使数据无穷精确也高度不稳定（Vandermonde 病态），
除非大量过采样。随机路线已知：无权 Legendre + 均匀 MC 需要
$M tilde N^2$（Cohen–Davenport–Leviatan）；Christoffel 加权最小二乘
加平衡测度采样可以达到接近最优的 $M tilde N$（Narayan–Jakeman–Zhou）。

Guo–Narayan–Yan–Zhou（2017）走确定性路线：贪心选点，
选出一组“条件数最优”的点集，并证明加权之后*贪心等于全局最优*。

== 构造：把 Christoffel 权重塞进基

对指标集 $Lambda$（大小 $N$），定义未归一化 Christoffel 函数
$K_Lambda (y) = sum_(alpha in Lambda) psi_alpha^2(y)$，然后定义*加权多项式空间*

$ Q = Q(Lambda) = "span"{ psi_alpha / sqrt(K_Lambda (y)) | alpha in Lambda } . $

直觉：$K_Lambda$ 大的地方（多项式容易取大值的地方），基函数被压低，
于是加权后的基在整个域上“能量均匀”。加权最小二乘为

$ hat(bold(f)) = arg min_(bold(v)) norm(bold(V)(A_M, Q) bold(v) - bold(W) bold(f))_2 , quad (W)_(m, m) = 1\/sqrt(K_Lambda (y_m)) , $

$V$ 是加权 #term("Vandermonde 矩阵", "元素为各基函数在各采样点取值的矩阵")：第 $(j, k)$ 元是第 $k$ 个基在第 $j$ 个点上的值，多项式插值的稳定性由它的条件数决定。四类点集的定义：Fekete 点
$A_N^F = arg max abs(det bold(V))$（行列式最大 = 插值意义下张得最开），
最优条件数点集 $A_N^C = arg min kappa(bold(V))$，以及两者的贪心版
（每次加一个点，让目标最优）。

== 主定理：贪心 = 全局最优

#proposition("定理（nd-optimal）")[
  设 $rho$ 是 $RR^d$ 上的概率密度（不要求紧支撑！），
  $Lambda$ 是大小为 $N$ 的任意指标集。点集 $A_N$ 满足

  $ abs(det bold(V)(A_N, Q)) = 1 quad "当且仅当" quad kappa(bold(V)(A_N, Q)) = 1 . $

  且若最优点集存在，贪心迭代有一条分支能取到它。
]

#proof[
  对每个加权行 $bold(psi)(y_j) = (psi_1(y_j)\/sqrt(K), dots)$ 有恒等式
  $norm(bold(psi)(y_j))_2 = 1$（定义直接验证：分母就是行自身的平方和）。
  由 Hadamard 不等式 $abs(det bold(A)) <= product_i norm(bold(a)_i)_2$ 且各行范数全是 1，
  所以 $abs(det bold(V)) <= 1$，等号当且仅当各行两两正交，
  即 $bold(V)$ 是正交矩阵、条件数 $kappa = 1$。
  贪心部分用归纳：$V^T$ 的 QR 分解给出
  $abs(det(V V^T)) = r^2$，$r = sqrt(1 - norm(bold(b))_2^2) <= 1$，
  等号当且仅当新点与已有的加权向量正交；只要最优集中的点可被贪心逐个达成
  （每次选到与已选集合正交的点），贪心分支就收敛到全局最优。$square$
]

一维情形的结构可以完全写出来。

#proposition("定理（1d-optimal）")[
  $Gamma = RR$、$Lambda = {0, ..., N-1}$ 时，对任意密度 $rho$，
  最优点集由比值 $r_N (y) = phi_N (y) \/ phi_(N-1)(y)$ 的水平集给出
  （$A_N (y) = r_N^(-1)(r_N (y))$）；它恰好是一条对次数 $<= 2N - 2$
  多项式精确的正求积法则的节点。取初值 $y_1$ 为 $phi_(N-1)$ 的零点时，

  $ A_N^C (Q) = A_N^C^*(Q) = A_N^F^*(Q) = A_N^F (Q) , quad abs(det bold(V)) = 1 = kappa . $

  即：一维任意分布下，贪心加权 Fekete 点几乎总是产生 Gauss 求积节点。
]

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a2_leg.png"), image("../fig_a/tmp_a2_herm25.png"))),
  caption: [条件数对比（论文 1708.01296，$M = 1.05 N$，50 次试验的分位数）。
  左：$d = 2$ Legendre（TD 与 HC 两种指标集），CFP（加权贪心）的条件数
  显著低于 AFP（无权贪心）与 MC 随机点。右：$d = 25$ Hermite，
  高维下改善幅度变小但排序不变：CFP $>=$ AFP $>=$ MC。],
)

#warn("这篇解决了什么、放弃了什么")[
  解决：AFP（经典近似 Fekete 点）需要紧支撑（平衡测度理论依赖紧性），
  加权版本 $1\/K_Lambda$ 让它在高斯这类非紧集上同样良定义；
  而且加权后贪心*有*全局最优性证明（无权版本没有）。
  放弃：确定性点集失去了随机采样的概率化误差分析；
  点数固定在 $N$ 附近，覆盖不了 $M << N$ 的重欠定区，那里仍是 $ell_1$ 的适用范围。
  两篇论文合起来恰好覆盖超定与欠定两个区段。
]

随机椭圆 PDE 实验（$kappa(x, y) = 1 + sigma sum_k cos(2 pi k x) y_k \/ (k^2 pi^2)$，
类 KL 展开形式）中，CFP 在 $d = 2$ 到 $25$ 的全部指标集上与 AFP 相当或更好、
远好于随机点，QoI 取 $u(0.5, y)$。

#figure(
  align(center, image("../fig_a/tmp_a2_ell.png", width: 78%)),
  caption: [扩散方程实验（论文 1708.01296，$d = 2$）：CFP 与 AFP 的逼近误差
  随样本量的下降曲线几乎重合且都远优于 MC。确定性点集在超定区
  是“免调参”的选择。],
)
