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
  按定义把离散内积摆出来：$psi^i_(j, n)(z^i_p) psi^i_(k, n)(z^i_p)
  = n lambda^i_n (z^i_p) phi^i_j (z^i_p) phi^i_k (z^i_p)$，于是

  $ 1/n sum_(p=1)^n psi^i_(j, n)(z^i_p) psi^i_(k, n)(z^i_p)
    = 1/n sum_(p=1)^n n lambda^i_n (z^i_p) phi^i_j (z^i_p) phi^i_k (z^i_p) . $

  下面把这行一步步变成 $delta_(j, k)$。

  *第一步（消去因子）。* 求和号外的 $1\/n$ 与每个被加项里的 $n$ 相消：

  $ 1/n sum_(p=1)^n n lambda^i_n (z^i_p) phi^i_j (z^i_p) phi^i_k (z^i_p)
    = sum_(p=1)^n lambda^i_n (z^i_p) phi^i_j (z^i_p) phi^i_k (z^i_p) . $

  *第二步（Christoffel 函数换成求积权）。* 预备章的引理“Gauss 权与
  Christoffel 函数”说：测度归一时，$n$ 点 Gauss 节点上的求积权等于
  用 $n$ 个基函数定义的 Christoffel 函数。这里的 $lambda^i_n$ 恰好用
  $phi^i_0, dots, phi^i_(n-1)$ 这 $n$ 个基定义，节点 $z^i_p$ 恰是
  $n$ 点 Gauss 节点，两个条件都对上，所以
  $lambda^i_n (z^i_p) = w^i_p$：

  $ sum_(p=1)^n lambda^i_n (z^i_p) phi^i_j (z^i_p) phi^i_k (z^i_p)
    = sum_(p=1)^n w^i_p phi^i_j (z^i_p) phi^i_k (z^i_p) . $

  这一步的适用条件要核对两点：采样点必须是 Gauss 节点（引理只对节点成立），
  测度必须归一（否则差一个总质量因子）。

  *第三步（求积换成积分）。* $phi^i_j phi^i_k$ 的次数不超过
  $(n - 1) + (n - 1) = 2n - 2$，落在 $n$ 点 Gauss 法则的精确范围
  $<= 2n - 1$ 内，所以和式等于连续积分：

  $ sum_(p=1)^n w^i_p phi^i_j (z^i_p) phi^i_k (z^i_p)
    = integral phi^i_j phi^i_k rho^i . $

  *第四步（用基的正交性）。* $phi^i_k$ 关于 $rho^i$ 正交归一，右端就是
  $delta_(j, k)$。四步连起来即所要的等式。

  多元版本：$psi_(bold(k), bold(n))(z) = product_i psi^i_(k_i, n_i)(z^i)$
  在张量积离散测度 $nu_(bold(n)) = times_i nu^i_n$ 下的内积按维分解成
  一元内积的乘积，每个因子是上面的 $delta_(j, k)$，只有全部维度
  指标都相同时乘积才为 1，故 $psi_(bold(k), bold(n))$ 在 $nu_(bold(n))$
  下同样正交归一。$square$
]

#insight("用日常语言读一遍")[
  这几步合起来只做一件事：把“连续测度下正交”这个资格复制到
  “网格上等权抽样”的离散内积下。翻译官是 Christoffel 函数：
  它在 Gauss 节点上的取值恰好等于求积权，于是“多项式乘
  $sqrt(n lambda_n)$”在网格上的平均等于“多项式乘求积权的平方根”
  在连续测度下的平均，正交性原样过河。后面一切 RIP 论证都站在这条引理上。
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
  是正交矩阵（第二步的引理张量化后的矩阵形式）。逐项估计每个元素的
  大小。先代入元素定义：多元基取乘积，所以矩阵元素按维分解为

  $ abs(A_(p, q)) = product_(i=1)^d abs(sqrt(w^i_(l_i (p))) phi^i_(k_i (q))(z^i_(l_i (p)))) , $

  再在每一维上乘进 $sqrt(n_i)$，并逐维套用第一步的恒等式
  $lambda^i_n (z^i_p) = w^i_p$（Gauss 节点上 Christoffel 函数等于求积权）：

  $ (product_(i=1)^d sqrt(n_i)) abs(A_(p, q))
    = product_(i=1)^d sqrt(n_i w^i_(l_i (p))) abs(phi^i_(k_i (q))(z^i_(l_i (p))))
    = product_(i=1)^d abs(psi^i_(k_i (q), n_i)(z^i_(l_i (p)))) , $

  最后一个等号就是 $psi$ 的定义 $psi = sqrt(n lambda) phi$（求积权换成
  Christoffel 函数后与定义式逐字相同）。于是第一步的逐点界
  $sup abs(psi)^2 <= L_i (n_i)$ 对每个因子适用：

  $ product_(i=1)^d abs(psi_(dot, n_i)) <= product_(i=1)^d sqrt(L_i (n_i)) = sqrt(product_(i=1)^d L_i (n_i)) = sqrt(L(bold(n))) . $

  矩阵每个元素都被 $sqrt(L(bold(n)) \/ product_i sqrt(n_i))$ 压住，再引用
  “正交矩阵无放回行子采样仍保持 RIP”的结果完成。$square$
]

#figure(
  align(center, image("../fig_a/tmp_a1_bounds.png", width: 90%)),
  caption: [上界 $L_i$ 的形状（论文 1602.00995）。左：Jacobi 情形的界随参数变化，
  与 $n$ 无关；右：Hermite/Laguerre 情形的 $n^(2\/3)$ 界。有界分布
  （Legendre/Chebyshev）是这条路的最爱：$L$ 是常数，样本量只比最优多
  $C^d$ 的维数代价（数值上接近 $2^d s log^3 s log N$）。],
)

#warn([无界分布的 $n^(2\/3)$ 是本质的])[
  Hermite（高斯输入）情形 $L_i tilde n^(2\/3)$ 在子采样分析策略下
  已经是这种策略的下确界（sharp），不是可以靠更细的估计抹掉的技术瑕疵。
  直觉：无界域上 Gauss 节点越跑越远，
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
  记加权行 $bold(psi)(y_j) = (psi_alpha (y_j) \/ sqrt(K_Lambda (y_j)))_(alpha in Lambda)$，
  即第 $j$ 个采样点上全体基函数取值除以同一个数 $sqrt(K_Lambda (y_j))$。

  *第一步（每行长度恒为 1）。* 算行向量的平方长度，把定义
  $K_Lambda (y) = sum_(alpha in Lambda) psi_alpha^2 (y)$ 用在分子上：

  $ norm(bold(psi)(y_j))_2^2
    = sum_(alpha in Lambda) ((psi_alpha (y_j))^2)/(K_Lambda (y_j))
    = (K_Lambda (y_j))/(K_Lambda (y_j)) = 1 . $

  合法性：每个分量除的是同一个标量，可以提到求和号外；分子恰好就是
  分母的定义式，所以不涉及任何不等式，是恒等式。

  *第二步（行列式的上界）。* Hadamard 不等式说：方阵的行列式绝对值
  不超过各行长度的乘积，$abs(det bold(V)) <= product_j norm(bold(psi)(y_j))_2$。
  第一步给出行长度全是 1，于是

  $ abs(det bold(V)(A_N, Q)) <= product_j 1 = 1 . $

  Hadamard 不等式对所有方阵成立，这里不需要矩阵有任何特殊结构。

  *第三步（等号条件翻译成条件数）。* Hadamard 不等式取等号当且仅当
  各行两两正交。行长度为 1 且两两正交意味着 $bold(V) bold(V)^T = bold(I)$，
  即 $bold(V)$ 是正交矩阵；正交矩阵的全部奇异值都是 1，于是
  $kappa(bold(V)) = 1$（条件数的理论下界，第 0 章的条件数定义）。
  这就证明了 $abs(det bold(V)) = 1$ 与 $kappa = 1$ 同时成立、互为等价。

  *第四步（贪心增量）。* 设已选 $n$ 个点，候选新点 $y$ 对应加权行
  $bold(psi)(y)$。对已选行的转置做 QR 分解，新行分解成
  “已选行空间内的投影 + 垂直分量”，垂直分量的长度平方记
  $1 - norm(bold(b))_2^2$（行总长为 1，投影吃掉 $norm(bold(b))_2^2$）。
  行列式按 QR 的上三角因子累乘，得到乘法更新

  $ abs(det(V_(n+1) V_(n+1)^T))
    = abs(det(V_n V_n^T)) dot (1 - norm(bold(b))_2^2) <= abs(det(V_n V_n^T)) , $

  等号当且仅当 $norm(bold(b))_2 = 0$，即新点的加权行与已选行全部正交。
  这条更新律说明：只要每一步都存在与已选集合正交的候选点
  （最优集中的点按构造可以逐个这样达成），沿这条分支贪心的
  行列式增量恒取最大值，贪心解就是全局最优解。$square$
]

#insight("用日常语言读一遍")[
  这段证明把两个优化目标接在一起：Hadamard 不等式把“行列式最大”
  翻译成“行两两正交”，而行正交恰好就是矩阵正交、条件数为 1。
  于是“选行列式最大的点”与“选条件数最好的点”在加权之后
  变成同一件事，贪心算法只要每步保住正交性，就不会掉出全局最优。
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
