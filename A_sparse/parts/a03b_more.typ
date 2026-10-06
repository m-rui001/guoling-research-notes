#import "../../preamble.typ": *

= 第一、二篇补充：实验画廊与一维最优性的完整推导

== 第一篇：六组实验的完整记录

主章节给出了主线结果（Legendre 恢复、PDE）。这里补齐其余五组，
它们共同回答一个问题：网格子采样在什么情况下*不*占优。

*Chebyshev 展开*（$d = 2$ 与 $d = 10$，$M = 85$）：
Chebyshev 多项式在 Chebyshev 网格上的表现与 Legendre 在
Gauss–Legendre 网格上的表现平行，“基与网格配套”的原则在第二族上同样成立。

*单项式基的反例*：同样的采样点换单项式基 $x^alpha$，恢复完全失败。
原因在预备章的框架下很清楚：单项式不正交，Gram 矩阵严重病态，
加权 $sqrt(n lambda_n)$ 也救不了它。这组实验的作用是划界：
方法的成功绑定在“正交基 + 配套测度”上，不是采样技巧单独的功劳。

*随机 ODE*（$d u \/ d t = -k(X) u$，$k(X) = beta X$，高斯输入，$n = 30$）：
一维 Hermite 情形，Gauss 网格子采样与按原密度采样相当，
验证无界分布的 $n^(2\/3)$ 代价在实践中不严重（$n = 30$ 时
$n^(2\/3) approx 9.7$，常数因子级别）。

*Hermite 函数*（$f_1 = 2^(-0.2 x_1^2 - 0.2 x_2^2)$、
$f_2 = e^(-0.6 x_1 - 0.6 x_2)$，$d = 2$）：小 $M$ 时各方法差距不明显，
大 $M$ 时 Gauss 子采样领先；与稀疏恢复实验的结论一致：
*恢复阶段*（系数完全稀疏）对采样最敏感，逼近阶段（系数衰减但非零）
各方法都被逼近误差主导。

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a1_cheb.png"), image("../fig_a/x_a1_mono.png"))),
  caption: [左：$d = 2$ 稀疏 Chebyshev 展开的恢复概率（$M = 85$），
  网格子采样（Gaussian 列）与对照方法在两族基上行为平行。
  右：单项式基的反例：同样的采样点配单项式基，误差不随 $M$ 下降。
  方法的成功绑定在“正交基 + 配套测度”上。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a1_ode.png"), image("../fig_a/x_a1_exp1.png"))),
  caption: [左：随机 ODE（高斯输入）的逼近误差随 $M$：一维 Hermite 情形
  网格子采样可用。右：Hermite 函数 $f_1$ 的逼近误差：小 $M$ 时各方法
  差距不大，大 $M$ 时网格子采样领先。],
)

== 第二篇：一维最优性的推导链

主章节引用了“贪心 = 全局最优”与“一维退化为 Gauss 点”两个定理。
这里把一维那条的推导走完，因为它短，而且每一步都用到了
预备章的材料。

设定：$Gamma = RR$、$rho$ 任意密度、$Lambda = {0, ..., N-1}$。
要证的：加权贪心点集 $A_N^F^*$ 是一条 $n$ 点正交求积法则的节点，
因而行列式与条件数同时最优。

*第一步：最优集的结构。* 对任意不在 $phi_(N-1)$ 零点集内的
参数 $y$，定义比值

$ r_N (y) = (phi_N (y)) / (phi_(N-1)(y)) , $

最优集是水平集 $A_N (y) = r_N^(-1)(r_N (y))$：$N$ 个点 $z$ 满足
$r_N (z) = r_N (y)$。直观：$r_N$ 是“下一基函数的预测比值”，
让所有点共享同一个比值的配置使 Gram 矩阵最平衡。

*第二步：行正交性。* 加权行 $bold(psi)(z) = (phi_j (z) \/ sqrt(K_Lambda (z)))_(j=0)^(N-1)$
的两两内积分三小步算出。

*2a. 先写 Christoffel–Darboux 恒等式。* #term("Christoffel–Darboux 恒等式", [预备章提到的核函数恒等式，
把 $sum_j phi_j (x) phi_j (y)$ 化成 $phi_N phi_(N-1)$ 的差商]) 说：核函数可以收成相邻两个多项式的差商

$ sum_(j=0)^(N-1) phi_j (x) phi_j (y) = C dot (phi_N (x) phi_(N-1)(y) - phi_(N-1)(x) phi_N (y)) / (x - y) , $

其中 $C$ 是只依赖递推系数的正常数（由三项递推的 $b_j$ 组成）。
适用条件：$x != y$；$x = y$ 时右端按差商的极限理解，取正值。

*2b. 把分子拆成比值函数的差。* 按定义 $r_N (z) = phi_N (z) \/ phi_(N-1)(z)$，
通分相减：

$ phi_N (x) phi_(N-1)(y) - phi_(N-1)(x) phi_N (y)
  = phi_(N-1)(x) phi_(N-1)(y) dot (r_N (x) - r_N (y)) , $

就是把分子两项同提出公因子 $phi_(N-1)(x) phi_(N-1)(y)$，纯粹的恒等变形。

*2c. 合并。* 加权行的内积是未加权核除以两个 $sqrt(K_Lambda)$：

$ bold(psi)^T (x) bold(psi)(y)
  = (C phi_(N-1)(x) phi_(N-1)(y)) / (sqrt(K_Lambda (x)) sqrt(K_Lambda (y)))
    dot (r_N (x) - r_N (y)) / (x - y) prop (r_N (x) - r_N (y)) / (x - y) . $

前面的比例因子里 $C$、$phi_(N-1)(x) phi_(N-1)(y)$、$sqrt(K_Lambda (x) K_Lambda (y))$
全都恒正（Christoffel 函数是平方和的倒数，处处为正），所以内积的
符号与零点完全由差商 $(r_N (x) - r_N (y)) \/ (x - y)$ 决定。于是
$r_N (z_j) = r_N (z_k)$（$j != k$）当且仅当行正交。
水平集构造让所有 $r_N (z_k)$ 相等，正交性成立，行列式达到上界 1
（Hadamard 论证，主章节已给）。

*第三步：与求积法则重合。* 行正交的 $N$ 个点配上权
$1\/K_Lambda (z_k)$ 恰好构成对次数 $<= 2N - 2$ 多项式精确的
正求积法则（预备章“求积权 = Christoffel 函数”的引理反过来用）；
而 $2N - 2$ 精度的 $N$ 点法则唯一，就是 Gauss 法则。
取参数 $y$ 落在 $phi_(N-1)$ 的零点时，水平集退化成
$phi_N$ 的零点集，正是 $N$ 点 Gauss 节点。

这条推导的解释力在于：*一维情形贪心点列的极限是 Gauss 求积节点*
有一条完整的等价链撑着：“行列式最优 $=>$ 行两两正交 $=>$ 求积法则”。
高维没有这个唯一性（Gauss 求积张量积不再是唯一最优集），
所以定理只在一维给出完整等式，高维只有不等式（贪心 $<=$ 全局最优）。

== 第二篇：条件数与精度的完整记录

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a2_leg6.png"), image("../fig_a/x_a2_leg10.png"))),
  caption: [Legendre 情形的条件数（$M = 1.05 N$，50 次试验，20%/80% 分位）。
  左：$d = 6$（TD 指标集）；右：$d = 10$（HC 指标集）。
  加权 Fekete（CFP）的条件数在最坏情形仍接近 1，随机采样（MC）
  随 $N$ 增长，无权贪心（AFP）居中。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a2_her10.png"), image("../fig_a/x_a2_ell25.png"))),
  caption: [左：$d = 10$ Hermite（HC 指标集）的条件数：三种点集的差距
  比 Legendre 情形小，与第一篇“无界域 $L$ 因子较大、各方法差距缩小”
  的观察一致。右：扩散方程 $d = 25$（HC 指标集）的 QoI 逼近误差：
  高维下 CFP 与 AFP 相当、远好于 MC，验证确定性点集可以走到 $d = 25$。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a2_ell6.png"), image("../fig_a/tmp_a2_ell.png"))),
  caption: [扩散方程 $d = 2$ 与 $d = 6$ 的 QoI 逼近误差（$u(0.5, y)$，
  随机系数为类 KL 展开的余弦级数）。CFP 与 AFP 的曲线几乎重合，
  都显著低于随机采样。确定性点集在超定区是免调参的选择。],
)
