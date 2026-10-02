#import "@preview/cetz:0.3.4": canvas, draw

#import "../../preamble.typ": *

= 预备知识（上）：正交多项式与 Gauss 求积

六篇论文的每一处推导都踩在正交多项式的三条性质上。这里把这三条备齐，
已经熟悉的读者可以跳到本章最后一节（两种约定警告），那里有一个容易踩的坑。

== 三项递推

区间 $Gamma$ 上的测度 $rho$ 给出一列正交多项式 $phi_0, phi_1, dots$
（$phi_n$ 恰好 $n$ 次，$integral phi_n phi_m rho = 0$ 对 $n != m$）。
任何这样的序列满足*三项递推*：

$ x phi_n (x) = b_(n+1) phi_(n+1)(x) + a_n phi_n (x) + b_n phi_(n-1)(x) , $

系数 $a_n, b_n$ 由测度决定。这一条是构基算法（第五篇的 Hankel–Cholesky
路线）的出口：矩序列决定 $a_n, b_n$，递推决定整列基。常用四族的参数：

#ttable(
  ([族], [权 $rho(x)$], [支撑], [对应分布]),
  (
    ([Legendre], [1], [$[-1, 1]$], [均匀]),
    ([Chebyshev], [$(1 - x^2)^(-1\/2)$], [$[-1, 1]$], [（无直接对应，平衡测度）]),
    ([Hermite], [$e^(-x^2)$], [$RR$], [高斯]),
    ([Laguerre], [$e^(-x) x^alpha$], [$[0, infinity)$], [指数]),
  ),
  widths: (auto, auto, auto, auto),
)

Jacobi 族把 Legendre 与 Chebyshev 统一为 $rho prop (1 - x)^gamma (1 + x)^delta$
的两个端点，第一篇的“任意 Beta 权”就是任意 Jacobi 参数。

== Gauss 求积与 2n-1 精度

#proposition("定理（Gauss 求积精度）")[
  设 $z_1, ..., z_n$ 是 $phi_n$ 的零点，$w_k$ 是对应求积权
  （使 $n$ 点法则对次数 $<= n - 1$ 的多项式精确），则该法则对
  一切次数 $<= 2n - 1$ 的多项式精确：

  $ integral_(Gamma) p(x) rho(x) d x = sum_(k=1)^n w_k p(z_k) quad "对" quad deg p <= 2n - 1 . $
]

#proof[
  取次数 $<= 2n - 1$ 的 $p$，用 $phi_n$ 做带余除法：$p = q phi_n + r$，
  其中 $deg q <= n - 1$、$deg r <= n - 1$。
  在节点上求和：$sum_k w_k p(z_k) = sum_k w_k q(z_k) underbrace(phi_n (z_k), 0) + sum_k w_k r(z_k) = sum_k w_k r(z_k)$。
  第一个求积（$q phi_n$，次数 $<= 2n - 1$）超出了法则的保证精度，
  但在每个节点上 $phi_n (z_k) = 0$ 使它逐项为零；
  第二个求积（$r$，次数 $<= n - 1$）在保证精度内，等于
  $integral r rho$。而正交性给出 $integral q phi_n rho = 0$，于是
  $integral r rho = integral (q phi_n + r) rho = integral p rho$。$square$
]

这就是第一篇“离散正交系统”引理里求积精确性的出处：
Christoffel 权重搬进取权和之后，恰好把正交性从连续测度搬到离散网格上。

== 求积权就是 Christoffel 函数

#lemma("引理（Gauss 权与 Christoffel 函数）")[
  $n$ 点 Gauss 求积在节点 $z_k$ 处的权与 $n - 1$ 阶 Christoffel 函数满足

  $ w_k = lambda_(n-1)(z_k) dot integral rho , quad lambda_(n-1)(x) = 1 / sum_(j=0)^(n-1) phi_j (x)^2 . $
]

#proof[
  取“指示插值”问题：对单个节点 $z_k$，构造次数 $<= n - 1$ 的多项式
  $ell_k$ 使 $ell_k (z_k) = 1$、$ell_k (z_j) = 0$（$j != k$）。
  对 $ell_k$ 用求积精确性：$w_k = integral ell_k rho$。
  而 $ell_k$ 可以由正交投影写出：$ell_k (x) = lambda_(n-1)(z_k) sum_j phi_j (z_k) phi_j (x)$
  （这是 #term("Christoffel–Darboux 恒等式", [把核函数 $sum_j phi_j (x) phi_j (y)$
  写成相邻次多项式的差商的标准恒等式]) 的直接推论），代入 $x = z_k$
  即得 $w_k = lambda_(n-1)(z_k) sum_j phi_j (z_k)^2 = lambda_(n-1)(z_k)$。
  测度归一（$integral rho = 1$）时两者相等。$square$
]

第一篇的预条件 $sqrt(n lambda_n^i)$ 里的 $lambda$ 因此有双重身份：
它既是离散网格上的求积权，又是连续测度下的 Christoffel 函数。
“求积权的平方根乘多项式在离散测度下正交”这条引理，
本质上是把连续正交性沿求积法则取样。


== 从递推系数到节点与权：Golub–Welsch 算法

三项递推系数 $a_n, b_n$ 不仅生成基函数，还直接给出 Gauss 节点与权：
把递推写成矩阵形式 $x bold(phi)(x) = bold(J) bold(phi)(x) + b_(n+1) phi_(n+1)(bold(e))_(n)$，
其中 $bold(J)$ 是对称三对角矩阵（对角元 $a_0, ..., a_(n-1)$，
次对角元 $sqrt(b_1 b_2 dots.c b_(n-1))$，即 *Jacobi 矩阵*）。
Gauss 节点是 $bold(J)$ 的特征值，求积权是特征向量的首分量平方
乘测度总质量。这就是 #term("Golub–Welsch 算法", "把求积节点与权的计算化归为对称三对角矩阵特征问题的标准算法")：一次对称三对角特征分解
（$O(n^2)$）同时给出节点与权，数值稳定性远好于显式构造正交多项式
再求根。

第五篇的构基走的就是这条路：Hankel 矩阵 Cholesky 分解给出
$a_n, b_n$（主章节的公式），Golub–Welsch 给出节点与权
（若需要配点），整条链没有一处需要显式写出多项式表达式。

== 一个可以手算的例子

Legendre（$rho = 1$ on $[-1, 1]$，$integral rho = 1$）：
$phi_0 = 1$。由正交性 $integral phi_1 d x = 0$ 与 $integral x phi_1^2 d x = 1$
解出 $phi_1 = sqrt(3) x$。再设 $phi_2 = alpha x^2 + beta x + gamma$，
正交于 1 与 $phi_1$ 给出 $gamma = -alpha\/3$、$beta = 0$，
归一化 $integral phi_2^2 = 1$ 给出 $alpha = sqrt(5)\/2$，
即 $phi_2 = sqrt(5) (3 x^2 - 1)\/2 \/ sqrt(2) ...$（具体常数取决于
归一化约定，形状不变）。递推系数 $a_0 = 0$（对称测度）、
$b_1 = 1\/3$、$b_2 = 4\/45$ 的倒数结构已经能看出：
Jacobi 矩阵的三对角元里编码了测度的全部信息。

这个例子还有一个用途：验证三项递推的次数关系。$x phi_1 = sqrt(3) x^2$
展开到 $phi_2$ 与 $phi_0$ 上，$phi_2$ 的系数是
$sqrt(3) \/ (sqrt(5)\/2 \/ ... )$。读者可以用上面的结果亲手核对，
两分钟内能走完，走完之后“递推系数决定一切”就不再是口号。

== 为什么 Hermite 的 $L$ 因子降不下来

第一篇的三类界里，Jacobi 是常数、Hermite/Laguerre 是 $C n^(2\/3)$。
差别来源可以在 Christoffel 函数的形状里看到。有界域上
（Legendre），$sum_j phi_j^2(x)$ 在端点最大、内部平缓，
Christoffel 函数 $K$ 在端点趋于零但*处处有正的下界型行为*，
加权 $sqrt(n lambda_n) phi_j$ 全域有界，即一致有界系统。

无界域上（Hermite），Gauss 节点随 $n$ 外跑，最外节点上的
$lambda_n$ 指数式地小，而多项式值只有代数式地大，
两者之比 $sqrt(n lambda_n) abs(phi_j)$ 的上界随 $n$ 增长
（幂律 $n^(2\/3)$，Mhaskar–Rahmanov–Saff 尺度 $a_n tilde n^(1\/2)$
与 Levin–Lubinsky 修正 $n^(-2\/3)$ 的差就是那个幂）。
这不是分析技术不够好：换任何“只在有限网格上采样”的策略，
最外层节点都会把一致界顶到幂律上。要绕开它只能换策略
（按原密度采样，第一篇的 $d = 10$ Hermite 数值正是这样更好），
或者接受幂律。这是“结构化网格”路线的边界，也是第六篇
induced 采样（网格无关）的动机之一。

== 一个容易踩的坑：两种约定

#definition("Christoffel 函数的两种写法")[
  文献里“Christoffel 函数”有两个互为倒数的版本：

  $ K(x) = N / sum_j Phi_j^2(x) quad quad "与" quad quad kappa(x) = 1/N sum_j Phi_j^2(x) . $

  $K$ 在多项式能量小的地方*大*（采样理论偏好这种点）；
  $kappa$ 在能量大的地方*大*（预条件用它的倒数压低能量）。
  本笔记第一、二、五篇用 $K$，第六篇用 $kappa$（与原文一致）。
  读任何一篇采样论文时，先查它的定义式再读采样方向的陈述，
  两个约定差一个倒数，方向完全相反。
]

= 预备知识（下）：稀疏恢复的几何图像

== 从 $ell_0$ 到 $ell_1$

恢复稀疏向量的“诚实”目标是数非零分量：$min norm(c)_0$ s.t. $Psi c = f$。
这是组合问题（NP 难）。凸松弛把计数换成绝对值和：

$ min norm(c)_1 quad "s.t." quad Psi c = f . $

为什么 $ell_1$ 会给出稀疏解？二维的几何图像：

#figure(
  align(center, canvas({
    import draw: *
    // ell1 菱形
    let d = 1.5
    draw.line((d,0), (0,d), (-d,0), (0,-d), close: true, fill: rgb("#dbe7f7"), stroke: 1pt + rgb("#1a4f9a"))
    // ell2 圆
    draw.circle((0,0), radius: d*0.85, fill: none, stroke: 0.8pt + rgb("#aaaaaa"))
    // 约束线 Psi c = f
    draw.line((-2.2, 1.9), (2.2, -0.35), stroke: 1.2pt + rgb("#b8860b"))
    // 解点：菱形顶点
    draw.circle((0,-d), radius: 0.06, fill: rgb("#a33"))
    draw.content((0.25, -d - 0.28), text(size: 8pt, fill: rgb("#a33"))[基追踪解：落在坐标轴顶点])
    draw.content((-2.0, 2.15), text(size: 8pt, fill: rgb("#b8860b"))[$Psi c = f$ 的解集])
    draw.content((1.15, 1.35), text(size: 8pt, fill: luma(90))[$ell_2$ 球（岭回归）])
    draw.content((-1.7, -1.15), text(size: 8pt, fill: rgb("#1a4f9a"))[$ell_1$ 球（基追踪）])
  })),
  caption: [$ell_1$ 促进稀疏的几何原因（示意）。扩大 $ell_1$ 球（菱形）时
  最先碰到约束线的是菱形的*顶点*，即只有一两个非零分量的点；
  $ell_2$ 球（圆）最先碰到约束线的地方一般在象限内部，得到的解是稠密的。
  高维时菱形的“顶点”占绝对多数，这个效应更强。],
)

== 两条判据的分工（与第二章衔接）

第二章给出的 RIP 与 MIP 是同一件事的两种读法。RIP 要求矩阵对*所有*
稀疏向量近似等距，是整体性质，随机采样下以高概率成立，
样本复杂度里有 $log$ 因子；MIP 只看*最坏的一对列*，条件强
（$mu < 1\/(2s - 1)$），但可以用逐点估计直接验证，不需要概率论证。
六篇论文的分工：走概率路线的第一、六篇用 RIP，
走确定性路线的第三篇（测量矩阵是手工拼的）用 MIP，
第四篇（非凸罚）沿用 RIP 但常数更好。

还有一个实用事实值得放在这里：RIP 检验本身是 NP 难的，
实践中没人去验一个给定矩阵是否满足 RIP；
理论的用法是“设计采样方式，让高概率成立”，而不是“拿到矩阵再检查”。
这也是为什么六篇论文的全部精力都花在*采样测度*上：
测度决定矩阵的系综，系综决定 RIP 概率。

== 最小二乘一侧的稳定性问题

第二、五篇走超定最小二乘路线，那里的病态来源与欠定情形不同：
即使 $M > N$（方程够多），Vandermonde 型矩阵的条件数仍可能爆炸。
一维的极端例子：等距节点上的单项式基，条件数随 $N$ 指数增长
（Lebesgue 常数发散），数据里的浮点噪声被放大到完全淹没信号。
第二篇开头的“等距 Newton–Cotes 点即使数据无穷精确也不稳定”说的就是它。

两条修复路线在谱系里都出现过：换*基*（Chebyshev 权、加权基
$phi \/ sqrt(K)$，把能量摊平）与换*点*（Fekete/Leja 式贪心选点）。
第二篇证明加权之后两条路线合流：贪心选出的点让加权矩阵的条件数
达到理论下界 1。
