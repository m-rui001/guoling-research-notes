#import "../../preamble.typ": *

= 第三篇：把导数信息塞进测量矩阵（2018）

== 要解决的问题

黑箱模拟的每次查询都昂贵，但很多问题里*梯度可以廉价获得*：
伴随方程（adjoint）方法用一次额外求解就能给出目标泛函对全部参数的梯度（伴随方程是原方程的对偶形式，解它的代价与解原方程同阶，这是它相对“逐参数差分”的根本优势）。
放着不用是浪费。Peng 等（2016）已经在 Hermite + 高斯采样的组合里
把梯度信息加进 $ell_1$；Guo–Narayan–Zhou（2018）做了两件事：
给出预条件器的*一般设计框架*（覆盖有界与无界域的全部经典多项式族），
并证明加入导数测量“几乎必然”改善稀疏恢复条件。

== 问题格式

把导数行拼进测量矩阵。函数值测量与梯度测量堆叠：

$ tilde(bold(f)) = mat(bold(f); bold(f)_partial), quad tilde(bold(Phi)) = mat(bold(Phi); bold(Phi)_partial), quad bold(Phi)_partial = mat((partial bold(Phi))/(partial x_1); dots.v; (partial bold(Phi))/(partial x_d)) , $

即每个采样点贡献 $d + 1$ 行（函数值 1 行 + 各方向导数 $d$ 行），
$tilde(bold(Phi)) in RR^(N (d+1) times M)$。梯度增强基追踪：

$ arg min_(bold(c)) norm(bold(c))_1 quad "s.t." quad bold(W) tilde(bold(Phi)) bold(P) bold(c) = bold(W) tilde(bold(f)) . $

#warn("直接拼接会坏事")[
  两个新矩阵各有职责，缺一不可。$bold(W)$ 是*行预条件*（按基与采样测度加权），
  $bold(P)$ 是*列归一*（按指标加权）。数值实验显示：不预条件的 $tilde(bold(Phi))$
  的互不相干常数比原矩阵 $bold(Phi)$ *更差*：导数行的引入破坏了
  平均各向同性（$E[1\/N hat(bold(Phi))^T hat(bold(Phi))] = bold(I)$），
  采样点一偏，导数行的范数就失衡。设计 $bold(W), bold(P)$ 使 widehat 矩阵
  恢复各向同性，是本篇的核心贡献。
]

== 预条件的设计：来自导数正交性

设计能显式写出，靠的是一条正交恒等式。以 Legendre 基 + Chebyshev 采样为例：
Legendre 多项式的导数关于权 $(1 - x^2)$ 正交，由此得

$ E^c[ 2^d/rho_c (bold(z)) psi_i (bold(z)) psi_j (bold(z)) + sum_(k=1)^d (1 - z_k^2)/rho_c (bold(z)) (partial psi_i)/(partial x_k)(bold(z)) (partial psi_j)/(partial x_k)(bold(z)) ] = delta_(i j) (1 + sum_(k=1)^d c_k i_k (i_k + 1)) , $

读法：函数值行与导数行在 Chebyshev 测度 $rho_c$ 下的加权内积*依然对角*，
而且对角元显式等于 $1 + sum_k c_k i_k (i_k + 1)$（$c_k$ 是第 $k$ 个导数行的权重系数，
$i_k$ 是该指标的次数）。于是取

$ (W^0)_(n, n) = ((4\/pi^2)(1 - (z^((n))_j)^2))^(d\/4) , quad (W^j)_(n, n) = (W^0)_(n, n) (1 - (z_j^((n)))^2)^(1\/2) \/ sqrt(2) , $

$ (P)_(i, i) = (1 + sum_k c_k i_k (i_k + 1))^(-1\/2) , $

三件套合起来恰好让 $hat(bold(Phi)) = bold(W) tilde(bold(Phi)) bold(P)$ 满足平均各向同性。
Jacobi + Chebyshev、Hermite + Gaussian 两个组合有平行公式
（Hermite 最简单：导数对同一测度正交，$bold(W) = bold(I)$，
$P_(i i) = (1 + sum_k i_k)^(-1\/2)$）。一个结构性事实限制了这个框架的适用范围：
*导数仍是同族正交多项式的一元族只有 Jacobi、Laguerre、Hermite 三族*
（Hahn 1935 的经典结果），所以这套显式预条件只在这三族上闭式可写；
其他基要换成导数 #term("Gramian 矩阵", "各列两两内积组成的对称矩阵")的逆平方根（数值预处理）。

== 恢复保证

#proposition("定理（coherence 界与零空间）")[
  记 $mu_L (bold(Phi)) = "sup"_(i, bold(z)) abs(bold(Phi)_i (bold(z)))_2^2$（逐列上界型 coherence，
  Hampton–Doostan 意义）。对 Jacobi 展开配 Chebyshev 采样：

  $ mu_L (bold(Phi)) <= product_(j=1)^d 2 e (2 + sqrt(alpha_j^2 + beta_j^2)) , quad beta_L (hat(bold(Phi))) <= C product_(j=1)^d 2 e (2 + sqrt(alpha_j^2 + beta_j^2)) , $

  其中 $1 <= C <= 1 + sqrt(2)\/2$（与维数无关）。并且

  $ cal(N)(hat(bold(Phi))) subset cal(N)(bold(Phi)) , $

  $bold(Phi)$ 欠定时几乎必然是真子集。
]

读法有三层。第一，两个 coherence 常数都是“与次数无关的常数乘 $d$ 次幂”，
梯度增强没有恶化样本复杂度的阶。第二，零空间包含说明梯度增强的解空间
严格更小：能在增强矩阵下有非零零空间向量的 $ell_1$ 问题更少，
恢复严格不差（$hat(bold(Phi))$ 是 $bold(Phi)$ 的可逆后乘意义下的子矩阵）。
第三，作者诚实说明：理论上没能证明 $beta_L <= mu_L$（严格更优的 coherence），
“更好”目前由零空间包含与数值实验支撑。

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a3_mic.png"), image("../fig_a/tmp_a3_rec.png"))),
  caption: [左：三个矩阵的 MIC 随样本数的变化（$d = 2$，$n = 30$，$s = 8$）：
  预条件后的梯度增强矩阵（GraLeCh）常数最小且稳定，未预条件的拼接矩阵反而最差，
  这张图是预条件必要性的最直观证据。右：恢复率随样本数（实线）与梯度可用比例
  （虚线）上升，$d = 2$ 与 $d = 10$ 的部分梯度场景同样受益。],
)

PDE 实验（随机椭圆方程，$d = 3, 10$，梯度由伴随方程获得）：
梯度增强的均值与方差误差全面低于同信息量的标准 $ell_1$，
也低于“把导数行换成同等数量的函数值行”（standard-double），
说明梯度行携带的信息类型本身有价值。

= 第四篇：换罚函数（TL1，2018）

== 要解决的问题

$ell_1$ 有两个已知弱点：它不一定给出*最*稀疏的解（$ell_1$ 球与 $ell_0$ 球
只在正交变换下重合）；测量矩阵相干（列相关性高）时表现变差。
非凸罚更贴近 $ell_0$：$ell_p (p in (0,1))$、$ell_(1-2)$ 都有文献，
Guo–Li–Liu（2018）采用 Zhang–Xin 提出的 #term("TL1", [transformed $ell_1$]) 罚，
配 DC 算法（DCA）求解，并改进了恢复理论。

== 罚函数

$ rho_a (abs(t)) = ((a + 1) abs(t))/(a + abs(t)) , quad P_a (x) = sum_i rho_a (abs(x_i)) . $

参数 $a$ 是一个插值旋钮：$a -> 0^+$ 时 $rho_a$ 趋近逐点计数（$ell_0$，
强稀疏促进）；$a -> infinity$ 时 $rho_a -> ell_1$（凸）。
$rho_a$ 单调递增且凹，导数在零点连续有界（这是 DCA 可行的条件）。
优化问题：$min P_a (x)$ s.t. $bold(A) bold(x) = bold(b)$，或带噪版本
$min P_a (x)$ s.t. $norm(B x - b)_2 <= epsilon$。

算法（DCA-TL1）把罚写成差之差（DC）分解
$min (a+1)\/a norm(x)_1 - {(a+1)\/a norm(x)_1 - P_a (x)}$，
外层迭代算超梯度、内层 ADMM 型迭代；*自适应*版本对候选集
$a in {0.05, 0.1, 0.2, 0.3, 1}$ 各解一次、保留最稀疏的解，
解决“非稀疏目标无法先验定 $a$”的问题（稀疏函数的经验最优是 $a = 0.3$）。

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a4_l1.png"), image("../fig_a/tmp_a4_tl1.png"))),
  caption: [罚函数的等值线（论文 1805.05416）。左：$ell_1$ 的菱形等值线；
  右：TL1（$a = 1$）的等值线，坐标轴附近更“尖”，压小系数的力度更强。
  $a$ 越小等值线越贴近坐标十字（越像 $ell_0$），$a$ 越大越接近左图。],
)

== 恢复理论（比前人更干净）

#proposition("定理（无噪恢复）")[
  设 $delta_(2s) < 1/(1 + (a + 1)\/a sqrt(2))$，则 TL1 的解 $hat(x)$ 满足

  $ norm(hat(x) - x)_2 <= C_0 s^(-1\/2) P_a (x - x_s) , quad C_0 = ((6 sqrt(2) a - 2 a + 2 sqrt(2)) delta_(2s) + 2 a)/(a - ((sqrt(2) + 1) a + sqrt(2)) delta_(2s)) , $

  $x$ 恰好 $s$-稀疏时恢复精确。
]

#proposition("定理（带噪恢复）")[
  同样的 RIP 条件、$norm(e)_2 <= epsilon$ 时：

  $ norm(hat(x) - x)_2 <= C_0 s^(-1\/2) P_a (x - x_s) + C_1 epsilon , quad C_1 = (2(2a + 1) sqrt(1 + delta_(2s)))/(a - ((sqrt(2) + 1) a + sqrt(2)) delta_(2s)) . $
]

证明链沿用 Candès 的 RIP 分析模板：$P_a (hat(x)) <= P_a (x)$ 出发，
把残差 $h$ 按 $s$ 个指标分块 $h_(T_0), h_(T_1), ...$，
用三角不等式型引理 $P_a (h_(T_0^c)) <= P_a (h_(T_0)) + 2 P_a (x_(T_0^c))$
与分块和的控制 $sum_(j >= 2) norm(h_(T_j))_2 <= s^(-1\/2) P_a (h_(T_0^c))$，
最后 RIP 配 Cauchy–Schwarz 收尾。

与前人结果的对比（选一）：$ell_(1-2)$ 的恢复条件是
$delta_(2s) < 1/(1 + sqrt(2) + (2 sqrt(2))\/(sqrt(s) - 1))$（依赖 $s$），
TL1 在 $a > (sqrt(s) - 1)\/2$ 时条件更宽松、误差常数更简洁；
Zhang–Xin 原始定理的条件含 $delta_R$ 与 $delta_(R + abs(T))$ 两项，
本篇化简为单一 $delta_(2s)$。

== 在 gPC 上的落点

#proposition("定理（稀疏 Legendre gPC 恢复）")[
  测量矩阵由 Legendre 多项式 + 均匀测度采样构成。若

  $ M >= C delta^(-2) 3^P s log^3 (2 s) log (N) , quad delta < 1/(1 + (a+1)\/a sqrt(2)) , $

  则以概率超过 $1 - N^(-gamma log^3 (2s))$：无噪
  $norm(hat(x) - x)_2 <= C_0 s^(-1\/2) P_a (x - x_s)$；带噪再加 $C_1 epsilon$。
]

证明链：Legendre 加权 $(1 - x^2)^(1\/4)$ 构成一致有界系统 → Rauhut–Ward RIP
→ 代入 TL1 恢复定理。注意采样侧完全沿用第一篇的框架：
*这篇换的是罚函数，不是采样*，两篇是正交的改进，可以叠加。

#figure(
  align(center, grid(columns: 2, gutter: 10pt, image("../fig_a/tmp_a4_rec.png"), image("../fig_a/tmp_a4_kdv.png"))),
  caption: [左：$d = 2$ 稀疏恢复概率随 $M$ 与 $s$ 的变化，排序为
  TL1 $>$ $ell_(1-2)$ $>$ $ell_1$。右：随机 KdV 方程（随机强迫 KL 展开）
  的 gPC 系数幅值图：TL1 恢复出的系数最稀疏，
  非零系数的位置与真值一致而幅值更集中。],
)
