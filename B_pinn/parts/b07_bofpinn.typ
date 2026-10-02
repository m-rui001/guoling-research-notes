#import "../../preamble.typ": *

= 第五篇：三难合一（BO-fPINN，2023）

== 问题：三个难点同时出现

现在把前面所有难点叠在同一个方程上。时间依赖的随机分数阶方程（2 维 Riesz 形式）：

$ (partial u(x, t; omega))/(partial t) = D_(|x_1|)^alpha u + D_(|x_2|)^beta u + K(x; omega) f(u) + g(x, t; omega) . $

三个成分各自的困难：随机系数与随机强迫（维数灾难，第 1 章）；
长时间积分（模态基要跟着演化，第 6 章）；
分数阶 Riesz 导数 $D_(|x_1|)^alpha$（非局部 + 奇异，不能自动微分，第 10 章）。

现有工具各管一段。经典 DO/BO 谱方法（Zhao 等 2021 的分数阶版本）
需要 QMC 或多项式混沌先算一段“起步段”来做 KL 分解，在特征值交叉或协方差奇异时不稳；
fPINN 的 GL 格式管分数阶但不管随机；
MC-PINN 管分数阶与随机但不用模态。Ma, Li, Zeng, Guo, Karniadakis（JCP 2023）
的 BO-fPINN 是一个直接拼装：*BO 弱形式（隐式约束）+ GL 离散（分数阶）*。

== Riesz 导数与 GL 离散

一维 Riesz 分数阶导数由左右 Riemann–Liouville 导数（阶同为 $alpha$）组合而成：

$ D_(|x_1|)^alpha u = -1/(2 cos(pi alpha\/2)) (D_(a, x_1)^"RL" u + D_(x_1, b)^"RL" u) , $

$ D_(a, x)^"RL" u(x) = 1/(Gamma(n - alpha)) [ d^n/(d x^n) integral_a^x (x - s)^(n - alpha - 1) u(s) d s ] , $

其中 $n$ 是满足 $n - 1 <= alpha < n$ 的整数。负号与 $1 \/ cos(pi alpha \/ 2)$
让算子在傅里叶意义下对应符号正确的耗散（$alpha = 2$ 时退化回 $-partial^2$，
$cos pi = -1$，读者可以自己验算）。

数值离散用二阶移位 Grünwald–Letnikov 格式。均匀网格 $Delta x$ 上：

$ D_(|x|)^alpha u(x) approx -1/(2 cos(pi alpha\/2)) [ alpha/2 delta_(Delta x, 1)^alpha u(x^j) + (1 - alpha/2) delta_(Delta x, 0)^alpha u(x^j) ] , $

$ delta_(Delta x, p)^alpha u(x^j) := 1/(Delta x)^alpha [ sum_(k=0)^j w_k^alpha u(x^j - (k - p) Delta x) + sum_(k=0)^(N - j) w_k^alpha u(x^j + (k - p) Delta x) ] , $

系数递推 $w_0^alpha = 1$，$w_k^alpha = (1 - (alpha + 1)\/k) w_(k-1)^alpha$。
两个移位参数 $p = 0, 1$ 的凸组合把一阶 GL 格式升到二阶精度。
关键点：离散算子 $cal(L)_(Delta x)^alpha$ 作用在网络输出上只需要*函数值*
（大量的历史网格点值），不需要自动微分；自动微分只负责整数阶项与时间导数。
这与第 10 章的 MC 路线是同一障碍的两条绕行路：MC 靠采样，GL 靠网格历史。

== 模态展开与损失（继承第 6 章）

展开形式与 BO 条件与第 6 章完全相同：

$ u_N (x, t; omega) = overline(u)(x, t) + sum_(i=1)^N u_i (x, t) Y_i (t; omega), quad ⟨u_i, u_j⟩ = lambda_i (t) delta_(i j), quad E[Y_i Y_j](t) = delta_(i j), quad E[Y_i] = 0 . $

弱形式的三条投影把 $cal(N)_x$ 换成 $cal(N)_x^(alpha, beta) := D_(|x_1|)^alpha + D_(|x_2|)^beta + K f(u)$：

$ E[partial_t u] = partial_t overline(u) = E[cal(N)_x^(alpha, beta)[u] + g] , quad ⟨partial_t u, u_i⟩ = ⟨cal(N)_x^(alpha, beta)[u] + g, u_i⟩ , quad E[partial_t u Y_i] = E[(cal(N)_x^(alpha, beta)[u] + g) Y_i] . $

BO 演化方程在分数阶情形的版本（附录定理 A.1，结构同第 6 章定理）：
定义 $S_(i j) = ⟨u_i, partial_t u_j⟩$，$M_(i j) = E[Y_i partial_t Y_j]$，
则

$ lambda_i (d Y_i)/(d t) = -sum_(j=1)^N S_(i j) Y_j + ⟨cal(N)_x^(alpha, beta)[u] + g - E[cal(N)_x^(alpha, beta)[u] + g], u_i⟩ , $

$ (partial u_i)/(partial t) = -sum_(j=1)^N M_(i j) u_j + E[(cal(N)_x^(alpha, beta)[u] + g) Y_i] , $

无特征值交叉时 $S, M$ 有与第 6 章相同形式的闭式
（$M_(i j) = (G_(i j) + G_(j i)) \/ (-lambda_i + lambda_j)$）。
经典谱方法（QMC-BO 混合算法）必须依赖这个闭式：
先跑 QMC 到某个时刻 $t_s$，做 KL 分解初始化，再切换到 BO 推进。
起步段贵，交叉处崩，这是 BO-fPINN 用隐式损失消除的两个依赖。

总损失把第 6 章的弱形式逐项搬过来，只是 $cal(N)_x$ 全部替换为 GL 离散算子
$cal(L)_(Delta x)^alpha$：

$ cal(L)(theta) = lambda_w "MSE"_w + lambda_"IC" "MSE"_"IC" + lambda_"BC" "MSE"_"BC" + lambda_"BO" "MSE"_"BO" + lambda_g "MSE"_g . $

损失权重不再手工固定：按 Wang–Perdikaris 的梯度病理方法动态学习；
振幅网络 $A_(n n)$ 拆成 $N$ 个独立小网络，因为各阶振幅跨尺度振荡。

== 两个实验亮点

*特征值交叉基准*：一维随机分数阶反应扩散方程，人工构造解使得振幅
$a_1(t), a_2(t)$ 在时域内多次交叉。经典 BO 方法在每个交叉点都需要特殊处理，
BO-fPINN 用隐式约束完全不受影响，全部分量相对误差不超过 1.3%。

#figure(
  align(center, image("../fig_b/tmp_crossing.png", width: 78%)),
  caption: [特征值交叉基准中的振幅演化（论文 2303.10913）。
  $a_1(t)$ 与 $a_2(t)$ 多次相撞又分开；经典 BO 的闭式解在每次相撞时分母为零，
  BO-fPINN 的隐式约束损失在整个时域上光滑。],
)

*分数阶迁移学习*：分数阶方程随阶数 $alpha$ 平滑变化，
因此低阶方程训练好的网络是高阶方程的优良初始化。
做法：在 $alpha = 1.8$ 上完整训练，把网络权重 $W^alpha$ 迁移给目标方程
$tilde(alpha) = 1.2$ 或 $1.5$，只用 L-BFGS-B 微调。
结果：完整训练（40 万步 Adam + L-BFGS-B）的相对时间是迁移微调的 58.2 倍，
精度相当。阶数平滑变化带来的这条捷径是分数阶问题特有的。

#figure(
  align(center, image("../fig_b/tmp_transfer.png", width: 82%)),
  caption: [分数阶迁移学习（论文 2303.10913）。源方程 $alpha = 1.8$ 的权重
  作为目标方程（$tilde(alpha) = 1.2, 1.5$）的初始化，只做局部优化微调。
  右侧曲线对比完整训练与迁移微调的误差收敛。],
)

其余实验：19 维随机输入的正向问题（相关长度短，需要 19 个 KL 模态），
6 至 8 个 BO 模态收敛，长时间积分 $T = 5$ 用时域分段训练，
gPC 对照组误差最大（基不随时间演化的问题在长时段暴露）；
反问题同时识别反应系数、噪声强度与分数阶阶数
（真值 $0.5, 0.3, 1.5$，识别值 $0.508, 0.291, 1.486$）；
二维双 Riesz 算子的反应扩散方程均值与方差全部吻合。

== 本篇的位置

BO-fPINN 本身没有引入全新的表示方法，它的价值是*组装*：
证明模态弱形式（第 6 章）与分数阶离散（第 10 章的另一条路）可以共存于一个损失函数，
并且互相不干扰。到此为止，这条研究线能解的方程越来越复杂，
但表示随机性的工具始终是“模态展开”。剩下的那个根本限制（模态个数随有效维数增长、空间维数高时代价爆炸）还没有被正面攻击。
下一篇处理它。
