#import "../../preamble.typ": *

= LVM-GP 的实现台：门控、正则与六个算例的数字

主章节讲了“把 IB-UQ 的先验换成高斯过程”这一句。这一章把它展开成
可以照着搭的规格：编码器的两个部件各自约束什么、解码器的积分核怎么归一、
正则项从独立版改到相关版时解析式长什么样（以及一处需要读者自己判断的地方）、
GP 求导的两条路、全部超参与基线配置，最后是六个算例的完整数字。

== 编码器：一个门控，一个带空间相关的先验

$ bold(z)(bold(x); theta_E) = "diag"(m(bold(x); theta_m)) overline(bold(z))(bold(x); theta_(overline(z)))
    + "diag"(bold(1) - m(bold(x); theta_m)) bold(z)_0 , $

$ bold(z)_0 tilde op("GP")(0, K(bold(x), bold(x)')) $，
$K(bold(x), bold(x)') = sigma_K^2 exp(- norm(bold(x) - bold(x)')^2 \/ (2 l^2)) bold(I)_(d_z)$，
各分量独立、共用同一个标量核，$sigma_K$ 固定为 1。
$m$ 的每个分量取值在 $[0,1]$。对固定 $bold(x)$ 求出的边际是

$ q_E (bold(z) | bold(x)) tilde cal(N)( m(bold(x)) overline(bold(z))(bold(x)),
    "diag"( (1 - m(bold(x)))^2 ) ) , $

所以 $op("Var")(bold(z)(bold(x))) = (1 - m(bold(x)))^2$，
门控本身就是方差参数。$m -> 1$ 时潜变量退化成确定性特征，$m -> 0$ 时
完全由 GP 先验给出，而 GP 先验带空间相关，这正是标准高斯潜变量给不了的。

#warn([源码里 $m$ 的两种说法])[
  式 (6) 后面一句写“$m(bold(x); theta_m)$ 是可解释为 $bold(x)$ 远离训练数据分布的概率”，
  紧接着的下一句写“$m(bold(x)) = 1$ 意味着 $bold(x)$ 靠近某个训练数据、$bold(z)$ 成为确定性潜表示”。
  两句话对 $m$ 的取向是反的。方差式 $(1 - m)^2$ 与后一句一致，
  与图里的行为（数据密集处带窄、边界与外推处带宽）也一致，本笔记按后一句读。
]

== 解码器：归一化的高斯核积分，以及一条“噪声不求导”的假设

潜随机场先由 $L$ 层积分算子处理（$bold(z)_1$ 就是上面的 $bold(z)$，
$omega_E tilde cal(N)(0, bold(I)_M)$ 提供 GP 的实现抽样）：

$ bold(z)_i (bold(x), omega_E) = sigma( W_i bold(z)_(i-1) + b_i
    + integral_Omega k_i (bold(x), bold(x)') bold(z)_(i-1)(bold(x)', omega_E) d bold(x)' ) ) ,
    quad i = 2, dots, L-1 , $

$ bold(z)_L (bold(x), omega_E) = W_L bold(z)_(L-1)(bold(x), omega_E) + b_L , $

核取位置变换核（Positional Transformer kernel）

$ k_i (bold(x), bold(x)') = exp(- alpha_i norm(bold(x) - bold(x)')^2)
    \/ integral_Omega exp(- alpha_i norm(bold(x) - bold(y))^2) d bold(y) dot V_i , $

分母把核在 $Omega$ 上归一，于是积分项是上一层特征的加权平均而不是加权求和，
尺度不随域的大小漂移；$alpha_i > 0$ 由 $alpha_i = tan(tilde(alpha)_i)$、
$0 < tilde(alpha)_i < pi \/ 2$ 参数化。输出 $mu_u := bold(z)_L$，
解本身写成条件高斯 $u(bold(x), omega_u) = mu_u + sigma_u (bold(x)) omega_D$。

三处设计要单独记。

*第一，噪声不参与求导。* 论文假设 $sigma_u (bold(x); theta_sigma^u)$ 不对 $bold(x)$ 微分，
于是 PDE 残差 $cal(N)_bold (x)[u] = cal(N)_bold (x)[mu_u]$，
自动微分只作用在均值上。好处是残差项不必穿过随机变量；
代价是残差对噪声场的空间结构完全没有梯度，
$sigma_u$ 只能靠数据似然项去学。

*第二，$f$ 与 $b$ 的均值不是独立网络。* $mu_f = cal(N)_bold (x)[mu_u]$、
$mu_b = cal(B)_bold (x)[mu_u]$，再各加自己的 aleatoric 项。
物理一致性焊在均值表达式里，不是罚在损失里，
所以“$f$ 的预测与 $u$ 的预测是否满足方程”这个问题在这个架构下恒成立。
三份数据似然共享同一次前向。

*第三，积分算子可换。* Remark 1 给出 DeepONet 型替代：
$mu_u = "Branch"(bold(z)(bold(x)_("grid",1)), dots, bold(z)(bold(x)_p)) dot "Trunk"(bold(z)(bold(x)))$，
Branch 输入维 $p times d_z$，$d_z$ 与 $p$ 都大时会撑爆；
FNO 型的积分是 $O(N^2)$，FFT 降到 $O(N log N)$。
附录 B 用与算子型完全相同的配置跑了 DeepONet 型，
均值与不确定度都可用，没有给出谁更好的数字。

== 正则项：独立版、相关版，以及一处核对

只用数据似然时，编码器会把观测点处的 $m(bold(x)_i)$ 推到 1，
潜 GP 在整个域上塌成确定性映射。正则项按 IB-UQ 的思路加：

$ cal(L)_"reg" (Z, X) = E_(bold(x) tilde cal(U)(Omega)) D_"KL"[ q_E (bold(z) | bold(x)) ‖ e_1(bold(z))) ] , $

$e_1$ 取标准高斯 $cal(N)(0, bold(I)_(d_z))$，
经验式是 $1/N sum_i [log q_E (bold(z)_i | bold(x)_i) - log e_1(bold(z)_i)]$，
其中 $bold(z)_i tilde q_E ( dot | bold(x)_i)$、$bold(x)_i tilde cal(U)(Omega)$。
总目标是 $cal(L) = cal(L)_"data" - beta cal(L)_"reg"$，
$cal(L)_"data"$ 写成 $log q_D$ 的平均，所以要*最大化*，
算法第 3 步的更新式 $W <- W + lr dot partial cal(L) \/ partial W$ 与此一致
（源码把 $cal(L)_"data"$ 称作 negative log-likelihood，与它自己的符号相反，按式子读）。

Remark 3 的相关版把一批 $B$ 个点的边际换成联合 GP：
$e_2(bold(z)_(1:B) | bold(x)_(1:B)) tilde cal(N)(0, bold(K))$，
$bold(K) = (K(bold(x)_i, bold(x)_j))_(i,j=1..B)$，解析式（源码式 14）为

$ tilde(cal(L))_"reg" = E_(bold(x)_(1:B)) 1/(2B) [ -2 log det(bold(D)_m) - B
    + op("tr")(bold(K)^(-1) bold(D)_m bold(K) bold(D)_m)
    + (bold(m overline(z)))^op ("T") bold(K)^(-1) (bold(m overline(z))) ] , $

$bold(D)_m = "diag"(1 - m(bold(x)_1), dots, 1 - m(bold(x)_B))$。
它的用处是让 $m$ 与 GP 核的超参数（长度尺度 $L_c$）能同时被优化。

#boxed(title: "对这一项的核对")[
  两个均值不同、协方差不同的高斯之间的标准式是
  $D_"KL" (cal(N)(mu, Sigma) ‖ cal(N)(0, K))
    = 1/2 [ op("tr")(K^(-1) Sigma) - B + log det K - log det Sigma + mu^op ("T") K^(-1) mu ]$。
  代入 $Sigma = bold(D)_m^2$、$mu = bold(m overline(z))$：
  $-log det Sigma = -2 log det bold(D)_m$ 与源码一致，
  末两项的形式也一致。差别在迹项：标准式给 $op("tr")(K^(-1) bold(D)_m^2)$，
  源码写的是 $op("tr")(K^(-1) bold(D)_m K bold(D)_m)$，
  只有 $K$ 与 $bold(D)_m$ 可交换时两者才相等；
  另外 $log det K$ 被略去，若 $L_c$ 真的要跟着一起优化，这一项不是常数。
  这一条是本笔记的核对，不是论文的结论；论文自己声明了
  $B -> oo$ 的统计性质与该项的数学适定性都留待分析。
]

== GP 怎么求导（附录 A 的两条路）

*Mercer / Karhunen–Loève 展开。* $k(bold(x), bold(x)') = sum_n lambda_n phi_n (bold(x)) phi_n (bold(x)')$，
$bold(z)_0(bold(x)) = sum_n sqrt(lambda_n) xi_n phi_n (bold(x))$，
$xi_n$ 独立标准高斯；特征函数足够光滑时可逐项求导，
$cal(N)_bold (x)[bold(z)_0] = sum_n sqrt(lambda_n) xi_n cal(N)_bold (x)[phi_n]$
仍是零均值高斯过程，数值上截断到 $N$ 项。

一维情形源码给了显式展开（$L_c$ 是物理相关长度，$L_p = max {1, 2 L_c}$、$L = L_c \/ L_p$）：

$ bold(z)_(0,i)(x, omega_i) = omega_(i,1) (sqrt(pi) L \/ 2)^(1/2)
    + sum_(j=2)^N zeta_j q_j (x) omega_(i,j) , $

$ zeta_j = (sqrt(pi) L)^(1/2) exp( - (floor(j/2) pi L)^2 \/ 8 ) , quad
q_j (x) = cases( sin(floor(j/2) pi x \/ L_p), "if " j " even", cos(floor(j/2) pi x \/ L_p), "if " j " odd") , $

$zeta_1$ 那一项被单独写成常数项，展开式里的 $q_j$ 分支条件源码写的是
“if $i$ even/odd”，按上下文应为 $j$。二维不写展开，改用乘积核
$k((x_1, x_2), (x_1', x_2')) = tilde(k)(x_1, x_1') tilde(k)(x_2, x_2')$。

*函数值与导数值联合抽样。* 一维、单输出、平方指数核时

$ mat(z_0(X_u); z_0(X_f); z_0^'(X_f); z_0^''(X_f)) tilde cal(N)(0, bold(K)) , $

$bold(K)$ 的 16 个分块全由核的混合偏导解析给出
（$k_x, k_y, k_(x y), k_(x x), k_(y y), k_(x x y y)$ 等，
约定 $k_x (x, x') = partial \/ partial x' k(x, x')$）。
于是非线性算子 $cal(N)_x[z_0](x) = F(z_0, z_0^', z_0^'')$
可以在每个实现上逐点算，输出一般不再是高斯，按样本统计。

== 训练协议与两个基线的配置

#table(
  columns: (auto, auto, auto),
  align: (left, left, left),
  stroke: 0.4pt + luma(190),
  [*项*], [*取值*], [*出处 / 备注*],
  [默认网络], [3 层全连接、每隐层 64、Mish], [编码器 $m$、$overline(z)$ 与解码器同规格],
  [潜变量维 $d_z$], [20], [全部算例],
  [$beta$], [0.01；外推算例 0.1], [$beta$ 只在这两处不同],
  [优化], [Adam，初始学习率 0.001，10000 步，每 1000 步乘 0.7], [默认协议],
  [两阶段], [前 5000 步只优化均值，后 5000 步同时优化标准差并微调均值], [默认协议],
  [多孔介质], [3 层 $times$ 20；两阶段为 2000 / 8000 步], [三种方法用同一尺寸网络],
  [二维反应扩散], [3 层 $times$ 128], [其余同默认],
  [六维反演], [编解码 3 层 $times$ 128；$lambda$ 网络 2 层 $times$ 128、tanh；20000 步（10000 / 10000），学习率恒 0.001], [$sigma_u$ 初始化 0.1],
  [一维反问题的 $lambda$ 网络], [单神经元隐层 + tanh], [其余算例未用],
  [噪声参数], [$sigma_u$、$sigma_f$ 取可学标量，初始化为 2 倍经验噪声], [$epsilon tilde cal(N)(0, 0.01^2)$ 时 $sigma_u = 0.02$],
  [核参数化], [$alpha_i = tan(tilde(alpha)_i)$，$tilde(alpha)_i in (0, pi \/ 2)$], [保证正],
)

两个基线的配置要按论文写的复现，不能凭印象：
*B-PINN-HMC* 用两层、每层 50 个神经元，权重与偏置的先验是各分量独立标准高斯，
质量矩阵 $M = bold(I)$，leapfrog 步数 $L = 50$，步长 $delta_t = 0.1$，
burn-in 2000 步、总样本 5000。
*深度集成* 用与 HMC 同尺寸的网络、权重衰减 4e-6、独立训练 20 次后取集成统计量。
除多孔介质一例明确写了“HMC 与深度集成采用与 LVM-GP 相同的网络结构”，
其余算例里三者的网络尺寸并不一致（LVM-GP 默认 3 $times$ 64，基线 2 $times$ 50）。

== 六个算例的数字

*一维 Poisson（正问题）。* $lambda partial_x^2 u = f$ 在 $[-0.7, 0.7]$，$lambda = 0.01$，
精确解 $u(x) = [sin(6 x)]^3$，$f$ 由方程给出；$f$ 有 32 个等距传感器，
$u$ 在两端各一个传感器给 Dirichlet 条件；噪声两档 $epsilon tilde cal(N)(0, 0.01^2)$
与 $cal(N)(0, 0.1^2)$。低噪时三种方法都合格；高噪时深度集成对 $f$ 的不确定度失效，
LVM-GP 的均值比 HMC 更准、带更紧。
把核长度尺度 $L_c$ 也放进优化后，它从初值约 1.0 单调下降，
约 8000 步后稳定在 0.58 附近。

*一维多孔介质。* $-nu_e \/ phi partial_x^2 u + nu u \/ K = g$ 在 $[0,1]$，
$nu_e = nu = 10^(-3)$、$phi = 0.4$、$K = 10^(-3)$、$g = 1$，
解析解 $u = g K \/ nu [ 1 - cosh(r (x - H \/ 2)) \/ cosh(r H \/ 2) ]$、
$r = sqrt(nu phi \/ (nu_e K))$；$g$ 有 16 个等距传感器、$u$ 两端各一个，噪声两档。
边界层解，$r H \/ 2$ 不小，两端陡。LVM-GP 与 HMC 相当，深度集成在高噪下不可靠。

*一维非线性 Poisson（反问题），识别 $lambda$。*
$k partial_x^2 u + lambda tanh(u) = f$，$k = 0.01$，$u = [sin(6x)]^3$，$lambda$ 真值 0.7；
32 个 $f$ 传感器、2 个边界 $u$ 传感器，另有 6 个内部 $u$ 传感器专门用来定 $lambda$。

#table(
  columns: 5,
  align: (left, center, center, center, center),
  stroke: 0.4pt + luma(190),
  [噪声], [], [LVM-GP], [B-PINN-HMC], [深度集成],
  [0.01], [均值], [0.6976], [0.6967], [0.6966],
  [0.01], [标准差], [$9.816 times 10^(-3)$], [$4.225 times 10^(-3)$], [$2.493 times 10^(-4)$],
  [0.01], [$|"偏差"| \/ $ 标准差], [0.24], [0.79], [13.6],
  [0.1], [均值], [0.6965], [0.6787], [0.6959],
  [0.1], [标准差], [$6.954 times 10^(-2)$], [$4.166 times 10^(-2)$], [$3.691 times 10^(-2)$],
  [0.1], [$|"偏差"| \/ $ 标准差], [0.05], [0.51], [0.11],
)

最后一行是本笔记加的：把 $|"均值" - 0.7|$ 除以该方法的预测标准差，
得到一个粗糙的 $z$ 值。低噪档里深度集成的标准差小到 $2.5 times 10^(-4)$，
偏差是它的 13.6 倍，这是“集成在数据少时给出假精确”的数字形态；
高噪档三者的 $z$ 都落在一左右，此时集成不再离谱。
这个换算只在“预测分布是单峰且标准差是其主要宽度参数”时才有意义，
它不是覆盖率检验。

*外推。* 只用 $x <= 0$ 一侧的 4 个 $u$ 观测，配 40 个均匀 $f$ 观测，噪声标准差 0.01，
靠 PDE 约束把 $u$ 外推到 $x > 0$。论文报的参数预测均值 0.7040、标准差 $1.361 times 10^(-2)$，
与 0.7 相符。这一段源码把待估参数写成 $k$，而该方程里 $k = 0.01$ 是已知的扩散系数、
0.7 是上一例的 $lambda$，属符号混用；按数值读它估的是 $lambda$。

*二维非线性扩散反应（反问题），识别 $lambda$。*
$k (partial_(x_1)^2 + partial_(x_2)^2) u + lambda u^2 = f$ 在 $[-1,1]^2$，$k = 0.01$，
$u = sin(pi x_1) sin(pi x_2)$，$lambda$ 真值 1；
$u$ 有 100 个域内传感器、$f$ 有 484 个，每条边界 25 个 $u$ 传感器。

#table(
  columns: 5,
  align: (left, center, center, center, center),
  stroke: 0.4pt + luma(190),
  [噪声], [], [LVM-GP], [B-PINN-HMC], [深度集成],
  [0.01], [均值], [1.0003], [1.0005], [1.0047],
  [0.01], [标准差], [$4.58 times 10^(-3)$], [$5.75 times 10^(-3)$], [$4.12 times 10^(-3)$],
  [0.01], [$|"偏差"| \/ $ 标准差], [0.07], [0.46], [1.14],
  [0.1], [均值], [0.9916], [0.9781], [0.9302],
  [0.1], [标准差], [$5.70 times 10^(-3)$], [$4.98 times 10^(-2)$], [$2.60 times 10^(-2)$],
  [0.1], [$|"偏差"| \/ $ 标准差], [1.47], [0.44], [2.68],
)

高噪档里 LVM-GP 的标准差比 HMC 小一个量级，偏差与自己的标准差同量级（1.47），
而 HMC 的标准差比自己的偏差大两倍（0.44）。一个偏紧、一个偏松；
深度集成的偏差 0.0698 是其标准差的 2.68 倍，是三家里唯一明显失准的。

#figure(
  align(center, grid(columns: 2, gutter: 8pt,
    image("../fig_c/x_c6_2derr.png", width: 100%),
    image("../fig_c/x_c6_2dstd.png", width: 100%))),
  caption: [二维非线性扩散反应算例的场级结果（论文 2507.22493）：
  左为预测误差场（色标 0 到 0.149），右为预测标准差场（色标 0.103 到 0.111）。
  右图值得多看一眼：标准差场几乎是平的，
  误差场却有明显的团块结构：带宽反映了观测噪声的水平，
  没有跟上误差实际集中的位置。],
)

*六维污染源反演。* $-lambda (partial_(x_1)^2 + partial_(x_2)^2) u - f_2 = f_1$ 在 $[0,1]^2$、
零边界，$lambda = 0.02$，已知源 $f_1 = 0.1 sin(pi x_1) sin(pi x_2)$，
未知源 $f_2 = sum_(i=1)^3 k_i exp[ -0.5 norm(x - x_(c,i))^2 \/ 0.15^2 ]$，
幅度 $bold(k) = (2, -3, 0.5)$ 已知、三个中心 $x_(c,i)$ 待估；
观测为 $u$ 的 1000 个点（$epsilon_u tilde cal(N)(0, 0.1^2)$）
与 $f_1$ 的 200 个点（$epsilon_(f_1) tilde cal(N)(0, 0.01^2)$）；
参照解由 FEniCS 有限元给出（源码此处写作 "finite method method"）。
真值中心 $(0.3, 0.3)$、$(0.75, 0.75)$、$(0.2, 0.7)$。

#table(
  columns: 4,
  align: (left, left, center, center),
  stroke: 0.4pt + luma(190),
  [], [], [$x_(c,1)$], [$x_(c,2)$], [$x_(c,3)$],
  [LVM-GP], [均值], [(0.2927, 0.3022)], [(0.7433, 0.7542)], [(0.2065, 0.7569)],
  [LVM-GP], [标准差（$times 10^(-2)$）], [(2.79, 6.25)], [(2.87, 2.36)], [(5.49, 6.34)],
  [LVM-GP], [$|"偏差"| \/ $ 标准差], [(0.26, 0.04)], [(0.23, 0.18)], [(0.12, 0.90)],
  [B-PINN-HMC], [均值], [(0.3014, 0.2883)], [(0.7473, 0.7496)], [(0.2268, 0.6519)],
  [B-PINN-HMC], [标准差（$times 10^(-3)$）], [(3.08, 3.45)], [(3.51, 2.52)], [(18.97, 11.47)],
  [B-PINN-HMC], [$|"偏差"| \/ $ 标准差], [(0.45, 0.34)], [(0.73, 0.16)], [(1.41, 4.19)],
)

两个标准差列的标度不同（LVM-GP 是 $10^(-2)$ 量级，HMC 是 $10^(-3)$ 量级），
表头分别写了 $times 10^2$ 与 $times 10^3$，读的时候要换算。
差别最大的是第三个中心的 $y$ 分量：LVM-GP 给 0.7569，偏 +0.057，
自己的带宽 0.0634，偏差落在带内（0.90）；
HMC 给 0.6519，偏 −0.048，带宽 0.0115，偏差是其 4.2 倍。
这个算例的结论不是“谁更准”，而是带宽相差一个量级时同一份数据能被解释成
“在带内”与“显著偏离”两种说法。

#figure(
  align(center, grid(columns: 2, gutter: 8pt,
    image("../fig_c/x_c6_invsc.png", width: 100%),
    image("../fig_c/x_c6_extrapf.png", width: 100%))),
  caption: [左：LVM-GP 反问题的结构图（论文 2507.22493），
  $u$、$f$、$lambda$ 三个解码支路共享同一个潜随机场与同一个 $omega_E$，
  总损失写成 $cal(L)_"data" - beta cal(L)_"reg"$。
  右：外推算例的 $f$ 面板，蓝圈是 40 个均匀观测、黑线精确解、红虚线预测均值；
  图例里的 2 倍标准差带在这个纵轴尺度上几乎看不见。],
)

== 这份源码里没有的东西

#boxed(title: "缺口清单")[
  - 没有定理。稳定性、逼近阶、正则项的适定性都在结论一节被列为待做。
  - 没有成本对比。HMC 要跑 2000 次 burn-in 加 5000 个样本、集成要训 20 次，
    论文没有报任何一侧的时间，效率主张无法核对。
  - 没有校准指标。表里只有预测均值与标准差，没有覆盖率、NLL 或区间宽度，
    上一节的 $z$ 值是本笔记用现成数字做的粗算。
  - 没有时间相关问题，作者自己列为待扩展。
  - 不是算子学习设定：每个算例是单一 PDE 实例的求解，
    与 UQ-SONet、IB-UQ 的算子版本不在同一比较口径上。
  - 附录 B 的图里还留着被注释掉的 “PI-IB” 图例，
    说明曾与 IB-UQ 的物理信息版本同图对比，正稿删掉了这一列，
    因此 LVM-GP 相对 IB-UQ 的增益没有可读的数字。
]
