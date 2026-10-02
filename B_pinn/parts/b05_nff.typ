#import "../../preamble.typ": *

= 第三篇：可逆变换与严格密度（NFF，2022）

== 新困难：模态展开的两处不严格

到第 6 章为止，随机性一直用“模态 + 正交基”表示。Guo, Wu, Zhou（JSC 2022）指出这条路的两个代价，
然后换了一条完全不同的路。

*代价一：基函数个数仍然指数增长。* aPC 的项数是 $(r + M)! \/ (r! M!)$，
模态方法的展开阶数受有效维数支配。相关长度短、非高斯的随机场会把阶数推得很高。

*代价二：模态方法只给统计量，不给密度。* 从模态重构可以算均值和方差，
但要回答“解在某点的完整概率分布长什么样”“给定一处测量后其他位置的条件分布是什么”，
模态展开没有直接答案。GAN 类方法（把随机场直接交给对抗训练）可以采样，
但训练不稳定，而且*完全不能评估密度*：无法算出一个样本出现的概率。

NFF（#term("Normalizing Field Flows", "正则化场流")）用的工具是 normalizing flow：
一条*可逆*变换链，把一个简单的参考分布精确地搬成目标分布。
它的两样东西都是精确的：采样（从参考分布采样再正向变换）和密度（变量替换公式逐层算）。

== 基础：变量替换公式

#definition("什么叫“normalizing flow”")[
  一串可逆变换 $k = cal(F)_"KZ" (z)$，把简单参考分布 $p_Z (z)$（通常是标准高斯）
  映成目标分布 $p_K (k)$。密度的变换规则是多元微积分的变量替换公式：

  $ p_K (k) = |det (partial cal(F)_"KZ" (k))/(partial k)| p_Z (z), quad z = cal(F)_"KZ"^(-1)(k) . $

  逆变换把 $k$ 拉回参考空间（“normalize” 就是这么来的：目标分布被拉平成高斯）。
  对数似然因此可算：

  $ log p_K (cal(D) | theta, phi) = sum_(i=1)^N log p_Z (z_i) - log |det (partial cal(F)_"ZK" (z_i))/(partial z)| . $

  $n$ 层复合 $cal(F)_"ZK" = cal(F)_n compose dots compose cal(F)_1$ 时，对数密度逐层累加：

  $ log p_K (k) = log p_(Z^0)(z^0) - sum_(j=1)^n log |det (partial cal(F)_j (z^(j-1)))/(partial z^(j-1))| . $

  全部的技术难度集中在一处：每一层的雅可比行列式必须*便宜*。任意神经网络不行，
  需要专门设计的可逆结构。
]

#term("RealNVP", "Real-valued Non-Volume Preserving") 耦合块是最常用的可逆结构。
把输入拆成两半 $(overline(z), tilde(z))$，只变换其中一半：

$ overline(k) = overline(z), quad tilde(k) = tilde(z) dot exp(s(overline(z))) + t(overline(z)) , $

其中 $s, t$ 是任意神经网络（不需要可逆）。雅可比矩阵是下三角的：

$ cal(J) = mat(I_m & 0; (partial tilde(k))/(partial overline(z)) & "diag"(exp(s(overline(z))))) , $

行列式只在对角线上：$det cal(J) = exp(sum s(overline(z)))$，一个求和就够。
逆变换直接解出：$tilde(z) = (tilde(k) - t(overline(k))) dot exp(-s(overline(k)))$。
奇数层变换前半、偶数层变换后半，交替进行，保证每个维度都被变换过。

#figure(
  align(center, image("../fig_b/b3_realnvp.png", width: 88%)),
  caption: [RealNVP 耦合块与多层 flow（论文 2108.12956）。每块只动一半变量，
  缩放因子 $exp(s)$ 与平移 $t$ 由另一半通过神经网络算出。
  因为被动的部分不经过任何网络，雅可比是三角阵，行列式只在对角线上。
  多层交替堆叠后，简单高斯被逐块塑造成任意复杂的目标分布。],
)

== NFF 模型：给随机场配上 flow

现在把 flow 用在随机场上。目标：建模条件密度 $p_K (k | x)$，
$k$ 是随机场在位置 $x$ 的取值。分两步。

*第一步：构造参考高斯随机场* $z(x, omega)$，它带 KL 结构但由网络参数化：

$ z(x, omega) = A(x) + B(x) xi(omega) + "diag"(C(x)) epsilon(x, omega) , $

$A in RR^D$（均值场）、$B in RR^(D times M)$（KL 模态，$M$ 是截断阶数）、
$C in RR^D$（对角噪声的标准差）都是神经网络，$xi tilde N(0, I_M)$，
$epsilon tilde N(0, I_D)$ 独立同分布。对角噪声 $epsilon$ 的作用是防止
任意传感器组合下的协方差矩阵奇异（退化成秩 $M$ 的矩阵不可逆）。

*第二步：用可逆变换把参考场搬成目标场：*

$ k(x, omega) = cal(F)_"ZK" (z(x, omega), x) . $

变换以 $x$ 为条件，所以同一个变换器可以在不同位置产生不同的局部行为。

这一步的巧妙之处在密度的*解析性*。给定一个 snapshot（$N$ 个位置的观测），
参考场的协方差矩阵是低秩加对角的：

$ bold(Sigma)^s = mat(B(x_1^s); dots.v; B(x_N^s)) mat(B(x_1^s); dots.v; B(x_N^s))^T + "diag" mat(C(x_1^s); dots.v; C(x_N^s))^2 , $

它的对数密度用 Woodbury 恒等式与行列式引理化成 $M times M$ 的矩阵运算，
$N$ 很大时也便宜。目标场的对数密度再叠加逐点的雅可比对数行列式：

$ log P(k_(1:N)^s | x_(1:N)^s) = log cal(N)(z_(1:N)^s | bold(mu)^s, bold(Sigma)^s) + sum_(i=1)^N log |det (partial F_"KZ" (k_i^s, x_i^s))/(partial k_i^s)| . $

数据损失就是对 snapshot 平均的负对数似然。整条链上没有任何近似：
高斯密度精确、雅可比精确、乘法律精确。这是它相对 GAN 的根本差别。

== 后验推断有闭式解

正向建模之外，反向问题需要“给定观测，隐变量 $xi$ 的后验分布”。
由于参考场是高斯的，这个后验*解析可写*：

$ bold(mu)_xi = B(bold(X))^T (B(bold(X)) B(bold(X))^T + "diag"(C(bold(X)))^2)^(-1) tilde(bold(Z)) , $

$ bold(Sigma)_xi = I_M - B(bold(X))^T (B(bold(X)) B(bold(X))^T + "diag"(C(bold(X)))^2)^(-1) B(bold(X)) , $

其中 $tilde(bold(Z))$ 是观测通过逆变换拉回参考空间的结果。
直觉读法：这是高斯线性回归的标准公式，$B$ 扮演设计矩阵，
$"diag"(C)^2$ 扮演噪声方差。观测越多、越靠近训练分布，
后验方差越小；反之 $bold(Sigma)_xi$ 趋近 $I_M$，即“什么都不知道”，退回先验。

== 物理损失：变分形式

把物理约束加进 flow 框架时有一个技术障碍：
flow 的输出是随机场的样本，方程残差要对它求导，而 flow 网络需要光滑激活函数，
训练昂贵。论文的解法是采用*变分（弱）形式*的方程损失：
把方程两边乘一个测试函数 $h$ 再积分，分部积分把导数转移到 $h$ 上
（$h$ 是已知的盒状函数，导数解析可得），于是 $u$ 只需要函数值，不需要导数：

$ cal(L)_"equ" = E_(omega, c) [ ( ⟨k nabla u, nabla h⟩ - ⟨f, h⟩ )^2 ] , $

测试函数取盒状 $h(x, c, r) = r^(-D_x) bold(1)_(norm(x - c)_infinity <= r\/2)$，
中心 $c$ 随机采样、支撑始终含在物理区域内。
平方项的期望直接用样本平均估计会有偏差（平方的期望不等于期望的平方），
论文用两组独立采样 $e_i, e_i'$ 的乘积平均 $hat(cal(L))_"equ" = -1/n sum e_i e_i'$
来保证估计无偏。这个“独立双采样”技巧在下一篇（MC-PINN）里会再次出现并成为主角。

#figure(
  align(center, image("../fig_b/b3_schematic.png", width: 92%)),
  caption: [PI-NFF 总框架（论文 2108.12956）。三个独立的 flow 分别建模随机系数 $k$、
  解 $u$、源项 $f$；损失由三部分组成：传感器数据的负对数似然（严格密度）、
  变分形式的方程残差、边界条件。正向、反向、混合问题用同一套代码：
  区别只在哪些量有观测数据。],
)

== 实验结论

模型能力实验：对指数高斯场（强非高斯）、混合双模态场（一半概率加正弦、一半概率减正弦），
生成的 $10^5$ 条路径的均值、标准差与协方差谱都与真值吻合；
单点测量的后验推断给出正确的条件分布。
方程求解实验：一维与二维随机椭圆方程的正向、反向、混合问题误差与蒙特卡洛参考解同量级。

#warn("NFF 换来的与放弃的")[
  换来的：密度严格可算、训练稳定（负对数似然是温和的目标）、
  基函数个数固定为 KL 截断 $M$ 不随非线性增长、传感器位置可以在 snapshot 之间变动。
  放弃的：模态展开“每个模态有物理含义”的可解释性；以及（下一节会看到）
  当物理方程本身是分数阶时，无论模态还是 flow，都还有另一个障碍在等着：
  分数阶导数没法自动微分。
]
