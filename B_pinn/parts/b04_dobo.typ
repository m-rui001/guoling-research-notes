#import "../../preamble.typ": *

= 第二篇：模态空间（NN-DO/BO，2020）

== 新困难：解要随时间演化

第 4 章的展开 $u(x, omega) = sum_alpha u_alpha (x) psi_alpha (bold(xi))$ 有一个隐藏假设：
模态 $u_alpha$ 只依赖空间，随机基 $psi_alpha$ 不随时间变。稳态方程没问题；
但时间依赖的随机偏微分方程（SPDE）里，随机结构的能量会重新分配：
一开始不重要的方向后来可能变成主导，两个方向的方差曲线还会相撞。
固定的基在长时间积分中需要非常多的模态才能跟上。

Zhang, Guo, Karniadakis（SISC 2020）采用 Sapsis 提出的*广义时变 KL 展开*：
让空间基和随机系数同时演化。截断到 $N$ 项：

$ u_N (x, t; omega) = overline(u)(x, t) + sum_(i=1)^N a_i (t) u_i (x, t) Y_i (t; omega) , $

其中 $overline(u)$ 是均值（时间依赖），$u_i (x, t)$ 是空间模态（归一化 $⟨u_i, u_i⟩ = 1$），
$Y_i (t; omega)$ 是随机模态（归一化 $E[Y_i^2] = 1$，零均值），
$a_i (t)$ 是标量振幅，吸收 $u_i$ 与 $Y_i$ 各自的尺度。

#warn("这个展开有一个自由度问题")[
  同一个 $u_N$ 可以由无穷多组 ${u_i, Y_i}$ 表出：
  把某个方向的空间模态缩小一半、随机模态放大一倍，和不变。
  没有约束的展开会把冗余摊到所有分量上，每个分量都学一半“半个真相”，
  任何分量都不再可解释。所以需要正交约束把自由度钉死。
]

== DO 与 BO：两种钉死自由度的约束

两类经典约束都能消除冗余，区别在约束加在哪一边。

- *DO（dynamically orthogonal，动力学正交）*：要求空间模态的演化方向
  垂直于它们自己张成的空间，$⟨partial_t u_i, u_j⟩ = 0$。
  直觉：$u_i$ 只在“现有模态表达不了的新方向”上演化，
  已有方向上的变化全部交给随机系数 $Y_i$。
- *BO（bi-orthogonal，双正交）*：要求两组基各自正交
  $⟨u_i, u_j⟩ = lambda_i delta_(i j)$，$E[Y_i Y_j] = delta_(i j)$。
  直觉：两组基在任意时刻都像一对配好的坐标轴，交叉内积矩阵是对角阵。

经典方法（Sapsis 2009 的 DO、Cheng–Hou–Zhang 2013 的 BO）会*显式推导*
模态的演化方程，然后数值积分。推导结果在这里值得完整写出，
因为后面 NN 版本的贡献正是绕开它们的两个致命假设。

#proposition("经典 DO 演化方程（Sapsis 2009）")[
  记 $C_(Y_i Y_j) = E[Y_i Y_j]$（随机系数协方差），
  投影 $product_(V_S^perp) F = F - sum_k ⟨F, u_k⟩ u_k$。DO 条件下模态满足

  $ (partial overline(u)(t, x))/(partial t) = E[cal(N)_x[u]] , $

  $ (d Y_i (t; omega))/(d t) = ⟨cal(N)_x[u] - E[cal(N)_x[u]], u_i⟩ , quad i = 1, ..., N , $

  $ sum_(i=1)^N C_(Y_i (t) Y_j (t)) (partial u_i (t, x))/(partial t) = product_(V_S^perp) E[cal(N)_x[u] Y_j] , quad j = 1, ..., N . $
]

#proposition("经典 BO 演化方程（Cheng–Hou–Zhang 2013）")[
  记 $S_(i j) = ⟨u_i, partial_t u_j⟩$，$M_(i j) = E[Y_i partial_t Y_j]$，
  $G_(i j) = ⟨E[cal(N)_x[u] Y_j], u_i⟩$。BO 条件下

  $ (partial overline(u))/(partial t) = E[cal(N)_x[u]] , $

  $ lambda_i (d Y_i)/(d t) = -sum_(j=1)^N S_(i j) Y_j + ⟨cal(N)_x[u] - E[cal(N)_x[u]], u_i⟩ , $

  $ (partial u_i)/(partial t) = -sum_(j=1)^N M_(i j) u_j + E[cal(N)_x[u] Y_i] , $

  并且当特征值互不相同时 $S, M$ 有闭式解

  $ M_(i j) = (G_(i j) + G_(j i))/(-lambda_i + lambda_j) quad (i != j) , quad S_(i j) = G_(i j) + lambda_i M_(i j) . $
]

#definition("两个失效模式")[
  *失效模式一（DO，协方差奇异）。* DO 方程的最后一式是一个以
  $C_(Y_i Y_j)$ 为系数矩阵的线性方程组，要求它可逆。
  而*确定性初值*（例一里最常见的情况：初始时刻没有任何随机性）下
  所有的 $Y_i$ 在 $t_0$ 都是零，协方差矩阵全零，方程组无解。经典 DO 直接失效。

  *失效模式二（BO，特征值交叉）。* BO 的闭式解分母是 $-lambda_i + lambda_j$。
  当两个模态的方差曲线 $lambda_i (t)$ 在某时刻相撞（*特征值交叉*），
  分母为零，闭式消失，显式方法在该时刻崩溃。
  长时间非线性积分里交叉几乎不可避免（6.6 节的 Burgers 实验会看到大量交叉）。
]

#figure(
  align(center, image("../fig_b/b2_bur_do_ai.png", width: 80%)),
  caption: [特征值交叉的实拍（论文 1905.01205，随机 Burgers 方程长时间积分）。
  图中 $a_1(t)$ 与 $a_2(t)$ 的曲线在时域内多次相撞：方差第一的模态和方差第二的
  模态互换身份。经典 BO 方法在每个交叉点都需要特殊处理，
  这正是 NN 版本要绕开的东西。],
)

== 关键思想：把约束隐式写进损失

NN-DO/BO 的做法：不做显式推导，把 DO/BO 条件本身当作损失函数。
网络只需要“在训练终点满足正交条件”，不需要“沿一条可积的演化方程走过去”，
于是协方差奇异与特征值交叉这两个障碍都不再出现。

推导分四步，这是本篇的理论核心。模型方程为
$partial_t u = cal(N)_x[u]$（把非齐次项收进 $cal(N)_x$）。

*第一步*：对 SPDE 两边取期望，得到均值方程

$ E[partial_t u] = partial_t overline(u) = E[cal(N)_x[u]] . $

*第二步*：两边与 $u_i$ 做物理空间内积 $⟨dot, u_i⟩$。

*第三步*：两边与 $Y_i$ 做概率空间内积 $E[dot Y_i]$。

*第四步*：把截断展开 $u_N = overline(u) + sum_i a_i u_i Y_i$ 代入，用 DO 或 BO 条件化简。

#insight("为什么这三条投影就是弱形式")[
  $u_i (x, t)$ 构成物理子空间 $V_S$ 的基，$Y_i (t, omega)$ 构成随机子空间
  $Omega_S$ 的基。对基函数做内积投影，正是偏微分方程理论里“弱形式”的标准构造：
  方程不必在每点逐点成立，只需在测试函数（这里是基函数）上的投影为零。
  第一步取期望相当于对随机空间里的“常数 1”（第 0 个基）做投影。
  三条投影合起来，就提供了在乘积空间 $V_S times Omega_S$ 中确定 $u_N$ 的全部信息。
]

网络用四个模块实现展开式的四个成分：

$ u_(n n)(x, t; xi) = overline(u)_(n n)(x, t) + sum_(i=1)^N A_(n n, i)(t) U_(n n, i)(x, t) Y_(n n, i)(xi, t) . $

$overline(u)_(n n)(x, t)$ 是均值网络；$U_(n n, i)(x, t)$ 是空间模态网络（$N$ 组输出）；
$Y_(n n, i)(xi, t)$ 是随机模态网络，输入 $xi$ 是标准随机变量的样本；
$A_(n n, i)(t)$ 是振幅网络。

#figure(
  align(center, image("../fig_b/b2_dobo_sketch.png", width: 92%)),
  caption: [NN-DO/BO 的计算图（论文 1905.01205）。左侧三个块是弱形式损失的三个投影
  （对 $omega$ 取期望、与 $U_i$ 做物理内积、与 $Y_i$ 做概率内积）；
  右侧的块是 DO 或 BO 的正交约束损失；底部是初值与边界损失。
  所有块都是可微的，联合训练。],
)

== 损失函数（完整）

弱形式损失：在配点 $(x_c^k, t_c^s)$ 与随机样本 $xi_c^l$ 上评估三条投影的均方误差：

$ "MSE"_w = 1/(n_x n_t) sum_(k, s) (epsilon_1^(k s))^2 + 1/(n_t n_xi) sum_(s, l) (epsilon_2^(s l))^2 + 1/(n_x n_t) sum_(k, s) (epsilon_3^(k s))^2 , $

$ epsilon_1^(k s) := E[ partial_t u_(n n)(x_c^k, t_c^s; xi) - cal(N)_x[u_(n n)(x_c^k, t_c^s; xi)] ] , $

$ epsilon_2^(s l) := ⟨partial_t u_(n n)(x, t_c^s; xi_c^l) - cal(N)_x[u_(n n)], U_(n n, i)(dot, t_c^s)⟩ , $

$ epsilon_3^(k s) := E[ (partial_t u_(n n) - cal(N)_x[u_(n n)]) Y_(n n, i)(t_c^s; xi) ] . $

正交约束损失（以 BO 为例；DO 版把中间项换成 $⟨dot U_i, U_j⟩$ 的对应形式）：

$ "MSE"_"BO" = 1/(N n_t) sum_(i, s) (E[Y_i])^2 + 1/(N^2 n_t) sum_(i, j, s) (⟨dot(U)_i, U_j⟩ + ⟨dot(U)_j, U_i⟩)^2 + 1/(N n_t) sum_(i, j, s) (E[Y_i dot(Y)_j] + E[Y_j dot(Y)_i])^2 . $

逐项读：第一项惩罚 $E[Y_i] = 0$（随机模态零均值）；第二项惩罚空间模态的 Gram 矩阵
是对称对角的；第三项惩罚随机模态的自相关矩阵是对称单位阵。
这些量全部用 mini-batch 里的样本平均来估计，因此整个损失可微。

总损失是四块的加权和：

$ cal(L) = "MSE"_w + 100 times ("MSE"_"IC" + "MSE"_"BC" + "MSE"_"DO/BO") + 0.1 times "MSE"_0 , $

其中 $"MSE"_0$ 是原方程强形式残差的正则项（加速训练），
初值与边界损失取弱形式（边界条件在随机空间里也投影到基上）。

== 确定性初值怎么处理（绕开失效模式一）

初值是确定性的（$u(x, t_0)$ 没有随机性）时，论文的初始化是：
空间模态 $u_i (x, t_0)$ 取满足边界条件的标准正交基 $v_i (x)$，
随机模态 $Y_i (t_0; xi)$ 取单位方差的 gPC 基，振幅 $a_i (t_0) = 0$。
初始时刻没有随机能量，所以振幅全零完全合法；
随机性随后由方程的随机系数（$cal(N)_x$ 里的 $k$）逐步注入，
振幅从零开始生长。经典 DO 在这个场景崩溃的原因是它要求用协方差矩阵除一次，
而 NN 版本从头到尾没有做过这个除法。

#figure(
  align(center, image("../fig_b/b2_adv_do_ai.png", width: 80%)),
  caption: [随机平流方程的 NN-DO 结果（论文 1905.01205）。方程 $u_t + xi u_x = 0$，
  $xi tilde N(0, 0.8^2)$，确定性初值 $u(x, 0) = -sin x$。振幅 $a_i (t)$ 从零开始
  增长：随机性由随机平流速度逐步注入。该例有精确的 DO 解
  （如 $u_1^"DO" = -1/sqrt(pi) cos x$，$Y_1^"DO" = -sqrt(pi) sin(xi t)$），
  NN-DO 各分量与精确解的相对 $L^2$ 误差在 2% 以内。],
)

== 实验与遗留问题

三个实验递进。平流方程（有精确解）验证正确性；随机 Burgers 方程积分到
$T = 10 pi$（时域切成 10 段逐段训练），人工构造解保证大量特征值交叉出现，
NN-DO 与 NN-BO 的全部分量相对误差不超过 1%，而经典 BO 在交叉处崩溃；
非线性扩散反应方程带 19 维随机输入（KL 19 阶），6 个 BO 模态就抓住大部分随机能量，
对 $sigma = 0.1$ 的噪声初值稳健，反向问题同时识别两个反应系数
（真值 0.5 与 0.3），10 万 epoch 内收敛。

#warn("模态空间方法的边界")[
  两个遗留问题。第一，训练成本高：$N$ 个模态、每个都要网络输出，
  绝对误差很难压到 $10^(-5)$ 以下（非凸优化的天花板）。第二，
  展开的基仍然受限于“模态”这个形式：对非高斯、多峰的随机场，
  正交多项式或正交模态的表达效率不高，模态个数还是要跟着有效维数走。
  如果我们不想要“模态”，只想要“解的精确分布”，需要另一类工具，见下一节。
]
