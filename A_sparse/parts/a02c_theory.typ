#import "../../preamble.typ": *

= 理论深处：三条主定理的证明思路

主章节把六篇的定理当黑盒引用。这一章把其中三条的证明骨架摊开，
摊到“每一步在做什么”的粒度，不展开纯技术引理。读完后，
概率路线的论文（第一、六篇）就不再是不可审计的公告。

== RIP 从哪来：一致界 + 浓缩

Rauhut–Ward 型定理（$M >= C delta^(-2) L s log^3 s log N$ 概率成立）
的证明分三层，第一篇的附录是这套流程在“Gauss 网格子采样”上的完整执行。

*第一层：逐点一致界。* 要让 $bold(D)\/sqrt(M)$ 近似等距，
先要求它的元素不至于爆：对加权的正交多项式系统，
存在与指标无关的逐点界

$ sup_x abs(psi_(k, n)(x))^2 <= L_i (n) , $

$L_i (n)$ 就是主定理里的分布常数（Jacobi：常数；指数型：$C n^(2\/3)$）。
第一篇附录的推导（Mhaskar–Rahmanov–Saff 数定位节点区间 + 加权多项式
一致界）就是在对三类分布分别估计这个 $L_i$。所有后续都是概率机器，
这一层是唯一的“硬分析”。

*第二层：矩阵 Concentration。* 把 $norm(bold(D) bold(v))_2^2$ 写成
独立采样量的和：每个样本贡献一个秩一矩阵 $psi(bold(x)_m) psi(bold(x)_m)^T$。
其期望恰好是单位阵（离散正交系统的定义），于是
$norm(bold(D) bold(v))_2^2$ 集中在 $M norm(bold(v))_2^2$ 附近。
集中的宽度由两个量决定：单个秩一项的范数界（第一层给的 $L_i$）
与有效维数 $N$。矩阵 Bernstein 或 Mendelson 的小球法给出：
样本量 $M$ 达到 $L s log^3 s log N$ 量级时，对所有 $s$-稀疏向量
（一个 组合规模的族（$binom(N, s)$ 个向量），所以有 $log N$ 与 $log s$ 因子）
一致集中的概率超过 $1 - N^(-gamma log^3 s)$。

*第三层：从行子采样到 RIP。* 满网格的加权矩阵是正交矩阵
（离散正交系统的张量化），从中无放回抽 $M$ 行：
“正交矩阵的随机行子集保持 RIP”是压缩感知的标准引理
（Rauhut 教材的 Corollary 12.38），其证明同样是第二层的集中论证
加上行范数一致界。

读这份骨架的收获：$L$ 因子出现在*第一层*，所以它对样本量的影响
是乘法性的、不可通过更多采样弥补的；$log$ 因子出现在*第二层*，
是“对所有稀疏向量一致集中”的代价，换更精细的论证（比如
 Mendelson 的小球法）只能改善常数与幂次。

== 加权最小二乘的无偏性：一笔与 induced 采样相同的账

第五篇的加权 LS（采样自 $kappa rho$、权重 $w = K = 1\/kappa$）
与第六篇的 induced 采样（采样自经验 $kappa omega_hat$、权重 $1\/kappa$）
是同一笔账的连续与离散版本。

#proposition("命题（加权 LS 的无偏性，差一个常数）")[
  设采样密度 $tilde(rho)(z) = kappa(z) rho(z) \/ c$，
  $c = integral kappa rho$（归一化常数），权重 $w(z) = 1\/kappa(z)$。
  对任意 $h$：

  $ E_(z tilde tilde(rho))[w(z) h(z)] = integral (kappa rho \/ c)(h \/ kappa) d z = (1/c) integral h rho d z . $

  即加权平均等于真期望乘 $1\/c$；常数与 $h$ 无关，
  被最小二乘的尺度不变性吸收。
]

#proof[
  两个密度相除：$w tilde(rho) = (1\/kappa)(kappa rho \/ c) = rho \/ c$。
  对 $h$ 积分即得。逐行成立，对矩阵形式
  $E[w(bold(z)) bold(psi)(bold(z)) bold(psi)(bold(z))^T] = (1\/c) bold(I)$
  同样成立（逐元素取 $h = phi_j phi_k$）。$square$
]

这条命题同时解释了两个设计选择。第一，为什么权重必须是
$1\/kappa$（而不是别的）：只有它能消掉采样偏置。第二，
为什么采样要偏向 $kappa$ 大（能量高）的区域：这样 $c$ 有限且
估计的方差被 $w$ 的有界性控制。若反过来从 $1\/kappa$ 小的区域
多采样，$w$ 会在那里爆炸，方差失控。
第五篇的“平衡测度”（$kappa rho$ 的极限）与第六篇的 induced 测度
（$kappa omega_hat$ 的离散版）因此是同一个提案分布，
一个取了位势论极限、一个停在有限样本。

== 稀疏性从哪来：随机 PDE 解的系数衰减

压缩感知配点法成立的前提是“解在多项式基下稀疏”。
这不是经验观察，对椭圆方程有定理级的依据
（Cohen–DeVore–Schwab 2010 起的分析；Babuška–Nobile–Tempone 的
早期数值证据）。

#proposition("定理（稀疏性，叙述形式）")[
  设扩散系数由仿射随机场给出：
  $kappa(x, y) = overline(kappa)(x) + sum_(j >= 1) kappa_j (x) y_j$，
  其中 $sum_j norm(kappa_j)_(L^infinity)$ 一致可和
  （随机场的相关长度短时自然满足）。则解的 gPC 系数满足指数型衰减：

  $ norm(u_(bold(nu)))_(L^2) <= C e^(-c sum_i nu_i) quad "对某个与" bold(nu) "无关的" c > 0 . $

  即按系数向量的“总次数”衰减，重要的系数集中在低次数、
  少数活跃维度上。
]

读法：衰减率 $c$ 由随机场的*能量谱*决定，$norm(kappa_j)$ 衰减越快
（相关长度越短，KL 特征值衰减越快），解越稀疏，压缩感知越省样本。
这也解释了第一篇到第六篇实验里反复出现的模式：
相关长度长（$l_c$ 大）的场难恢复，相关长度短（$l_c$ 小）的场好恢复。
反向推论同样有用：如果实际问题里输入场的相关结构未知且可能平缓
（$sum norm(kappa_j)$ 接近临界），稀疏假设本身要先用数据检验
（第六篇的“作者如实说明高维低次时 induced 与 MC 相当”就是这种情形）。

== 互不相干路线的逐点界

第三篇走的 MIP 路线不需要概率集中，它的全部内容是逐点估计。
对加权后的测量矩阵 $hat(bold(Phi)) = bold(W) tilde(bold(Phi)) bold(P)$，
互不相干常数由逐列上界控制：

$ mu(hat(bold(Phi))) <= max_i sup_z norm(hat(bold(Phi))_(i)(z))_2^2 / (min_i norm(hat(bold(Phi))_(i))_2) . $

分子的逐点界由“Jacobi 导数恒等式 + 加权一致界”给出
（第三篇定理里的 $2 e (2 + sqrt(alpha^2 + beta^2))$ 因子），
分母由列归一矩阵 $bold(P)$ 抬高。这套论证与第一层一致界是同一技术
（加权多项式的逐点估计）在不同矩阵上的应用。
这也是为什么六篇论文共享同一个数学工具箱。
