#import "../../preamble.typ": *

= 第一篇：把随机性放进 PINN（NN-aPC，2019）

== 它要解决的问题

PINN 原始版本只处理确定性方程。Zhang, Lu, Guo, Karniadakis（JCP 2019）第一次问：
方程本身带随机性时，PINN 怎么办？他们把“总不确定度”拆成两层：

- *参数不确定度*：方程系数 $k(x, omega)$ 是随机场，解 $u(x, omega)$ 因此也是随机场；
- *近似不确定度*：神经网络本身没拟合好，带来预测误差。

第一层靠下面的模态展开处理，第二层靠 dropout 处理（见本章末尾）。

直接把随机变量 $xi$ 和物理坐标 $x$ 一起塞进网络输入是行不通的：
这样做不给出“解的随机性”的显式表达（你没法从网络输出读出各阶统计量），
而且随机样本没法像坐标那样在配点上标注。论文的做法是先把随机场降维，
再把多项式混沌展开搭进网络结构。

== 从传感器数据出发的 KL 降维

设有 $N$ 条历史记录（snapshot），每条给出 $k$ 在传感器位置 $x_k^((i))$ 的读数。
把传感器读数堆成向量 $bold(k)_s$，计算样本协方差矩阵

$ K_(i j) = op("Cov")(k^((i)), k^((j))), $

$K$ 的第 $(i, j)$ 个元素是位置 $i$ 与位置 $j$ 的读数在 $N$ 条记录上的协方差；
样本平均代替真实期望的合法性由大数定律给出（第 0 章的蒙特卡洛段落）。
这个三步展开把从 $K$ 到展开式的全部距离走完。

*第一步（对角化）*：$K$ 是对称半正定矩阵（对称性来自协方差的定义；
半正定性来自 $bold(c)^T K bold(c) = op("Var")[sum_i c_i k^((i))] >= 0$），
所以可以正交对角化 $K = Phi^T Lambda Phi$，$Phi$ 的行是单位特征向量、
互相正交，$Lambda$ 的对角元 $lambda_1 >= lambda_2 >= ... >= 0$ 是特征值。

*第二步（白化坐标）*：定义

$ bold(xi)_s = Phi^T sqrt(Lambda)^(-1) bold(k)_s , $

理由算一下新坐标的协方差就清楚：
$E[bold(xi)_s bold(xi)_s^T] = Phi^T sqrt(Lambda)^(-1) K sqrt(Lambda)^(-1) Phi = Phi^T Phi = bold(I)$
（中间把 $K = Phi^T Lambda Phi$ 代入，$sqrt(Lambda)^(-1) Lambda sqrt(Lambda)^(-1) = bold(I)$，
最后用 $Phi$ 的行正交）。也就是说 $xi_(s, 1), ..., xi_(s, D)$ 互不相关、
各自方差为 1：原来的相关坐标被换成了不相关坐标。

*第三步（反解重构）*：把第二步的定义两边左乘 $sqrt(Lambda) Phi$
（用 $Phi Phi^T = bold(I)$）解出 $bold(k)_s = Phi sqrt(Lambda) bold(xi)_s$，
写成分量并把样本均值部分 $k_0$ 加回去：

$ k(x_k^((i)); omega_s) approx k_0(x_k^((i))) + sum_(l=1)^M sqrt(lambda_l) k_l (x_k^((i))) xi_(s, l) . $

近似号来自截断：求和只保留特征值最大的前 $M$ 项，丢掉的都是方差很小的方向。
这一步是纯线性代数：不需要假设协方差函数的形式，样本协方差就是全部输入。
$M$ 取截断阶数，随机维数从“无穷”降到 $M$。

== aPC：不需要分布假设的正交基

接下来要把解 $u(x, omega)$ 展开成 $xi$ 的正交多项式。经典 gPC 的正交基由
$xi$ 的分布决定（高斯分布配 Hermite 多项式、均匀分布配 Legendre 多项式）。
但数据驱动场景里 $xi$ 的分布是未知的：传感器数据不保证是高斯的。

#term("aPC", "任意多项式混沌") 的做法是：
不假设分布，直接用*经验测度*构造正交基。把 $N$ 条记录的隐变量
${bold(xi)_s _(s=1)^N}$ 当作经验分布

$ rho(bold(xi)) approx nu_S (bold(xi)) = 1/N sum_(s=1)^N delta_(bold(xi)_s)(bold(xi)) , $

然后在离散测度 $nu_S$ 上做 Gram–Schmidt 正交化。所有运算只用到函数在
$N$ 个数据点上的取值，所以先定义数据上的内积（下文反复出现）：

$ ⟨f, g⟩_S := 1/N sum_(s=1)^N f(bold(xi)_s) g(bold(xi)_s) , $

即“在经验测度上的平均乘积”。单项式先排好序（按总次数从低到高，
同次内部顺序任意但固定），$psi_alpha^* (bold(xi)) = product_(i=1)^M xi_i^(alpha_i)$
是候选单项式，$alpha = (alpha_1, ..., alpha_M)$ 是非负整数组，
$beta prec alpha$ 表示 $beta$ 排在 $alpha$ 之前。正交化按次序逐个进行：
处理第 $alpha$ 个候选时，$prec alpha$ 的那些已经正交化完毕。递推为

$ psi_alpha (bold(xi)) = w_alpha^alpha psi_alpha^*(bold(xi)) - sum_(beta prec alpha) w_beta^alpha psi_beta (bold(xi)) , $

系数不用解方程组，逐个算内积就够。推导如下。

*第一步（正交条件）*：要求 $⟨psi_alpha, psi_gamma⟩_S = 0$ 对每个已入选的
$gamma prec alpha$。把递推式代入，用已完成的正交性
$⟨psi_beta, psi_gamma⟩_S = delta_(beta gamma)$：

$ 0 = w_alpha^alpha ⟨psi_alpha^*, psi_gamma⟩_S - sum_(beta prec alpha) w_beta^alpha delta_(beta gamma) = w_alpha^alpha ⟨psi_alpha^*, psi_gamma⟩_S - w_gamma^alpha , $

所以 $w_gamma^alpha = w_alpha^alpha ⟨psi_alpha^*, psi_gamma⟩_S$：每个系数是一个
已经能算的内积乘公共因子。（原文说“解一个小线性方程组”，是因为把
全部正交条件联立起来看；按递推次序看，方程组是三角的，回代就是上式。）

*第二步（归一化）*：要求 $⟨psi_alpha, psi_alpha⟩_S = 1$ 定出公共因子。
把第一步的系数代回递推式展开（三次用到 $⟨psi_beta, psi_gamma⟩_S = delta_(beta gamma)$
把双重求和塌缩成单个求和）：

$ ⟨psi_alpha, psi_alpha⟩_S = (w_alpha^alpha)^2 [ ⟨psi_alpha^*, psi_alpha^*⟩_S - sum_(gamma prec alpha) ⟨psi_alpha^*, psi_gamma⟩_S^2 ] = 1 , $

开方取正号即 $w_alpha^alpha$。两步合起来，正交条件

$ 1/N sum_(s=1)^N psi_alpha (bold(xi)_s) psi_beta (bold(xi)_s) = delta_(alpha, beta) $

对*这批数据*精确成立：正交化只触及函数在 $N$ 个点上的取值，
不涉及分布的任何假设。基函数个数是 $P + 1 = (r + M)! \/ (r! M!)$，$r$ 是阶数。

#insight("用日常语言读一遍")[
  aPC 的正交化在做一件事：把相关的单项式特征改造成互不相关的特征，
  改造标准完全由这批数据自己决定。第一步保证新特征与所有旧特征
  在数据上不相关，第二步把尺度统一。因为整个过程只查数据点上的
  取值表，分布是高斯、多峰还是别的，对算法没有区别；
  代价是正交性只对这批数据成立，换一批数据要重算。
]

#keypoint("aPC 换来了什么，代价是什么")[
  换来的：正交基精确匹配数据的实际分布（哪怕非高斯、多峰），
  展开系数的计算稳定（正交性保证模态互不干扰）。
  代价没有消失：基函数个数仍随 $M$ 和 $r$ 指数增长。
  这个限制到第 8 章（NFF）才被解除，那里会把基的规模固定在截断阶 $M$ 上。
]

== 网络结构与损失函数

解的展开写成

$ u(x, omega) = sum_(alpha=0)^P u_alpha (x) psi_alpha (bold(xi)) , $

$u_alpha (x)$ 是*模态*：第 $alpha$ 个多项式方向的系数，是物理坐标的确定性函数。
NN-aPC 用两个独立的深度网络分别输出 $k$ 的模态和 $u$ 的模态：

- 网络 $wide(k_i)$ 输出 $k$ 的 $M + 1$ 个模态 $k_0, sqrt(lambda_1) k_1, ..., sqrt(lambda_M) k_M$；
- 网络 $wide(u_alpha)$ 输出 $u$ 的 $P + 1$ 个模态 $u_0, u_1, ..., u_P$。

重构规则是

$ tilde(k)(x; omega_s) = wide(k_0)(x) + sum_(i=1)^M sqrt(lambda_i) wide(k_i)(x) xi_(s, i) , $

$ tilde(u)(x; omega_s) = sum_(alpha=0)^P wide(u_alpha)(x) psi_alpha (bold(xi)_s) . $

#figure(
  align(center, image("../fig_b/b1_nn_inverse_sketch.png", width: 88%)),
  caption: [NN-aPC 解随机椭圆方程的结构（论文 1809.08327）。上下两支网络分别
  输出 $k$ 模态与 $u$ 模态；$xi_s$ 来自观测数据的 KL 降维；
  重构 $tilde(k), tilde(u)$ 代入随机方程得到残差 $f tilde(f)$，
  与三类数据的失配一起构成损失。一次训练同时覆盖所有 snapshot。],
)

损失函数对每个 snapshot 平均：

$ cal(L)(cal(S)_t) = "MSE"_u + "MSE"_k + "MSE"_f , $

$ "MSE"_u = 1/(N N_u) sum_(s=1)^N sum_(i=1)^N_u [tilde(u)(x_u^((i)); omega_s) - u(x_u^((i)); omega_s)]^2 , $

$ "MSE"_f = 1/(N N_f) sum_(s=1)^N sum_(i=1)^N_f [cal(N)_x[tilde(u)(x_f^((i)); omega_s); tilde(k)(x_f^((i)); omega_s)]]^2 . $

$"MSE"_k$ 与 $"MSE"_u$ 同构，惩罚 $k$ 的观测失配。三个术语的含义：
$"MSE"_u$ 要求重构的解贴合观测；$"MSE"_k$ 要求重构的参数贴合观测；
$"MSE"_f$ 要求*随机方程本身在每个 snapshot 上成立*，这是物理约束。

两个实践技巧：第一，均值模态 $u_0$ 单独用一个小网络学，
因为均值和其他模态的尺度差很多。第二，同阶的 aPC 模态分进同一个网络，
因为同阶模态尺度接近，跨尺度混在一起会让梯度被大尺度模态主导。

#figure(
  align(center, image("../fig_b/b1_nn_details.png", width: 88%)),
  caption: [模态网络的分组结构（论文 1809.08327）。$k$ 的各阶模态与 $u$ 的各阶模态
  分别分组输出。分组的依据是模态振幅的尺度：把尺度相近的模态放进同一个网络，
  避免优化时大尺度模态淹没小尺度模态。],
)

== 近似不确定度：dropout

参数不确定度（$k$ 的随机性）由模态展开显式表达。
近似不确定度（网络没学好）用 #term("MC dropout", "蒙特卡洛丢弃") 估计：
训练时随机丢弃一部分神经元，推断时也保持丢弃并采样 $T$ 次，统计预测的均值与方差：

$ E(y) approx 1/T sum_(t=1)^T cal(N)_t (x) . $

它的理论依据是 Gal 与 Ghahramani 的结果：带 dropout 的网络在数学上
近似于一个贝叶斯神经网络（权重带先验、后验取变分近似）。
采样的波动就是“网络对自己的预测有多确定”的读数。

dropout 在这里有一个超出“报误差棒”的用法：*主动学习*。
网络在哪个位置不确定度最大，就说明那里最缺数据；
把新的传感器（或新的方程求解）投放到那个位置，信息收益最大。
论文的实验里，用这个准则迭代 15 步只新增 8 个传感器，
反问题的参数预测误差就从 5% 以上压到 1% 以下。

#figure(
  align(center, grid(columns: 3, gutter: 8pt, image("../fig_b/b1_al_first.png"), image("../fig_b/b1_al_last.png"), image("../fig_b/b1_al_error.png"))),
  caption: [主动学习的过程（论文 1809.08327）。左：初始 4 个 $k$ 传感器的预测，
  外推区（两端）不确定度大；中：按 dropout 不确定度最大处补传感器后的预测；
  右：误差随迭代步的下降。不确定度大的位置恰好就是误差大的位置，
  这是“不确定度可以当作误差的代理”的直接证据。],
)

== 实验结论与遗留问题

正向（随机 Poisson，$f$ 为高斯随机场）与反向（识别非高斯随机场 $k$）都验证了框架；
二阶 aPC 在反向问题里把 $k$ 的标准差误差从 5.26% 压到 1.87%。
另一个观察：即使每个 $k$ 传感器的均值和标准差几乎相同，
训练出的 $k$ 网络仍能重构出全场非平凡的波形。
这说明三类训练数据与随机方程之间的耦合是*信息融合*，
网络学到的东西超出了插值器的范畴：它回答“什么参数能产生这样的解”。

#warn("NN-aPC 留下的两个问题")[
  *第一，基函数个数随有效维数指数增长。* $M = 19$ 时一阶 aPC 就有 20 个模态，
  $r = 2$ 时是 210 个。相关长度短的随机场会把 $M$ 推得更大。
  *第二，gPC 类展开的基不随时间演化*，所以它只能处理稳态方程或固定分布。
  时间依赖的随机方程（例一）需要基函数本身跟着时间走。
  这两个问题分别由第 8 章的 NFF 和第 6 章的模态空间方法接手。
]
