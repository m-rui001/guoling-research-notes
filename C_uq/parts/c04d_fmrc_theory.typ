#import "../../preamble.typ": *

= FMRC 深处：算子误差、可计算准则与一次验证

FMRC 这一篇在四篇后继里最“数学”，也最短。它把一个降维问题
写成两个条件密度模型的拟合问题，再证明拟合损失控制算子误差；
实验只有一个七井势。这个配比决定了读法：定理链条要逐步走完，
实验部分则要清楚它到底验证到了哪一步。

== 设定：允许非平衡的转移算子

相空间 $cal(X) subset RR^D$ 带 Lebesgue 测度 $lambda$。
给定 lag $tau$，转移概率与转移密度

$ P_tau (x, B) = op("Prob") { X_(t+tau) in B | X_t = x } , quad P_tau (x, d y) = p_tau (x, y) d y , $

手上有的数据是 lag-$tau$ 的对 ${ (x_n, y_n) }_(n=1)^N$，
其联合密度写成

$ rho(x, y) = rho_0(x) p_tau (x, y) , quad rho_1(y) = integral rho_0(x) p_tau (x, y) d x . $

这里 $rho_0$ 与 $rho_1$ 允许不同，也就是*不假设过程处于平衡态*。
平衡是特例 $rho_0 = rho_1 = pi$。为什么要费这个力气：
采样数据往往来自一段有限轨迹，起点的分布并不等于不变测度，
按平衡假设写的理论在这样的数据上不成立。

两个算子分别按密度前传与观测量前传：

$ (cal(T)_tau u_t)(y) = integral rho_0(x) \/ rho_1(y) p_tau (x, y) u_t (x) d x
    = E[ u_t (X_t) | X_(t+tau) = y ] , $

$ (cal(K)_tau f)(x) = integral p_tau (x, y) f(y) d y
    = E[ f(X_(t+tau)) | X_t = x ] , $

前者是 Perron–Frobenius 算子，作用在 $L^1_(rho_0) -> L^1_(rho_1)$；
后者是 Koopman 算子，在 $L^oo$ 上。
作者特别指出与 Bittracher 等的定义差别：
那里 Perron–Frobenius 是相对不变测度写的，这里相对经验初值测度 $rho_0$ 写。
这一个改动决定了后面所有范数里出现的权重是 $rho_0, rho_1$。

反应坐标 $r : cal(X) -> RR^d$（$d << D$，$C^1$），
水平集 $Sigma_r (z)$ 是 $D - d$ 维子流形。

== 两个准则，以及它们为什么成对出现

#definition("lumpability（可归并）")[
  存在 $p_tau^L : cal(Z) times cal(X) -> RR^+$ 使 $p_tau (x, y) = p_tau^L (r(x), y)$ 对所有 $x, y$。
  同一个水平集 $r^(-1)(r(x))$ 里的所有点有完全相同的未来分布。
]

#definition("decomposability（可分解）")[
  存在 $p_tau^D : cal(X) times cal(Z) -> RR^+$ 与 $p_"local" : cal(Z) times cal(X) -> RR^+$
  使 $p_tau (x, y) = p_tau^D (x, r(y)) dot p_"local" (r(y), y)$。
  转移分两步：先用 $p_tau^D$ 决定未来落在哪个水平集，
  再用只依赖目标集与目标点的 $p_"local"$ 在集内定位。
]

两条不是同一个性质的正反说法，它们各自管一半信息。
lumpability 说“把过去压缩到 $r(x)$ 不丢预测未来的能力”；
decomposability 说“未来状态可以由 $r(y)$ 加一层与过去无关的局部噪声生成”。
只有两条同时成立，$r$ 才是对动力学无损的坐标。

反向密度是 $p_(-tau)(y, x) = rho_0(x) p_tau (x, y) \/ rho_1(y)$，
它是 Bayes 反演意义下的“从 $y$ 回到 $x$”的转移密度。
论文用一条引理把 decomposability 换成反向的陈述：

#lemma("引理（反向等价）")[
  $p_tau$ 关于 $r$ decomposable 当且仅当存在
  $p_(-tau)^D : cal(Z) times cal(X) -> RR^+$ 使 $p_(-tau)(y, x) = p_(-tau)^D (r(y), x)$。
]

#proof[
  必要性，分三步。
  *第一步：代入反向密度的定义。* 反向密度是 Bayes 反演：
  $p_(-tau)(y, x) = rho_0(x) p_tau (x, y) \/ rho_1(y)$
  （联合密度 $rho(x, y) = rho_0(x) p_tau (x, y)$ 除以 $y$ 的边缘）。
  代入 decomposability 的分解式：

  $ p_(-tau)(y, x) = (rho_0(x) p_tau^D (x, r(y)) p_"local" (r(y), y)) \/ rho_1(y) . $

  *第二步：把分母算开。* 按 $rho_1$ 的定义对 $x$ 积分：

  $ rho_1(y) = integral rho_0(x) p_tau^D (x, r(y)) p_"local" (r(y), y) d x
    = p_"local" (r(y), y) dot integral rho_0(x) p_tau^D (x, r(y)) d x , $

  第二个等号的依据：$p_"local" (r(y), y)$ 不含积分变量 $x$，
  提出积分号。这就是后文"分子里的 $p_"local"$ 与归一化里的它相消"中
  归一化里的那一份。

  *第三步：约分。* 把第二步代回第一步，$p_"local"$ 上下消去：

  $ p_(-tau)(y, x) = (rho_0(x) p_tau^D (x, r(y))) \/ (integral rho_0(x) p_tau^D (x, r(y)) d x) , $

  右端只通过 $r(y)$ 依赖 $y$，把它取作 $p_(-tau)^D (r(y), x)$ 即为所求
  （若分母还差一个只依赖 $r(y)$ 的常数因子，并进 $p_(-tau)^D$ 的定义即可）。

  充分性：由 $p_(-tau)^D (r(y), x)$ 反解
  $p_tau (x, y) = p_(-tau)(y, x) rho_1(y) \/ rho_0(x)
    = [ p_(-tau)^D (r(y), x) \/ integral rho_1(y) p_(-tau)^D (r(y), x) d y ] rho_1(y)$，
  其中分母用了 $integral rho_1(y) p_(-tau)^D (r(y), x) d y = rho_0(x)$
  （把 $p_(-tau)(y, x) = rho_0(x) p_tau (x, y) \/ rho_1(y)$ 乘 $rho_1(y)$
  后对 $y$ 积分，右端 $rho_0(x) integral p_tau (x, y) d y = rho_0(x)$，
  转移密度对 $y$ 归一）。方括号里第一项只依赖
  $(x, r(y))$，取它为 $p_tau^D (x, r(y))$；剩下的 $rho_1(y)$ 取为
  $p_"local"$（$p_"local"$ 的定义本就允许任意依赖目标点 $y$），
  正好是 $p_tau^D dot p_"local"$ 的形状。$square$
]

这条引理的用处是工程性的：它说明做反向拟合时，
只要把反向速度场也瓶颈到 $r(y)$ 上，就是在检验 decomposability。
下一节的损失函数正是这样构造的。

把两个准则翻译成算子语言，是论文给出的第二个等价：

#definition("降阶算子")[
  $hat(cal(K))_L h(x) = E[ h(X_(t+tau)) | r(X_t) = r(x) ]$，
  $hat(cal(T))_D h(x) = E[ h(X_t) | r(X_(t+tau)) = r(x) ]$。
]

于是“无损”有第三种说法：$hat(cal(K))_L = cal(K)_tau$ 等价于 lumpability，
$hat(cal(T))_D = cal(T)_tau$ 等价于 decomposability。
这三种说法（密度等式、反向密度等式、算子等式）在论文里互为引理，
它们把“评价一个反应坐标”这件事变成“比较两个算子”。

== 最优性定理：损失值就是坐标质量

flow matching 的部分按 rectified flow 写。虚拟时间 $s in [0,1]$，
插值路径 $Y^s = (1 - s) Y^0 + s Y^1$，
最优速度场是条件期望 $v^*(s, y) = E[ Y^1 - Y^0 | Y^s = y ]$，
生成靠解 ODE $d Y^s \/ d s = v(s, Y^s)$，训练损失是
$integral_0^1 E[ norm(v(s, Y^s) - (Y^1 - Y^0))^2 ] d s$。

论文同时拟合正向与反向两个条件密度。记
$X tilde rho_0$、$X_1 tilde p_tau (X, dot)$、$X' tilde cal(N)(0, I)$（用于反向传播的辅助高斯），
两个全速度场 $v_0, v_1 : RR times RR^D times RR^D -> RR^D$ 的最小化目标是
$cal(L) = cal(L)_0 + cal(L)_1$，其中

$ cal(L)_0[v] = E[ norm(v(s, X_1^s, X) - (X_1 - X'))^2 ] , $
$ cal(L)_1[v] = E[ norm(v(s, X^s, X_1) - (X - X'))^2 ] , $

第一项把第三个参数当“过去”，拟合 $p_tau (dot | X)$；
第二项把第三个参数当“未来”，拟合反向密度。

瓶颈就是在第三个参数上做的：把 $X$ 换成 $r(X)$、$X_1$ 换成 $r(X_1)$，得到
$cal(L)^"FMRC"[v_0, v_1, r] = cal(L)_0[v_0^"RC" (dot, dot, r(dot))]
  + cal(L)_1[v_1^"RC" (dot, dot, r(dot))]$，
然后对 $(r, v_0, v_1)$ 共享参数 $theta$ 联合最小化。

#proposition("定理（最优性）")[
  对给定的反应坐标 $r$，
  $ min_(v_0^"RC", v_1^"RC") cal(L)^"FMRC" >= min_(v_0, v_1) cal(L) $，
  且等号成立当且仅当 $p_tau$ 关于 $r$ 同时 lumpable 与 decomposable。
]

#proof[
  左边的可行集是右边可行集的子集（任何只读 $r$ 的速度场也是全速度场），
  所以不等式平凡。

  等号 $=>$ 两准则成立：取右端最优场 $v_0^*, v_1^*$ 与左端瓶颈最优场
  $v_0^("RC",*), v_1^("RC",*)$。等号意味着两者几乎处处相等，
  于是正向流 $d Y^s \/ d s = v_0^*(s, dot, X)$ 的起点只通过 $r(X)$ 进入方程，
  它诱导的密度可写成 $p_tau^L (r(x), dot)$，这正是 lumpability；
  反向流同理给出 $p_(-tau)^D (r(y), dot)$，按上一条引理换成 decomposability。

  两准则成立 $=>$ 等号：lumpable 时
  $E[ X_1 - X' | X_1^s = y, X = x ]$ 在固定 $r(x)$ 的水平集上是常数，
  所以 $v_0^*(s, y, x) = v_0^*(s, y, r(x))$，最优场本身就已经瓶颈化；
  decomposable 时对反向场做同样的论证。两者合起来，
  瓶颈问题的最优值不小于全问题的最优值（平凡方向）又不大于它（取到同一个场），
  故相等。$square$
]

这个定理的实际含义是把“验证一个反应坐标”变成一次前向计算：
在固定的 $r$ 上训到底，$cal(L)^"FMRC"$ 的收敛值与全模型的收敛值之差
就是这个坐标的代价。它不需要真样本，也不需要预先知道亚稳态划分。

== 误差链条的四块砖

最优性是定性的（等号何时成立），定量版本要说“差一点等于损失值多大”。
链条有四块，每块都是一个引理。

*砖一：$W_2$ 与齐次 Sobolev 范数的比较。*
在 Riemannian 流形 $M$ 上
$W_2(mu, nu) = inf_(gamma in Pi(mu, nu)) [ integral d(x, y)^2 d gamma ]^(1/2)$，
半范数 $norm(f)_(dot(H)^1(m)) = ( integral |nabla f|^2 d m )^(1/2)$，
其对偶定义 $norm(nu)_(dot(H)^(-1)(m)) = sup { |angle(f, nu)| : norm(g)_(dot(H)^1(m)) <= 1 }$。
Benamou–Brenier 线性化说 $W_2(mu, mu + d mu) = norm (d mu) |_dot (H^(-1))(mu) + o(d mu)$。
需要的引理（Loeper 型）：若 $mu, nu <= beta dot lambda$（两个测度的密度都被
$beta$ 倍体积测度控制），则

$ norm(mu - nu)_(dot(H)^(-1)) <= sqrt(beta) W_2(mu, nu) . $

方向要说清：是 $W_2$ 控制 $dot(H)^(-1)$ 范数，不是反过来。
直觉是 $dot(H)^(-1)$ 是“最粗的”输运度量之一，
把质量搬动一点点的代价在它下面被低估，而 $W_2$ 按距离计价。

*砖二：速度场误差到分布误差。*
若
$integral_0^1 E[ norm(v_0^("RC",*) (s, X_1^s, X) - v_(0, theta)^("RC")(s, X_1^s, X))^2 ] d s <= epsilon^2$
且 $v_(0, theta)^("RC")$ 对 $(y^s, x)$ 是 $L_s$-Lipschitz，则
联合律 $(X, X_1)$ 与 $(X, hat(X)_1)$ 之间的 $W_2$ 不超过
$epsilon dot e^(integral_0^1 L_s d s)$。
证明的做法是把速度场补一个零块，让 $X$ 分量不动，
然后套用 Benton–Deligiannidis–Doucet 的 flow matching 稳定性结果
（Gronwall 型：路径误差随速度误差指数放大）。

*砖三：镜像的一支。* 对反向那半有同样的结论，误差记 $tilde(epsilon)$、
Lipschitz 常数记 $tilde(L)_s$，控制 $(hat(X), X_1)$ 与 $(X, X_1)$ 的 $W_2$。
论文对这一支只写了一句“证明同前，略”。

*砖四：密度误差到算子误差。*
弱算子范数定义为

$ norm(cal(K)_tau - hat(cal(K))_L)_H
    = sup { angle(g, cal(K)_tau h - hat(cal(K))_L h)_(rho_0) :
        norm(g)_(H^1(rho_0)) <= 1, norm(h)_(H^1(rho_1)) <= 1} , $

把联合密度的 $dot(H)^(-1)$ 误差限制在形如 $g(x) h(y)$ 的张量测试函数上，
正好就是上式。这个"正好"拆成三步：

+ *把算子作用写成对联合密度的配对。* Koopman 算子按密度前传：
  $(cal(K)_tau h)(x) = integral p_tau (x, y) h(y) d y$。代入内积的定义：

  $ angle(g, cal(K)_tau h)_(rho_0)
    = integral g(x) [ integral p_tau (x, y) h(y) d y ] rho_0(x) d x
    = integral integral g(x) h(y) rho_0(x) p_tau (x, y) d y d x , $

  读法：末端的被积对象是联合密度 $rho(x, y) = rho_0(x) p_tau (x, y)$
  与测试函数 $g(x) h(y)$ 的乘积。
+ *两个算子相减，认出密度差。* $hat(cal(K))_L$ 同样按某个联合密度前传
  （水平集内的平均，对应密度 $hat(rho)$）。两个配对相减时
  测试函数 $g h$ 不变，剩下的全部差别在被积的密度里：
  $angle(g, (cal(K)_tau - hat(cal(K))_L) h)_(rho_0)
    = angle(g h, rho - hat(rho))$，
  右端是联合密度差 $delta rho$ 对测试函数 $g h$ 的作用。
+ *部分上确界不超过全上确界。* 算子范数的上确界只遍历
  能写成 $g(x) h(y)$ 乘积形状的测试函数；$dot(H)^(-1)$ 范数的上确界
  遍历全部单位球测试函数，张量积是其中的子集。
  子集上的最大读数不超过全集上的最大读数，
  于是 $norm(cal(K)_tau - hat(cal(K))_L)_H <= norm(rho - hat(rho))_(dot(H)^(-1))$。

于是主定理：

#proposition("定理（主定理）")[
  在砖二、砖三成立且砖一的密度支配条件以常数 $beta, tilde(beta)$ 成立时，

  $ norm(cal(K)_tau - hat(cal(K))_L)_H <= sqrt(beta) dot epsilon dot e^(integral_0^1 L_s d s) , $

  $ norm(cal(T)_tau - hat(cal(T))_D)_H <= sqrt(tilde(beta)) dot tilde(epsilon) dot e^(integral_0^1 tilde(L)_s d s) . $
]

（源码里这一条的第二式在定理陈述处把 $hat(cal(T))_D$ 写成了 $hat(cal(T))_L$，
附录证明里用的是 $D$。这是笔误，不影响结论。）

#warn("这条链上真正强的与真正弱的部分")[
  强的是形状：损失值（砖二的左端）确实直接进了算子误差的右端，
  中间没有“存在某常数使得”之类的非构造环节。
  弱的是常数：$e^(integral L_s)$ 里的 Lipschitz 常数在实践中拿不到上界；
  密度支配 $beta$ 要求联合密度不出现极低的沟谷，
  多峰系统里峰间的低密度区正是要跨的东西。
  所以这条链给的是“误差随损失消失”的定性保证，
  不能拿来给具体数值定界。
]

论文里还有第五块砖被整段注释掉了：平衡情形（$rho_0 = rho_1 = mu$）下的
$L^2_(mu)$ 误差界。它的思路是用连续性方程展开 $delta f = f - hat(f)$，
把 $norm delta v <= epsilon_1^2$ 与 $norm nabla dot delta v <= epsilon_2^2$
代进 $u' <= alpha u + beta$ 型微分不等式，再走 Gronwall，
目标是 $norm(cal(K)_tau - hat(cal(K))_L)_op <= norm dot _"HS"$，
接上谱误差（引 Bittracher 等的引理 3.5）。
附录里这段推导停在一个未完成的不等式
$integral_(Omega_Y) mu nabla dot (v_X delta f) d y <=$ 上。
正稿的取舍是：保留 $H$ 范数那条完整链条，
把没写完的谱误差留在注释里。

== 算法与它没说的事

训练循环六步，写得很具体：取小批量 ${ (x_b, y_b) }_(b=1)^B$；
抽 $x'_b, y'_b ~ cal(N)(0, I)$ 与 $s ~ cal(U)(0, 1)$；
构造插值 $x_b^s = I(s, x'_b, x_b)$、$y_b^s = I(s, y'_b, y_b)$；
算

$ hat(L) = 1/B sum_b norm(v_0^"RC" (s, y_b^s, r_theta (x_b)) - (y_b - y'_b))^2
        + 1/B sum_b norm(v_1^"RC" (s, x_b^s, r_theta (y_b)) - (x_b - x'_b))^2 $

；随机梯度步 $theta <- theta - eta dot partial hat(L) \/ partial theta$；重复到收敛。

有三处“没说”要如实记下。
一，网络只写成多层感知机，没有层数、宽度、批量、迭代次数、学习率、求解器容差。
二，推断阶段要生成样本就得解那两个 flow ODE，论文把它列为可行但没有实验。
三，源码里瓶颈速度场的类型写成
$v_0^"RC", v_1^"RC" : RR times RR^D times RR^d -> RR^d$，
值域应当是 $RR^D$（要输出下一状态），是笔误。
另外损失里没有散度项，也没有 Hutchinson 估计器；
散度只出现在被注释掉的 $L^2$ 推导里。

== 唯一的实验说明了什么，没说明什么

数据构造：$RR^3$ 上的过阻尼 Langevin
$d X_t = -nabla V(X_t) d t + sqrt(2 beta^(-1)) d W_t$，势能

$ V(x_1, x_2, x_3) = V'(x_1, x_2) + 10 x_3^2 , $
$ V'(x_1, x_2) = cos(7 dot arctan(x_2, x_1)) + 10 ( sqrt(x_1^2 + x_2^2) - 1 )^2 , $

也就是圆周上七个井，外加一个把 $x_3$ 牢牢按在 0 附近的强约束，
其平衡分布是 $cal(N)(0, 0.05)$，混合时间远小于 $(x_1, x_2)$ 方向。
用 Euler–Maruyama 模拟，步长 $Delta t = 0.001$，
再套一个把平面卷成三维流形的 Swiss roll 映射，
使观测空间中的线性降维注定失败。
反应坐标维度取 $d = 1$。

验证方式是把学到的 $r(x)$ 投回二维平面看颜色分层，
再用 PCCA+ 在投影数据上做亚稳态聚类，检查两件事：
每个 PCCA+ 簇内的 RC 值分布是否窄，簇之间是否有间隙。

#figure(
  align(center, grid(columns: 2, gutter: 8pt,
    image("../fig_c/x_c4_pot.png"), image("../fig_c/x_c4_traj.png"))),
  caption: [数据构造（论文 2408.15981 图 1 的两格）：
  左为圆周上七个井的势能面，右为加上 $x_3$ 方向强约束后的三维采样点。
  真实动力学本身是二维的，三维只是被刻意“搅浑”的观测。],
)

#figure(
  align(center, grid(columns: 3, gutter: 6pt,
    image("../fig_c/x_c4_cluster.png"), image("../fig_c/c4_rc.png"),
    image("../fig_c/x_c4_hist.png"))),
  caption: [三步验证（论文 2408.15981）：左为平面投影上的聚类真值；
  中为学到的 $r(x)$ 着色在同一投影上，七口井各自连续成段；
  右为 PCCA+ 各簇内 RC 值的直方图，簇内窄、簇间有间隙。
  这正是 lumpability 的可视图证：同一簇（近似同一水平集）里的点
  有相似的未来，不同簇被 $r$ 分开。],
)

这篇没有的东西同样要说清：没有数值表格，没有与 TICA、VAMPnets、
time-lagged 自编码器的定量对照，没有 alanine dipeptide 之类的真实分子数据，
没有训练超参数，引用的补充材料文件不在源码包里。
它给出的是一件结构性的事（损失与算子误差的可微联系）加一个演示。
所以合理的使用方式是：把 $cal(L)^"FMRC"$ 当成可优化的目标与可比较的坐标质量分数，
把它的双向瓶颈结构当成设计模板，
而不是把它的图示当成性能基准。

#figure(
  align(center, image("../fig_c/x_c4_meta.png", width: 82%)),
  caption: [PCCA+ 亚稳态与 RC 的并置图（论文 2408.15981 图 2b）：
  簇的边界与 RC 的间隙对齐。作者用的 PCCA+ 是事后聚类工具
  （Roblitz–Schütte 版本），不是模型的一部分。],
)

== 它站的位置

往回看，Bittracher–Mollenhauer–Koltai–Schütte（2023）给出了这两个准则
以及“满足准则即无损”的理论，但他们的数值流程要估计并积分转移密度，
在高维上跑不动。FMRC 的贡献不在提出准则，
而在把准则改写成*两个条件密度模型的拟合损失*，
于是可以用 flow matching 的现成工具箱端到端训练，
顺手把“优化时也在评价”这两件事合成了一个目标。
往前看，它留下的空格也很清楚：
多峰低密度沟谷下的密度支配、$L^2$ 与谱误差那一支、
以及任何意义上的定量基线比较。
