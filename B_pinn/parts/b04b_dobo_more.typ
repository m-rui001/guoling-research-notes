#import "../../preamble.typ": *

= 手推一个模态分解：随机平流方程的 DO 结构

主章节引用了“平流方程有精确 DO 解”。这一节把它完整推出来，
推完之后 DO 的四个组件（均值、空间模态、振幅、随机模态）
各自“长什么样、为什么长这样”就都有了具体的脸。

== 方程与解

$ u_t + xi u_x = 0, quad u(x, 0) = -sin x, quad xi tilde cal(N)(0, sigma^2), quad x in [-pi, pi], $

$xi$ 是随机平流速度。特征线解法：$u(x, t) = -sin(x - xi t)$。
这是整个推导的锚点：任何声称的模态分解都可以与它逐项核对。

== 均值

$ overline(u)(x, t) = E[-sin(x - xi t)] = -sin x dot E[cos xi t] + cos x dot E[sin xi t] . $

这个等式链有两步，拆开看。

*第一步（拆开角度）*：和角公式 $sin(x - xi t) = sin x cos xi t - cos x sin xi t$
对每个固定样本值逐点成立，所以对随机变量也成立。

*第二步（期望的线性）*：$sin x$ 与 $cos x$ 不依赖 $omega$，是常数因子，
移出期望即得上式。剩下的问题是算两个三角矩。设
$I(s) := E[cos xi s]$，对 $s$ 求导并在积分号下求导（$xi$ 的密度
衰减很快，交换合法），再用分部积分。记 $xi$ 的密度为 $phi$，
它满足 $x phi(x) = -sigma^2 phi'(x)$（直接对密度表达式求导可验证）：

$ I'(s) = -E[xi sin xi s] = -integral sin(s x) x phi(x) d x = sigma^2 integral sin(s x) phi'(x) d x = -sigma^2 s I(s) , $

第三处用 $x phi = -sigma^2 phi'$，第四处分部积分（边界项为零，
$phi$ 在无穷远处衰减到零），得到常微分方程 $I'(s) = -sigma^2 s I(s)$，
初值 $I(0) = 1$，解为 $I(s) = e^(-sigma^2 s^2 \/ 2)$（这就是高斯分布
特征函数的模长，概率论的标准结果）。对称性给出另一个矩：
$E[sin xi s] = 0$，因为被积函数 $sin(s x) phi(x)$ 是奇函数，
在对称区间上积分为零。取 $s = t$，得

$ overline(u)(x, t) = -sin x dot e^(-sigma^2 t^2 \/ 2) . $

均值按 $e^(-sigma^2 t^2 \/ 2)$ 衰减：随机平流把初始波形逐渐
“搅散”进波动项里。物理直觉：不同样本的波以不同速度移动，
位置对齐的样本越来越少，逐点平均自然缩小，但能量没丢
（见下面的振幅）。

== 分解随机部分

把上一步算出的 $overline(u)$ 从 $u = -sin(x - xi t)$ 里减去，
再用一次和角公式：

$ u - overline(u) = -sin x (cos xi t - e^(-sigma^2 t^2 \/ 2)) + cos x dot sin xi t . $

两项的随机系数都零均值，且恰好由两个*互相正交的空间形状*承载：
$cos x$ 项与 $sin x$ 项（在 $[-pi, pi]$ 上 $integral sin x cos x d x = 0$，
被积函数是奇函数）。这提示二模态截断是精确的：随机结构张成的空间就是
$"span"{cos x, sin x}$。

== 组装 DO 形式

取归一化空间模态 $u_1 = -cos x \/ sqrt(pi)$、$u_2 = -sin x \/ sqrt(pi)$。
归一化与正交各验一遍：$cos^2 x = (1 + cos 2 x) \/ 2$ 在整周期上积分时
$cos 2 x$ 项积分为零，所以 $integral_(-pi)^pi cos^2 x d x = pi$，
$⟨u_1, u_1⟩ = 1$，$u_2$ 同理；$⟨u_1, u_2⟩$ 的被积函数 $sin x cos x$
是奇函数，对称区间积分为零。

*第一步（投影到 $u_1$）*：记 $m := e^(-sigma^2 t^2 \/ 2)$（均值里的衰减因子，
对 $omega$ 是常数）。把 $u - overline(u)$ 与 $u_1$ 的内积逐项积分：

$ ⟨u - overline(u), u_1⟩ = 1/sqrt(pi) integral_(-pi)^pi [-sin x (cos xi t - m) + cos x sin xi t] (-cos x) d x = (cos xi t - m)\/sqrt(pi) integral sin x cos x d x - sin xi t \/sqrt(pi) integral cos^2 x d x . $

第一处积分为零（奇函数），第二处等于 $pi$，于是

$ ⟨u - overline(u), u_1⟩ = -sqrt(pi) sin xi t . $

$x$ 的积分把空间形状挑出来，
随机系数 $cos xi t - m$、$sin xi t$ 原样留在系数位置。

*第二步（投影到 $u_2$）*：同样的算法，这次 $sin^2 x$ 的积分留下、
$sin x cos x$ 的积分为零：

$ ⟨u - overline(u), u_2⟩ = (cos xi t - m)\/sqrt(pi) integral sin^2 x d x = sqrt(pi)(cos xi t - m) . $

因为 $u_1, u_2$ 标准正交，投影系数就是展开系数：
$u - overline(u) = c_1 u_1 + c_2 u_2$，$c_1 = -sqrt(pi) sin xi t$，
$c_2 = sqrt(pi)(cos xi t - m)$。

*第三步（拆成振幅乘标准随机模态）*：DO 约定
$u = overline(u) + sum_i a_i u_i Y_i$ 且 $E[Y_i^2] = 1$。比较得
$c_i = a_i Y_i$，于是振幅取系数的标准差、随机模态取标准化后的系数：
$a_i^2 = E[c_i^2]$，$Y_i = c_i \/ a_i$。

*第四步（算两个方差）*。第一个：
$E[sin^2 xi t] = (1 - E[cos 2 xi t]) \/ 2$（倍角公式 $sin^2 theta = (1 - cos 2 theta) \/ 2$），
三角矩公式（上一节，取 $s = 2 t$）给 $E[cos 2 xi t] = e^(-2 sigma^2 t^2)$，所以

$ a_1(t)^2 = pi dot E[sin^2 xi t] = pi (1 - e^(-2 sigma^2 t^2))\/2 . $

第二个用平方展开加期望的线性：

$ E[(cos xi t - m)^2] = E[cos^2 xi t] - 2 m E[cos xi t] + m^2 . $

三块各自有：$E[cos^2 xi t] = (1 + E[cos 2 xi t]) \/ 2 = (1 + e^(-2 sigma^2 t^2)) \/ 2$
（倍角公式 $cos^2 theta = (1 + cos 2 theta) \/ 2$ 加三角矩公式）；
$E[cos xi t] = m$（三角矩公式取 $s = t$）；$m^2 = e^(-sigma^2 t^2)$。
合并通分：

$ (1 + e^(-2 sigma^2 t^2))\/2 - 2 e^(-sigma^2 t^2) + e^(-sigma^2 t^2) = (1 - 2 e^(-sigma^2 t^2) + e^(-2 sigma^2 t^2))\/2 = (1 - e^(-sigma^2 t^2))^2\/2 , $

所以

$ a_2(t)^2 = pi dot E[(cos xi t - e^(-sigma^2 t^2 \/ 2))^2] = pi (1 - e^(-sigma^2 t^2))^2 \/ 2 . $

随机模态按第三步的公式写出，并验证归一化
（$E[Y_1^2] = pi E[sin^2 xi t] \/ a_1^2 = 1$，$Y_2$ 同理）：

$ Y_1(t) = -sqrt(pi) sin xi t \/ a_1(t) , quad Y_2(t) = sqrt(pi)(cos xi t - e^(-sigma^2 t^2 \/ 2)) \/ a_2(t) . $

系数 $c_1 = -sqrt(pi) sin xi t$ 就是主章节图注引用的精确解
$-sqrt(pi) sin (xi t)$ 的来源（那里把乘积 $a_1 Y_1$ 整体写作随机模态）。
这套分解在任意时刻都精确，且：
*振幅从零开始增长*（$a_1(0) = a_2(0) = 0$：确定性初值没有随机能量，
随机性由随机平流速度逐步注入）；
*空间模态恒定*（这个例子里动力学不产生新的空间方向，
所以 DO 条件 $⟨partial_t u_i, u_j⟩ = 0$ 自动满足）；
*经典 DO 方法在这个例子上能工作*：协方差矩阵
$"diag"(a_1^2, a_2^2)$ 只在 $t = 0$ 奇异，而数值积分从 $t > 0$ 起
就可以推进（或者像 NN-DO 一样从 $t = 0$ 用合法的零振幅初始化）。

#insight("用日常语言读一遍")[
  这几步合起来做了一次完整的分账。先算均值：对随机速度取期望，
  波形按 $e^(-sigma^2 t^2 \/ 2)$ 缩小。再把均值之外的部分投影到
  $cos x$ 与 $sin x$ 两个固定形状上：涨落全部落在这两个方向，
  没有第三个方向参与。最后把每个投影系数拆成“尺度乘标准随机量”：
  尺度是振幅 $a_i$，标准随机量是 $Y_i$。拆完之后，DO 展开的四个成分
  各有一个闭式，网络训练出的每个输出都可以与公式对照。
]

== 这个例子检验了什么

NN-DO 在此问题上的全部分量与上面手推的公式对照，
相对 $L^2$ 误差 2% 以内（主章节的图）。它的价值是“可核对性”：
每一个网络输出都有一个闭式真值。更复杂的实验
（Burgers 特征值交叉、19 维扩散反应）没有这种真值，
只能对照蒙特卡洛；但平流例子的核对保证了“网络学的分解结构
与理论定义一致”，这一致性是后面所有实验可信的前提。

== DO 与 BO 的等价性

Choi 等（2014）证明了 DO 与 BO 经微分变换相互等价：
给定任一 DO 分解，可以构造对应的 BO 分解（反之亦然），
两者对解的表示完全相同。等价性的根据是冗余的具体形态，
先在最简单的单模态情形把它看清楚。$N = 1$ 时展开是
$u = overline(u) + a_1 u_1 Y_1$，把 $u_1$ 换成 $2 u_1$、$a_1$ 换成 $a_1 \/ 2$、
$Y_1$ 换成 $Y_1$，乘积 $a_1 u_1 Y_1$ 与 $u$ 都不变，但空间模态的范数
翻了倍。多模态时同样的自由度变成一族随时间变化的旋转与重标度：
解 $u$ 不变，分量怎么分账可以变。DO 与 BO 是把这份自由度
钉死的两种约定：DO 用微分条件（$⟨partial_t u_i, u_j⟩ = 0$，
空间基只在新增方向上演化），BO 用代数条件（$⟨u_i, u_j⟩ = lambda_i delta_(i j)$
与 $E[Y_i Y_j] = delta_(i j)$，两组基时刻正交）。Choi 等的结果是：
沿时间构造那个保乘积的重标度变换，可以把 DO 分解的分量
逐步重整成满足 BO 条件的分量，且变换可逆，所以两个约束集合
约束的是同一批分解里的不同代表元，表示能力没有差别。
NN 版本把两种约定都写成损失后，
这个等价性表现为：NN-DO 与 NN-BO 在 Burgers 实验里的精度相同
（主章节表格），选择哪一个只是看哪种约束的损失更容易优化。

== DO/BO 实验画廊

#figure(
  align(center, image("../fig_b/x_b2_advmean.png", width: 78%)),
  caption: [平流方程的均值场：NN-DO 的 $overline(u)_(n n)$（点）与闭式
  $-sin x dot e^(-sigma^2 t^2 \/ 2)$（线）对照，衰减形状一致。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_b/x_b2_noisy.png"), image("../fig_b/x_b2_nlvar.png"))),
  caption: [左：噪声初值（30 个传感器读数加 $sigma = 0.1$ 高斯噪声）
  对末端预测几乎无影响：模态分解天然滤掉了高频观测噪声。
  右：19 维扩散反应方程的方差场：6 个 BO 模态捕捉大部分随机能量。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_b/x_b2_loss.png"), image("../fig_b/x_b2_invai.png"))),
  caption: [左：各损失分量的训练曲线：弱形式、约束、初值边界各自收敛，
  权重设计（约束项放大）的依据是它们的相对量级。
  右：反问题中各阶振幅 $a_i$ 的学习轨迹：反应系数识别
  （真值 0.5 与 0.3）在振幅结构上留下可读的印记。],
)
