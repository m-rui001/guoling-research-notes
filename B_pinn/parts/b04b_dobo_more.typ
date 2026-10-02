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

高斯随机变量的特征函数给出 $E[cos xi t] = e^(-sigma^2 t^2 \/ 2)$、
$E[sin xi t] = 0$（对称分布的一阶三角矩为零），于是

$ overline(u)(x, t) = -sin x dot e^(-sigma^2 t^2 \/ 2) . $

均值按 $e^(-sigma^2 t^2 \/ 2)$ 衰减：随机平流把初始波形逐渐
“搅散”进波动项里。物理直觉：不同样本的波以不同速度移动，
位置对齐的样本越来越少，逐点平均自然缩小，但能量没丢
（见下面的振幅）。

== 分解随机部分

$ u - overline(u) = -sin x (cos xi t - e^(-sigma^2 t^2 \/ 2)) + cos x dot sin xi t . $

两项的随机系数都零均值，且恰好由两个*互相正交的空间形状*承载：
$cos x$ 项与 $sin x$ 项（在 $[-pi, pi]$ 上 $integral sin x cos x d x = 0$）。
这提示二模态截断是精确的：随机结构张成的空间就是
$"span"{cos x, sin x}$。

== 组装 DO 形式

取归一化空间模态 $u_1 = -cos x \/ sqrt(pi)$、$u_2 = -sin x \/ sqrt(pi)$
（各自 $⟨u_i, u_i⟩ = integral_(-pi)^pi u_i^2 d x = 1$，且互正交）。
投影随机部分：

$ ⟨u - overline(u), u_1⟩ = sqrt(pi) dot sin xi t , quad ⟨u - overline(u), u_2⟩ = sqrt(pi)(cos xi t - e^(-sigma^2 t^2 \/ 2)) . $

DO 约定 $u = overline(u) + a_i u_i Y_i$ 且 $E[Y_i^2] = 1$，于是振幅取
随机系数的标准差、随机模态取标准化后的系数：

$ a_1(t)^2 = pi dot E[sin^2 xi t] = pi (1 - e^(-2 sigma^2 t^2))\/2 , quad Y_1(t) = sin xi t \/ a_1(t) , $

$ a_2(t)^2 = pi dot E[(cos xi t - e^(-sigma^2 t^2 \/ 2))^2] = pi (1 - e^(-sigma^2 t^2))^2 \/ 2 , quad Y_2(t) = (cos xi t - e^(-sigma^2 t^2 \/ 2)) \/ a_2(t) . $

（$E[cos xi t] = e^(-sigma^2 t^2)$ 的平方配 $E[cos^2] = (1 + e^(-2 sigma^2 t^2))\/2$
给出第二个方差。）这套分解在任意时刻都精确，且：
*振幅从零开始增长*（$a_1(0) = a_2(0) = 0$：确定性初值没有随机能量，
随机性由随机平流速度逐步注入）；
*空间模态恒定*（这个例子里动力学不产生新的空间方向，
所以 DO 条件 $⟨partial_t u_i, u_j⟩ = 0$ 自动满足）；
*经典 DO 方法在这个例子上能工作*：协方差矩阵
$"diag"(a_1^2, a_2^2)$ 只在 $t = 0$ 奇异，而数值积分从 $t > 0$ 起
就可以推进（或者像 NN-DO 一样从 $t = 0$ 用合法的零振幅初始化）。

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
两者对解的表示完全相同。直觉：两类约束只是“把冗余摊到哪里”
的不同约定：DO 让空间基的演化方向垂直于既有张成空间
（动态约定），BO 让两组基在任意时刻互相正交（静态约定）；
它们消除的是同一个自由度。NN 版本把两种约定都写成损失后，
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
