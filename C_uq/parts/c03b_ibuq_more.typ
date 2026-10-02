#import "../../preamble.typ": *

= IB-UQ 深入：率失真视角、算法细节与实验画廊

== 用率失真的语言重读 IB-UQ

第 5 章的推导可以再用一遍信息论的经典对照来读。
经典率失真问题：给定源 $X$ 与失真度量，求“失真不超过 $D$ 时
最少要多少码率”的曲线 $R(D)$。IB 的变体：失真不是平方误差，
而是“丢掉了多少关于 $Y$ 的信息”$I(X;Y|Z)$，码率是 $I(X;Z)$。
第 5 章第四步的改写 $cal(L)_"IB" = I(X;Z) + beta I(X;Y|Z) - beta I(X;Y)$
正是“码率 + $beta times$ 失真（差常数）”，$beta$ 是两条曲线的
兑换率。

IB-UQ 在这个框架上换了两样东西，对照如下：

#ttable(
  ([成分], [经典 IB], [IB-UQ]),
  (
    ([码率], [$I(X;Z)$], [$I(tilde(Z); tilde(X))$，在加温分布上算]),
    ([失真], [$I(X;Y|Z)$], [同（经变分下界近似）]),
    ([兑换率], [$beta$，扫参定], [$beta$，实验最优 $[0.2, 0.3]$]),
    ([无失真极限], [$Z$ 保留全部 $X$ 的信息], [$Z$ 只保留训练分布内的信息；OOD 处 $Z$ 与 $X$ 独立]),
  ),
  widths: (auto, 1fr, 1fr),
)

最后一行是关键差别：经典 IB 的“最优编码”会把 $X$ 的全部信息
（包括噪声）都录进 $Z$，遇到分布外的 $X$ 就瞎编；
IB-UQ 的最优编码被训练成“只在训练分布内依赖输入”，
OOD 输入自动落到先验。这就是“宽分布压缩”在率失真语言里的样子：
*把失真的定义域扩大到了训练分布的邻域*。

== Mixup 与 $alpha = 0.005$

算法里 Mixup 的系数 $alpha = 0.005$ 远小于原论文的 0.1 到 0.4，
值得解释。原 Mixup 用大 $alpha$ 做数据增强（标签平滑）；
这里 $lambda tilde "Beta"(alpha, alpha)$ 在 $alpha$ 很小时几乎总是
取 0 或 1（分布集中在两端），Mixup 退化为“偶尔把两个样本连线
中间的点加入训练”。它扮演的角色不是增强而是*平滑经验分布*：
数据少时 $q_D (y|z)$ 会过拟合到有限个点，偶尔的插值样本
把解码器的条件高斯拉宽。这是防 $q_D$ 过拟合的轻量手段，
与换更强的密度模型（normalizing flow 解码器，论文留作选项）
是两条互补的路。

== 算子版算法的完整流程

DeepONet 版（IB-ONet）每个训练步：

+ 抽 mini-batch 的输入函数-解对 $(o_b, s_b)$；
+ （可选）Mixup；
+ GIN 抽加温输入 $tilde(o)_b$；
+ 编码 $z_b = "diag"(m(o_b)) overline(z)(o_b) + "diag"(1 - m(o_b)) z_(0,b)$；
+ 对每个 $o_b$ 另抽 $B_M$ 个查询点 $y_j^b$（trunk 网络输入）；
+ 估计
  $hat(cal(L)) = 1/(B B_M) sum_(b, j) log q_D (s(y_j^b) | y_j^b, z_b) - beta/B sum_b [log q_E (tilde(z)_b | tilde(o)_b) - log e(tilde(z)_b)]$；
+ 梯度上升更新 $q_D, q_E, e$（全部网络）。

与回归版唯一的结构差异在第 5 步：解码器的输入是
(查询点，潜变量) 二元组，所以每个样本要配多个查询点才能
约束整个函数。$B_M$ 的选择是精度与代价的折中。

== 实验画廊

#figure(
  align(center, image("../fig_c/tmp_c2_four.png", width: 88%)),
  caption: [不连续函数的四种方法对比（论文 2302.03271）：
  GP、HMC、深度集成、IB-UQ。预测均值（线）与 2 倍标准差带（阴影）。
  GP/HMC/IB-UQ 在外推区（两端与中间断点附近）不确定度都增大；
  深度集成的带在外推区没有张开，印证综述的“集成盲区”。],
)

#figure(
  align(center, image("../fig_c/tmp_c2_don.png", width: 84%)),
  caption: [算子版结构（论文 2302.03271）：瓶颈加在 DeepONet 的
  branch 网络（传感器测量一侧），置信门 $m(o)$ 决定每个维度
  走确定性特征还是先验噪声；trunk 网络保持确定性。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_c/x_c2_1d.png"), image("../fig_c/x_c2_2d.png"))),
  caption: [反应扩散与对流方程的算子学习 UQ（论文 2302.03271）。
  左：一维反应扩散在 OOD 输入（不同相关长度）下的预测与不确定带；
  右：二维情形的代表样本。误差大的区域与不确定度大的区域重合。],
)

#figure(
  align(center, image("../fig_c/x_c2_climate.png", width: 88%)),
  caption: [气候模型实验（温度场到气压场）：删除部分季度数据做 OOD。
  上：完整数据；中、下：分别删除第二季度与二、三季度。
  IB-UQ 的不确定度在数据缺失的区域增大，而普通 DeepONet
  不输出任何置信度信息，误差无人报告。],
)
