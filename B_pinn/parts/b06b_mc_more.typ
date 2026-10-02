#import "../../preamble.typ": *

= 分数阶算子的三笔推导

主章节给出了分数阶 Laplacian 与 Caputo 导数的 MC 估计公式。
这里把公式怎么来的走完：三步对称化、内外分段、Beta 采样密度的来历。
另外补一段算子本身的性质，解释它在物理里对应什么。

== 分数阶 Laplacian 的性质速览

#definition("分数阶 Laplacian 的三个基本事实")[
  *傅里叶符号*：$cal(F)[(-Delta)^(alpha\/2) u](omega) = norm(omega)^alpha cal(F)[u](omega)$。
  即它在频率空间乘 $norm(omega)^alpha$：$alpha = 2$ 时是 $-Delta$
  （乘 $norm(omega)^2$），$alpha$ 越小高频被放大的越少、
  算子越“弥散”。

  *$alpha = 2$ 的退化*：常数 $C_(d, 2)$ 的分母 $|Gamma(-1)| -> infinity$，
  主值积分收敛到经典 Laplacian。$alpha -> 0$ 时退化为恒等算子。

  *概率对应*：它是稳定 Lévy 过程的生成元，$alpha$ 越小，
  对应的随机游走越“跳跃”（重尾步长分布）。
  分数阶方程描述的异常扩散，物理上就是粒子做长程跳跃。
]

这三个事实解释了第 10 章开头的两个坏性质从哪来：核
$norm(x - y)^(-d - alpha)$ 在大距离衰减慢（重尾跳跃），
在近距离奇异（生成元的积分表示）。

== 内段公式的推导

目标：把主值积分的球内部分写成期望。固定 $x$，用球坐标
$y = x + r xi$（$xi$ 均匀于单位球面 $S^(d-1)$，$d y = r^(d-1) d r d xi$）：

$ integral_(B_(r_0)(x)) (u(x) - u(y))/(norm(x - y)^(d + alpha)) d y = integral_0^(r_0) r^(d - 1 - d - alpha) integral_(S^(d-1)) (u(x) - u(x + r xi)) d xi d r . $

球面平均的对称化是关键一步：对 $xi$ 与 $-xi$ 求平均
（球面均匀分布对 $-xi$ 不变），

$ integral_(S^(d-1)) (u(x) - u(x + r xi)) d xi = 1/2 integral_(S^(d-1)) [2 u(x) - u(x + r xi) - u(x - r xi)] d xi , $

方括号里的二阶差商在 $r -> 0$ 时是 $O(r^2)$（一阶项相消），
除以 $r^2$ 后有界，这就是把“主值发散”变成“可积期望”的全部秘密。
代回并整理 $r$ 的幂：

$ = (|S^(d-1)|)/(2) integral_0^(r_0) (2 u(x) - u(x + r xi) - u(x - r xi)) r^(-1 - alpha) d r dot (r_0^(2 - alpha))/(r_0^(2 - alpha)) , $

把 $r^(-1-alpha)$ 吸进采样密度：令 $f_I (r) ∝ r^(1 - alpha)$
（在 $[0, r_0]$ 上归一化），则

$ integral_0^(r_0) g(r) r^(-1 - alpha) d r = (r_0^(2 - alpha))/(2 - alpha) E_(r tilde f_I)[ g(r) \/ r^2 ] , quad f_I (r) = (2 - alpha)/(r_0^(2 - alpha)) r^(1 - alpha) . $

合并常数即得主章节的内段公式（$|S^(d-1)| r_0^(2 - alpha) \/ (2(2 - alpha))$
乘期望）。*采样密度正比于积分核的径向部分*：被估量 $g(r)\/r^2$
在 $f_I$ 下有界（$g(r) = O(r^2)$ 恰好抵消 $r^2$ 分母），
方差有限。Beta 分布采样 $r\/r_0 tilde "Beta"(2 - alpha, 1)$
给出这条幂律密度，一次 Beta 样本换一次网络评估。

== 外段与 Caputo 的推导思路

外段（$r in [r_0, infinity)$）同样对称化，差异只在权重常数：
$r^(-1 - alpha)$ 在 $[r_0, infinity)$ 的积分是 $r_0^(-alpha)\/alpha$，
采样密度 $f_O (r) = alpha r_0^alpha r^(-1 - alpha)$（重尾：
远处贡献小但覆盖全域），对应 $r_0\/r tilde "Beta"(alpha, 1)$。
被估量不再除 $r^2$（远处二阶差商本身 $O(r^2)$ 衰减，
不需要除法去奇异），所以外段期望里的量是
$2 u(x) - u(x + r xi) - u(x - r xi)$ 本身。

Caputo 导数的 MC 化同理：把
$integral_0^t (t - tau)^(-gamma) partial_tau u d tau$ 分部积分
（$u(t) - u(0)$ 项析出），剩余积分的核 $tau^(-gamma)$
用 $f_(I, t)(tau) = (1 - gamma) tau^(-gamma)$ 吸收，
得到主章节的公式。三处推导共用一个模式：
*找出积分核的径向/时间部分，令其归一化后作为采样密度，
被估量换成就地有界的差商*。这是重要性采样在算子估计里的标准打法。

== MC-PINN 实验画廊

#figure(
  align(center, grid(columns: 3, gutter: 8pt,
    image("../fig_b/x_b6_fwd2d.png"), image("../fig_b/x_b6_f2r.png"), image("../fig_b/x_b6_ade3d.png"))),
  caption: [左：2 维分数阶 Laplacian 方程的精确解、预测与绝对误差
  （论文 2203.08501）：人工解 $(1 - norm(x)^2)^(1 + alpha\/2)$，
  误差集中在边界附近。右：时空分数阶 ADE 反问题在 3 维的参数
  收敛（$alpha, gamma, c, v$ 四参数同时识别）。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_b/x_b6_rand5d.png"), image("../fig_b/x_b6_dmu.png"))),
  caption: [左：随机参数分数阶扩散方程（$d = 5$）的代理预测与误差。
  右：ABC 筛选得到的 $mu$ 后验分布：真值落在后验峰附近，
  验证代理 + ABC 组合可以做完整的不确定性量化。],
)
