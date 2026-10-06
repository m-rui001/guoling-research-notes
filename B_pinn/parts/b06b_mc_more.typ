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
$y = x + r xi$：$r$ 是 $y$ 到球心的距离，$xi$ 是单位向量，
均匀分布于单位球面 $S^(d-1)$。体积元分解为径向与角向两部分：
$d y = r^(d - 1) d r d xi$（半径方向的 $d r$ 乘球面元 $r^(d-1) d xi$）。
分母 $norm(x - y)^(d + alpha) = r^(d + alpha)$，与 $r^(d-1)$ 合并：

$ integral_(B_(r_0)(x)) (u(x) - u(y))/(norm(x - y)^(d + alpha)) d y = integral_0^(r_0) r^(-1 - alpha) integral_(S^(d-1)) (u(x) - u(x + r xi)) d xi d r . $

球面平均的对称化是关键一步。$xi$ 服从球面均匀分布，
而均匀分布旋转不变，所以换元 $xi -> -xi$ 不改变积分：
$integral_(S^(d-1)) F(xi) d xi = integral_(S^(d-1)) F(-xi) d xi$ 对任何 $F$ 成立。
把原式与换元后的式子相加除二：

$ integral_(S^(d-1)) (u(x) - u(x + r xi)) d xi = 1/2 integral_(S^(d-1)) [2 u(x) - u(x + r xi) - u(x - r xi)] d xi . $

为什么这一步重要，看小 $r$ 时方括号的行为。对 $u$ 做泰勒展开到二阶
（$xi$ 是单位向量，各方向导数有界）：

$ u(x plus.minus r xi) = u(x) plus.minus r nabla u(x) dot xi + r^2/2 xi^T nabla^2 u(x) xi + O(r^3) . $

两个展开相加，一阶项 $plus.minus r nabla u dot xi$ 变号相消，
剩下 $u(x + r xi) + u(x - r xi) = 2 u(x) + r^2 xi^T nabla^2 u(x) xi + O(r^3)$。
于是方括号 $2 u(x) - u(x + r xi) - u(x - r xi) = O(r^2)$，
除以 $r^2$ 后有界。对比不进行对称化的情形：
原被积函数 $u(x) - u(x + r xi)$ 的一阶项 $-r nabla u dot xi$ 乘上
权重 $r^(-1-alpha)$ 后按 $r^(-alpha)$ 发散，积分只能在对称极限
（柯西主值）下理解；对称化把产生发散的一阶项整体消掉，
剩下的量在 $r = 0$ 附近有界，期望存在，不再需要极限手续。
这就是把“主值发散”变成“可积期望”的机制。

代回并把球面积分换成球面平均（$integral_(S^(d-1)) dot d xi = |S^(d-1)| dot E_(xi)$，
均匀平均的定义）：

$ integral_(B_(r_0)(x)) (u(x) - u(y))/(norm(x - y)^(d + alpha)) d y = (|S^(d-1)|)/2 integral_0^(r_0) g(r) r^(-1 - alpha) d r , quad g(r) := E_(xi)[ 2 u(x) - u(x + r xi) - u(x - r xi) ] . $

把 $r$ 的幂律权重挪进采样分布。令 $f_I (r) prop r^(1 - alpha)$
并归一化（归一化积分
$integral_0^(r_0) r^(1 - alpha) d r = r_0^(2 - alpha)\/(2 - alpha)$）：

$ f_I (r) = (2 - alpha)/(r_0^(2 - alpha)) r^(1 - alpha) , quad integral_0^(r_0) g(r) r^(-1 - alpha) d r = (r_0^(2 - alpha))/(2 - alpha) E_(r tilde f_I)[ g(r) \/ r^2 ] . $

右边由期望的定义直接展开验证：
$E_(f_I)[g(r) \/ r^2] = integral_0^(r_0) (g(r) \/ r^2) dot (2 - alpha) r^(1 - alpha) \/ r_0^(2 - alpha) d r$，
分子里 $r^(1 - alpha)$ 与分母 $r^2$ 合并成 $r^(-1 - alpha)$，正是左边的积分。
合并常数，$|S^(d-1)| \/ 2$ 乘 $r_0^(2 - alpha) \/ (2 - alpha)$，
即得主章节的内段公式（系数 $|S^(d-1)| r_0^(2 - alpha) \/ (2(2 - alpha))$
乘期望）。

*采样密度为什么这样选*：上一步的泰勒展开给
$g(r) = O(r^2)$，所以被估量 $g(r) \/ r^2$ 在 $r -> 0$ 时趋于
径向二阶导（有界），除以 $r^2$ 的除法专门用来压住原点的奇性；
密度里的 $r^(1 - alpha)$ 增长则把权重 $r^(-1 - alpha)$ 全部吸收。
被估的量在整个 $[0, r_0]$ 上有界，方差有限。
Beta 分布采样给出这条幂律：
$R = r \/ r_0 tilde "Beta"(2 - alpha, 1)$ 的密度是
$(2 - alpha) R^(1 - alpha)$（$[0, 1]$ 上），换元回 $r$ 恰是 $f_I$，
一次 Beta 样本换一次网络评估。

#insight("用日常语言读一遍")[
  内段公式做了三件事。第一，把积分变量换成“方向加距离”，
  奇异性全部集中在距离这一个变量上。第二，用对称化把
  造成发散的一阶项消掉：差 $u$ 在两个相对方向的取值相加，
  线性成分抵消，剩下的二阶差商除以 $r^2$ 后有界。
  第三，把距离上的幂律权重挪进采样分布，让被估计的量
  变成一个有界函数的普通平均。三件事合起来，
  一个需要取对称极限的奇异积分，变成了对 $(xi, r)$ 的期望，
  可以用有限次采样估计。
]

== 外段与 Caputo 的推导思路

外段（$r in [r_0, infinity)$）同样做对称化，差异只在权重常数。
对称化后的积分是 $|S^(d-1)| \/ 2 integral_(r_0)^oo g(r) r^(-1 - alpha) d r$，
这个广义积分算得出来：

$ integral_(r_0)^oo r^(-1 - alpha) d r = [r^(-alpha)\/(-alpha)]_(r_0)^oo = r_0^(-alpha) \/ alpha , $

收敛性由 $alpha > 0$ 保证。采样密度取 $f_O (r) = alpha r_0^alpha r^(-1 - alpha)$
（归一化检查：$integral_(r_0)^oo alpha r_0^alpha r^(-1 - alpha) d r = alpha r_0^alpha dot r_0^(-alpha) \/ alpha = 1$），
代入期望的定义可以验证
$integral_(r_0)^oo g(r) r^(-1 - alpha) d r = E_(f_O)[g(r)] \/ (alpha r_0^alpha)$，
于是外段公式为系数 $|S^(d-1)| r_0^(-alpha) \/ (2 alpha)$ 乘期望。
换元 $R = r_0 \/ r$ 可验证 $R tilde "Beta"(alpha, 1)$，
与内段用同一类 Beta 采样（主章节的写法）。
这条密度重尾：远处贡献被 $r^(-1 - alpha)$ 压小，但采样覆盖到无穷远。
被估量是 $g(r)$ 本身，不再除 $r^2$：远处 $r >= r_0$ 时
$2 u(x) - u(x + r xi) - u(x - r xi)$ 有界（$u$ 有界时被 $4 sup norm u$ 压住），
原点的奇性在积分分段时已经留在内段，除法没有必要；
密度恰好吸收全部 $r$ 幂，被估量有界，方差有限。

Caputo 导数的 MC 化的完整推导如下。目标积分是不含
$1 \/ Gamma(1 - gamma)$ 因子的部分（该因子是与网络无关的正常数，
论文的估计器沿用不含它的形式，它只缩放这一项损失的尺度）：

$ I := integral_0^t (t - tau)^(-gamma) (partial u)/(partial tau)(x, tau) d tau . $

*第一步（选对原函数）*：直接给 $partial_tau u$ 配原函数 $u$ 会在
$tau = t$ 处产生发散的边界项。改给差 $v(tau) := u(tau) - u(t)$ 配，
它满足 $v(t) = 0$，恰好压住奇异端点。分部积分
（$d/(d tau) (t - tau)^(-gamma) = gamma (t - tau)^(-gamma - 1)$）：

$ I = [(t - tau)^(-gamma) v(tau)]_(tau = 0)^(tau = t) - gamma integral_0^t v(tau) (t - tau)^(-gamma - 1) d tau . $

*第二步（边界项）*：$tau = t$ 端为零：$v(t) = 0$，更仔细地看
$v(tau) approx (tau - t) partial_tau u(x, t)$ 乘 $(t - tau)^(-gamma)$
后按 $(t - tau)^(1 - gamma) -> 0$ 消失（$gamma < 1$ 保证）。
$tau = 0$ 端给出 $t^(-gamma) v(0) = t^(-gamma) (u(x, 0) - u(x, t))$，
带着前面的负号进入 $I$，合计 $t^(-gamma) (u(x, t) - u(x, 0))$。

*第三步（剩余积分换元）*：把 $v = u(tau) - u(t)$ 代回并变号，
再换元 $s = t - tau$（$d tau = -d s$）：

$ gamma integral_0^t (u(t) - u(tau)) (t - tau)^(-gamma - 1) d tau = gamma integral_0^t (u(t) - u(t - s)) s^(-gamma - 1) d s = gamma integral_0^t (u(t) - u(t - s)) \/ s dot s^(-gamma) d s . $

*第四步（吸收核）*：换元 $s = tau t$ 把区间搬到 $[0, 1]$
（$d s = t d tau$，$s^(-gamma) = t^(-gamma) tau^(-gamma)$）：

$ gamma t^(1 - gamma) integral_0^1 (u(t) - u(t - tau t)) \/ (tau t) dot tau^(-gamma) d tau . $

密度 $f_(I, t)(tau) = (1 - gamma) tau^(-gamma)$ 在 $[0, 1]$ 上归一化
（$integral_0^1 tau^(-gamma) d tau = 1 \/ (1 - gamma)$），所以积分等于
$1 \/ (1 - gamma)$ 乘 $E_(tau tilde f_(I, t))[(u(t) - u(t - tau t)) \/ (tau t)]$。
合并三步：

$ I = gamma/(1 - gamma) t^(1 - gamma) E_(tau tilde f_(I, t)) [ (u(x, t) - u(x, t - tau t)) \/ (tau t) ] + (u(x, t) - u(x, 0)) \/ t^gamma , $

即主章节的公式。数值实现把 $tau$ 换成
$tau_epsilon = max{tau, epsilon_t \/ t}$（主章节的 $tau_epsilon t$）：
$tau -> 0$ 时差商的分子分母同时变小，两个几乎相等的数相减
会吃掉有效数字，理由与内段的 $r_epsilon$ 相同。

三处推导共用一个模式：
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
