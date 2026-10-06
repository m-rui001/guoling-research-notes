#import "../../preamble.typ": *

= 第三、四篇补充：梯度信息的完整记录与非凸罚的证明链

== 第三篇：三族测试函数与高维部分梯度

主章节给了 MIC 对比与恢复率。这里补齐函数逼近与 PDE 的完整记录，
外加一个重要的推广讨论。

*三个测试函数*（$d = 2$，$n = 20$，标准 $ell_1$ 与梯度增强的对比）：
球面函数 $f_1 = sum_i x_i^2$（梯度线性，信息量最大）、
高斯脉冲 $f_2 = exp(-sum_i 0.01 (1\/2 (x_i + 1) - 0.375)^2)$
（梯度局部化）、多频正弦 $f_3 = sum_i 0.3 + sin(16\/15 x_i - 0.7) + sin^2(dot)$
（梯度振荡）。三者的共同结论：梯度增强把 RMSE 压低一半以上，
对梯度结构最不敏感的是 $f_1$、最敏感的是 $f_3$。

*梯度信息的两种“部分可用”*：实际工程里往往只有一部分参数的梯度
（比如只对边界条件有伴随），或只有方向导数
$D_(bold(v)_t) f(bold(z)^((j))) = ⟨nabla f(bold(x)), bold(v)_t ⟩ |_(bold(x) = bold(z)_j)$
（某个物理方向的敏感度）。论文的框架对两者都开放：
测量矩阵的“导数行”可以只拼入可用方向，预条件矩阵按可用方向缩放。
$d = 10$（$n = 3$，$s = 6$）的部分梯度实验显示恢复率随可用方向的
比例单调上升，全梯度与 10% 梯度之间的差距主要在小样本区。

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a3_sphere.png"), image("../fig_a/x_a3_sin.png"))),
  caption: [测试函数逼近误差（$d = 2$）。左：球面函数（梯度线性）；
  右：多频正弦（梯度振荡，梯度增强收益最大的情形）。
  两图中梯度增强（GraLeCh）都显著低于同信息量的标准 $ell_1$。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a3_micpce.png"), image("../fig_a/x_a3_pde10.png"))),
  caption: [左：MIC 随 gPC 项数的变化（$d = 3$，$n = 80$）：预条件后的
  梯度增强矩阵在全部项数上保持最低常数，验证预条件的尺度无关性。
  右：随机椭圆 PDE（$d = 10$）的均值误差：梯度由伴随方程获得，
  梯度增强全面占优。],
)

== 第四篇：TL1 带噪恢复的证明链

主章节给的是无噪定理。带噪版本（定理 3.6）的证明链值得完整记录，
因为它是 Candès RIP 分析模板在非凸罚上的移植，读一遍就能看出
“哪一步依赖凸性、哪一步不依赖”。

#proof[
  记 $bold(h) = hat(bold(x)) - bold(x)$，目标是 $norm(bold(h))_2 <= C_0 s^(-1\/2) P_a (bold(x) - bold(x)_s) + C_1 epsilon$。

  *第一步（罚的单调性）。* $hat(bold(x))$ 是最优解，故
  $P_a (hat(bold(x))) <= P_a (bold(x))$。用三角不等式型引理
  $rho_a (abs(u)) - rho_a (abs(v)) <= rho_a (abs(u - v))$ 把真解的贡献拆走，
  得到对残差的约束 $P_a (bold(h)_(T_0^c)) <= P_a (bold(h)_(T_0)) + 2 P_a (bold(x)_(T_0^c))$，
  其中 $T_0$ 是 $bold(x)$ 最大 $s$ 个系数的指标集。
  （这一步用到 $rho_a$ 的凹性与次可加性，不依赖凸性。）

  *第二步（尾部分块和）。* 把 $T_0^c$ 的系数按大小分成
  $s$ 个一组的块 $T_1, T_2, dots$，逐块比较得
  $sum_(j >= 2) norm(bold(h)_(T_j))_2 <= s^(-1\/2) P_a (bold(h)_(T_0^c))$。
  （这一步对任何满足该型不等式的罚都成立，$ell_1$ 时退化为
  Candès 的经典分块。）

  *第三步（RIP 控制块交互）。* 约束 $B hat(x) - b$ 的范数不超过 $epsilon$
  给出 $norm(B bold(h))_2 <= epsilon$；把 $bold(h)_(T_0)$ 与
  $bold(h)_(T_0^c)$ 用 RIP 展开、Cauchy–Schwarz 配合
  $delta_(2s) < 1/(1 + (a+1)\/a sqrt(2))$ 吸收交叉项，
  解出 $norm(bold(h)_(T_0))_2$ 的上界，代回第二步合并。

  *第四步（误差项）。* 约束右端的 $epsilon$ 沿第三步的链路传到
  最终界，系数 $C_1$ 来自 RIP 下界对 $norm(B bold(x)_(T_0^c))_2$ 的放大。
  $square$
]

这个证明链的教学价值：非凸罚没有改写恢复理论的骨架，
改写的只是第一步与第二步里“范数不等式”的具体形状
（$P_a$ 替代 $norm(dot)_1$）。所以 TL1 与 $ell_1$ 的差别
在常数与条件，不在证明结构。

== 第四篇：参数与非稀疏目标的完整记录

#figure(
  align(center, grid(columns: 3, gutter: 8pt,
    image("../fig_a/tmp_a4_tl1.png"), image("../fig_a/x_a4_tl1op.png"), image("../fig_a/x_a4_tl1big.png"))),
  caption: [TL1 罚的等值线随参数 $a$ 的变化：$a = 1$（左）、
  $a = 0.1$（中，坐标轴附近更尖，接近 $ell_0$）、$a = 100$（右，
  几乎退化为 $ell_1$ 菱形）。$a$ 是从 $ell_1$ 通向 $ell_0$ 的连续旋钮。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a4_viaa.png"), image("../fig_a/x_a4_kdv.png"))),
  caption: [左：恢复率随 $a$ 的扫描（$d = 2$，$k = 20$，$N = 231$，
  $M = 45$）：$a = 0.3$ 附近最优，过小（过近 $ell_0$，DCA 易陷局部）
  与过大（退化为 $ell_1$）都变差。右：随机 KdV 方程的相对 RMSE：
  TL1 与 $ell_(1-2)$、$ell_1$ 的对比，TL1 在两个维数上都最低。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a4_f1.png"), image("../fig_a/x_a4_cp.png"))),
  caption: [非稀疏解析函数的逼近误差：左 $f_1 = 1 \/ sum_i (0.5 + 0.1 z_i)$，
  右 corner-peak 函数。目标本身不稀疏时，非凸罚的优势缩小但排序不变
  （TL1 $>$ $ell_(1-2)$ $>$ $ell_1$）；这正是自适应选择 $a$
  （对候选集各解一次、保留最稀疏解）想要覆盖的场景。],
)
