#import "../../preamble.typ": *

= 第五、六篇补充：数据驱动路线的完整记录与无偏性的推导

== 第五篇：矩问题什么时候有解

主章节的 Hankel–Cholesky 构基有一个前提没有展开：什么样的测度
能保证正交多项式在 $L^2$ 中稠密（展开式收敛到目标函数）？
Ernst 等（2012）给了五个充分条件，满足其一即可：

1. *紧支撑*：$Gamma$ 有界；
2. *矩的增长受控*：$liminf_(k -> infinity) root(2k, mu_(2k)) \/ (2k) < infinity$；
3. *指数可积*：$E[e^(a abs(eta))] < infinity$ 对某个 $a > 0$；
4. *Carleman 条件*：$sum_(k=0)^infinity 1 \/ root(2k, mu_(2k)) = infinity$；
5. *Lin 条件*：密度对称正定且 $-x f'(x)\/f(x) -> infinity$。

直觉读法：五个条件都在排除“矩增长太快以至于测度不由矩决定”的
病态情形（例如对数正态分布在某些参数下矩不唯一）。第五篇的实验
刻意包含二项、Poisson 这类离散分布：它们的支撑有限，条件 1 直接
满足，这也是数据驱动基对离散分布好用的理论原因。

== 第五篇：条件数与精度的完整记录

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a5_uni.png"), image("../fig_a/x_a5_raw.png"))),
  caption: [左：2 维均匀输入（纯连续分布）的条件数：对数线性采样率
  （log-linear）平稳。右：直方图原始数据的条件数：数据只有一份
  经验直方图时，构基与采样照常工作，条件数行为与合成分布一致。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a5_f1.png"), image("../fig_a/x_a5_res.png"))),
  caption: [左：函数 $f_1 = exp(sum xi_k)$（Poisson 输入）的离散
  $ell_2$ 误差。右：电阻网络（$d = 2$，均匀分布电阻）的逼近误差：
  1000 个样本估矩、构基、加权最小二乘，误差随采样率收敛。],
)

== 第六篇：无偏性的完整推导

主章节给了“采样按 $kappa$、加权按 $1\/kappa$、乘积还原经验分布”的读法。
这里把这笔账算清楚，因为它是 induced 采样正确性的全部。

#proposition("命题（induced 采样的无偏性，差一个全局常数）")[
  设经验样本集 $Xi = {z^((1)), ..., z^((Q))}$，induced 测度
  $mu(z) = sum_j tilde(kappa)_j delta_(z^((j)))$，
  $tilde(kappa)_j = kappa(z^((j))) \/ sum_(q=1)^Q kappa(z^((q)))$，
  权重 $W_j = 1\/kappa(z^((j)))$。对任意函数 $h$：

  $ E_(z tilde mu)[W(z) h(z)] = sum_(j=1)^Q tilde(kappa)_j (h(z^((j))))/(kappa(z^((j)))) = (1)/(sum_(q=1)^Q kappa(z^((q)))) sum_(j=1)^Q h(z^((j))) = (1)/(overline(kappa)) dot (1/Q sum_(j=1)^Q h(z^((j)))) , $

  其中 $overline(kappa) = 1/Q sum_q kappa(z^((q)))$ 是 $kappa$ 在样本集上的均值。
]

#proof[
  第一个等号按定义展开期望（有限支撑上求和）；
  第二个等号把 $tilde(kappa)_j$ 的分子 $kappa(z^((j)))$ 与权重里的
  $1\/kappa(z^((j)))$ 消去；第三个等号把常数提到求和外，
  拆成 $kappa$ 均值的倒数乘经验均值。$square$
]

#insight("这笔账在说什么")[
  加权后的估计等于经验平均乘 $1\/overline(kappa)$：*形状完全正确，
  只差一个与 $h$ 无关的常数*。回归问题对这个常数不敏感，
  把 $sqrt(W) A$ 与 $sqrt(W) b$ 同乘一个常数，基追踪的解不变
  （约束与目标同时缩放）。这就是“无偏”在压缩感知语境下的准确含义：
  不是估计量等于真值，而是*偏差与目标无关，被问题的尺度不变性吸收*。

  与重要性采样的对照也更清楚了：若把 $mu$ 看成提案分布、
  经验分布看成目标，Radon–Nikodym 导数 $d omega_hat \/ d mu = 1\/kappa$
  乘上一个归一化常数。常数在蒙特卡洛里需要显式除回去，
  在回归里免费。
]

== 第六篇：稀疏恢复与函数逼近的完整记录

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a6_raw2d.png"), image("../fig_a/x_a6_expsum.png"))),
  caption: [左：直方图原始数据的恢复率随 $M$（$d = 2$）：induced
  对 CSA 与 MC 的优势在数据完全“无分布信息”时同样成立。
  右：$f = exp(sum z_i)$（二项 × Poisson 数据，$K = 20$）：
  大 $M$ 时 induced 明显占优。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_a/x_a6_cp.png"), image("../fig_a/x_a6_pde10.png"))),
  caption: [左：corner-peak 函数（二项 × Poisson 数据）：induced 的优势
  在小 $M$ 区最明显。右：Kirchhoff 板问题的 $ell_2$ 误差（$d = 10$，
  3 个 KL 模态）：此时 $mu approx omega$，induced 与 MC 相当，
  正是方法自己指出的适用边界。],
)
