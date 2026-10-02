#import "../../preamble.typ": *

= NFF 的两笔推导：低秩高斯密度与后验

主章节用了“Woodbury 恒等式化成 $M times M$ 运算”与“后验有闭式”
两句话。这里把两笔账都算完。它们都是标准线性代数与高斯公式，
但值得亲手走一遍，因为 NFF 的全部计算优点都在里面。

== 低秩加对角协方差的对数密度

参考场在 snapshot 的 $N$ 个传感器位置上给出向量
$bold(z) tilde cal(N)(bold(mu), bold(Sigma))$，其中

$ bold(Sigma) = bold(B) bold(B)^T + bold(D)^2 , quad bold(B) = mat(B(x_1); dots.v; B(x_N)) in RR^(N times M), quad bold(D) = "diag"(C(x_1), dots.c, C(x_N)) . $

直接求 $bold(Sigma)^(-1)$ 是 $O(N^3)$，$N$ 大时不可行。
设 $N >> M$。两个恒等式：

#lemma("Woodbury 恒等式与行列式引理")[
  $ (bold(D)^2 + bold(B) bold(B)^T)^(-1) = bold(D)^(-2) - bold(D)^(-2) bold(B) (bold(I)_M + bold(B)^T bold(D)^(-2) bold(B))^(-1) bold(B)^T bold(D)^(-2) , $

  $ det(bold(D)^2 + bold(B) bold(B)^T) = det(bold(D)^2) dot det(bold(I)_M + bold(B)^T bold(D)^(-2) bold(B)) . $
]

第一个等式两边右乘 $bold(D)^2 + bold(B) bold(B)^T$ 直接展开验证；
第二个由 $det(bold(I) + bold(A) bold(B)) = det(bold(I) + bold(B) bold(A))$
（Sylvester 行列式恒等式）取
$bold(A) = bold(D)^(-2) bold(B)$、$bold(B) = bold(B)^T bold(D)^2$ 即得。

#proof[
  对数密度的最终形态：

  $ log cal(N)(bold(k) | bold(mu), bold(Sigma)) = -1/2 [ N log 2 pi + sum_i 2 log C(x_i) + bold(tilde(e))^T bold(tilde(e)) - log det(bold(I)_M + bold(B)^T bold(D)^(-2) bold(B)) ] , $

  其中 $bold(tilde(e)) = bold(D)^(-1)(bold(k) - bold(mu)) - bold(D)^(-1) bold(B) bold(w)$，
  $bold(w) = (bold(I)_M + bold(B)^T bold(D)^(-2) bold(B))^(-1) bold(B)^T bold(D)^(-2) bold(D)^(-1)(bold(k) - bold(mu))$。
  全部矩阵运算的规模：$bold(B)^T bold(D)^(-2) bold(B)$ 是 $M times M$
  （$O(N M^2)$），它的逆与行列式是 $O(M^3)$。
  总代价 $O(N M^2 + M^3)$，与 $N$ 线性。$square$
]

直觉读法：对角噪声 $bold(D)^2$ 让“白化残差”$bold(D)^(-1)(bold(k) - bold(mu))$
逐点独立；低秩部分 $bold(B)$ 只在 $M$ 维子空间里修正相关性，
Woodbury 把修正量限制在 $M$ 维内求解。KL 截断的 $M$ 固定，
$N$ 增长时代价线性，这是 NFF 敢用大 snapshot 的原因。

== 后验公式从哪里来

给定观测 $bold(k)_(1:n)$，要推断 $bold(xi) tilde cal(N)(0, bold(I)_M)$。
参考场在这些位置的均值是 $bold(B)_n bold(xi)$（把 $bold(A), bold(C)$
暂记为零的简化情形），于是观测模型是线性高斯：

$ bold(k)_(1:n) = bold(B)_n bold(xi) + bold(D)_n bold(epsilon), quad bold(epsilon) tilde cal(N)(0, bold(I)) . $

#proposition("命题（线性高斯模型的后验）")[
  $ bold(xi) | bold(k)_(1:n) tilde cal(N)(bold(mu)_xi, bold(Sigma)_xi) , quad bold(Sigma)_xi = (bold(I)_M + bold(B)_n^T bold(D)_n^(-2) bold(B)_n)^(-1) , $

  $ bold(mu)_xi = bold(Sigma)_xi bold(B)_n^T bold(D)_n^(-2) bold(k)_(1:n) . $
]

#proof[
  联合分布 $(bold(xi), bold(k))$ 是高斯，条件分布公式
  （贝叶斯公式的线性高斯版本）给出
  $bold(mu)_(bold(xi)|bold(k)) = bold(Sigma)_bold (xi bold(k)) bold(Sigma)_(bold(k) bold(k))^(-1) bold(k)$
  与 $bold(Sigma)_(bold(xi)|bold(k)) = bold(I) - bold(Sigma)_bold (xi bold(k)) bold(Sigma)_(bold(k) bold(k))^(-1) bold(Sigma)_(bold(k) bold(xi))$。
  代入 $bold(Sigma)_bold (xi bold(k)) = bold(B)_n^T$、
  $bold(Sigma)_(bold(k) bold(k)) = bold(B)_n bold(B)_n^T + bold(D)_n^2$，
  再对第二个量用 Woodbury 恒等式化简：
  $(bold(B)_n bold(B)_n^T + bold(D)_n^2)^(-1) bold(B)_n = bold(D)_n^(-2) bold(B)_n (bold(I) + bold(B)_n^T bold(D)_n^(-2) bold(B)_n)^(-1)$，
  整理即得。完整观测模型（含均值场 $bold(A)$ 与对角噪声以外的部分）
  把 $bold(k)$ 换成逆变换后的 $tilde(bold(Z))$，形式不变。$square$
]

#insight("后验公式在说什么")[
  $bold(Sigma)_xi$ 的形状是“先验单位阵减去一个半正定修正”：
  每条观测从先验里减掉一点不确定性，减多少由
  $bold(B)_n^T bold(D)_n^(-2) bold(B)_n$（信息矩阵）决定。
  观测越多、噪声越小，修正越大，后验越收缩到点估计；
  没有观测时退回 $bold(I)_M$，即先验。这个“自动收缩”行为
  就是主章节说“退回先验”的公式形态。
]

== 标量场怎么办

RealNVP 耦合块要求维度 $>= 2$（要把输入拆成两半）。
一维随机场（$D = 1$）的处理：给场值附加一个独立的标准高斯
“陪跑”变量 $v tilde cal(N)(0, 1)$，把 $(k, v)$ 拼成二维做 flow，
推断时 $v$ 的边缘解析已知，可以直接积掉。
这不是 elegance 问题而是实现要点：一维是最常见的标量 QoI 场景。

== NFF 实验画廊

#figure(
  align(center, image("../fig_b/x_b3_mgp.png", width: 88%)),
  caption: [混合双模态非高斯场（一半概率加正弦、一半概率减正弦）：
  生成路径的均值、标准差与协方差谱（论文 2108.12956）。
  按两个条件模态分别对比，双峰结构被 flow 精确表达。],
)

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_b/x_b3_infer.png"), image("../fig_b/x_b3_invhyb.png"))),
  caption: [左：单点测量的后验推断：后验均值与方差包住真值，
  验证闭式后验公式。右：一维随机椭圆方程反向与混合问题的误差
  （$N_k = 1$ 或 $3$ 个 $k$ 传感器 + 若干 $u$ 传感器），
  与蒙特卡洛参考同量级。],
)
