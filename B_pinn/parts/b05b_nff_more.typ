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

两个等式都值得把中间矩阵运算写全。记
$bold(K) := bold(B)^T bold(D)^(-2) bold(B)$（$M times M$，这是后面反复出现的量）。

*逆公式的验证*：把第一个等式的右端记为 $bold(M)_0$，验证
$bold(M)_0 (bold(D)^2 + bold(B) bold(B)^T) = bold(I)_N$。逐项算乘积。
第一项：

$ bold(D)^(-2) (bold(D)^2 + bold(B) bold(B)^T) = bold(I)_N + bold(D)^(-2) bold(B) bold(B)^T . $

第二项（先把 $bold(D)^(-2)$ 分配进去，再凑出 $bold(K)$）：

$ bold(D)^(-2) bold(B) (bold(I)_M + bold(K))^(-1) bold(B)^T bold(D)^(-2) (bold(D)^2 + bold(B) bold(B)^T) = bold(D)^(-2) bold(B) (bold(I)_M + bold(K))^(-1) (bold(B)^T + bold(B)^T bold(D)^(-2) bold(B) bold(B)^T) . $

第二处把 $bold(B)^T bold(D)^(-2) bold(B) = bold(K)$ 从括号里认出，
两个括号相乘消掉：

$ = bold(D)^(-2) bold(B) (bold(I)_M + bold(K))^(-1) (bold(I)_M + bold(K)) bold(B)^T = bold(D)^(-2) bold(B) bold(B)^T . $

第三处用了 $bold(B)^T bold(D)^(-2) bold(B) bold(B)^T = bold(K) bold(B)^T$
（$bold(K)$ 从左侧乘入）。两项相减：
$bold(I)_N + bold(D)^(-2) bold(B) bold(B)^T - bold(D)^(-2) bold(B) bold(B)^T = bold(I)_N$。
再验证左乘也得到 $bold(I)$，或者直接注意 $bold(Sigma)$ 对称时右逆即左逆。
逆公式成立。

*行列式引理的验证*：Sylvester 行列式恒等式说
$det(bold(I)_N + bold(P) bold(Q)) = det(bold(I)_M + bold(Q) bold(P))$，
其中 $bold(P)$ 是 $N times M$、$bold(Q)$ 是 $M times N$
（证明思路是把 $bold(I)_N + bold(P)bold(Q)$ 与分块矩阵消元联系起来，这里引用）。
取 $bold(P) := bold(B)$（$N times M$）、$bold(Q) := bold(B)^T bold(D)^(-2)$（$M times N$）。
一侧是 $bold(I)_M + bold(Q) bold(P) = bold(I)_M + bold(K)$；另一侧提出右因子
$bold(D)^(-2)$（$bold(D)^2 bold(D)^(-2) = bold(I)$）：

$ bold(I)_N + bold(P) bold(Q) = bold(I)_N + bold(B) bold(B)^T bold(D)^(-2) = (bold(D)^2 + bold(B) bold(B)^T) bold(D)^(-2) . $

两边取行列式并用 $det(bold(X) bold(Y)) = det(bold(X)) det(bold(Y))$、
$det(bold(D)^(-2)) = 1 \/ det(bold(D)^2)$：

$ det(bold(I)_M + bold(K)) = det(bold(D)^2 + bold(B) bold(B)^T) \/ det(bold(D)^2) , $

移项即得第二个等式。

#proof[
  分三块走：对数密度拆成三部分，每部分各自降到 $M$ 维运算。

  *第一步（拆块）*：多维高斯密度的定义给出

  $ log cal(N)(bold(k) | bold(mu), bold(Sigma)) = -1/2 [ N log 2 pi + log det bold(Sigma) + (bold(k) - bold(mu))^T bold(Sigma)^(-1) (bold(k) - bold(mu)) ] , $

  三部分分别是归一化常数、行列式、二次型。最后一项是唯一含数据的块。

  *第二步（行列式块）*：用行列式引理与 $det bold(D)^2 = product_i C(x_i)^2$：

  $ log det bold(Sigma) = log det bold(D)^2 + log det bold(I)_M + bold(K) = sum_i 2 log C(x_i) + log det bold(I)_M + bold(K) . $

  *第三步（二次型块，先代入 Woodbury）*：记
  $bold(r) := bold(D)^(-1)(bold(k) - bold(mu))$（“白化残差”）、
  $bold(A) := bold(D)^(-1) bold(B)$（于是 $bold(A)^T bold(A) = bold(B)^T bold(D)^(-2) bold(B) = bold(K)$）。
  把 Woodbury 公式的 $bold(Sigma)^(-1)$ 代入二次型并逐项识别：

  $ (bold(k) - bold(mu))^T bold(Sigma)^(-1) (bold(k) - bold(mu)) = bold(r)^T bold(r) - bold(r)^T bold(A) (bold(I)_M + bold(K))^(-1) bold(A)^T bold(r) , $

  因为修正项的三个因子依次是
  $bold(r)^T bold(A)$、$(bold(I)_M + bold(K))^(-1)$、$bold(A)^T bold(r)$
  （$bold(D)^(-1)$ 对称，转置时左右各提一个出来：
  $(bold(k) - bold(mu))^T bold(D)^(-2) bold(B) = (bold(D)^(-1)(bold(k) - bold(mu)))^T (bold(D)^(-1) bold(B)) = bold(r)^T bold(A)$，
  右端因子 $bold(B)^T bold(D)^(-2) (bold(k) - bold(mu)) = bold(A)^T bold(r)$ 同理）。

  *第四步（配方）*：令 $bold(w) = (bold(I)_M + bold(K))^(-1) bold(A)^T bold(r)$，
  断言二次型等于 $bold(tilde(e))^T bold(tilde(e)) + bold(w)^T bold(w)$，
  其中 $bold(tilde(e)) := bold(r) - bold(A) bold(w)$。验证：展开平方项，

  $ bold(tilde(e))^T bold(tilde(e)) + bold(w)^T bold(w) = bold(r)^T bold(r) - 2 bold(w)^T bold(A)^T bold(r) + bold(w)^T (bold(A)^T bold(A)) bold(w) + bold(w)^T bold(w) = bold(r)^T bold(r) - bold(w)^T bold(A)^T bold(r) , $

  中间一步把 $-2 bold(w)^T bold(A)^T bold(r) = -2 bold(w)^T (bold(I)_M + bold(K)) bold(w)$
  （$bold(w)$ 的定义移项）与 $bold(w)^T (bold(I)_M + bold(K)) bold(w) + bold(w)^T bold(K) bold(w)$
  合并。最后用 $(bold(I)_M + bold(K)) bold(w) = bold(A)^T bold(r)$ 再移一次项，
  得到的正是第三步的二次型。

  *第五步（合并）*：把两块代回第一步：

  $ log cal(N)(bold(k) | bold(mu), bold(Sigma)) = -1/2 [ N log 2 pi + sum_i 2 log C(x_i) + bold(tilde(e))^T bold(tilde(e)) + bold(w)^T bold(w) + log det(bold(I)_M + bold(B)^T bold(D)^(-2) bold(B)) ] , $

  其中 $bold(tilde(e)) = bold(D)^(-1)(bold(k) - bold(mu)) - bold(D)^(-1) bold(B) bold(w)$，
  $bold(w) = (bold(I)_M + bold(B)^T bold(D)^(-2) bold(B))^(-1) bold(B)^T bold(D)^(-2) (bold(k) - bold(mu))$。
  翻译：对数密度等于“白化残差减去低秩修正后的平方”加“修正量自身的平方”，
  再加对角噪声与低秩修正各自的行列式贡献。

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
  用贝叶斯公式加配方，两步走完，不需要引用条件分布公式。

  *第一步（写出后验密度）*：贝叶斯公式说后验密度正比于先验乘似然：

  $ p(bold(xi) | bold(k)_(1:n)) prop underbrace(e^(-1/2 bold(xi)^T bold(xi)), "先验" cal(N)(bold(0), bold(I)_M)) dot underbrace(e^(-1/2 (bold(k)_(1:n) - bold(B)_n bold(xi))^T bold(D)_n^(-2) (bold(k)_(1:n) - bold(B)_n bold(xi))), "似然") . $

  似然的来历：观测模型 $bold(k)_(1:n) = bold(B)_n bold(xi) + bold(D)_n bold(epsilon)$
  里噪声项是 $bold(epsilon) = bold(D)_n^(-1)(bold(k)_(1:n) - bold(B)_n bold(xi))$，
  它服从 $cal(N)(bold(0), bold(I))$，写出这个标准高斯的密度就是第二个因子。
  比例常数不含 $bold(xi)$，配方时先不管。

  *第二步（按 $bold(xi)$ 展开）*：指数上的和只保留含 $bold(xi)$ 的项
  （平方展开产生三类：$bold(xi)^T bold(xi)$、交叉项、与 $bold(xi)$ 无关的常数）：

  $ bold(xi)^T bold(xi) + (bold(k) - bold(B)_n bold(xi))^T bold(D)_n^(-2) (bold(k) - bold(B)_n bold(xi)) = bold(xi)^T (bold(I)_M + bold(B)_n^T bold(D)_n^(-2) bold(B)_n) bold(xi) - 2 bold(xi)^T bold(B)_n^T bold(D)_n^(-2) bold(k)_(1:n) + bold("const.") , $

  其中展开 $bold(xi)^T bold(B)_n^T bold(D)_n^(-2) bold(B)_n bold(xi)$ 与
  $-2 bold(xi)^T bold(B)_n^T bold(D)_n^(-2) bold(k)$ 来自交叉项
  （$bold(D)_n^(-2)$ 对称）。记 $bold(K)_n := bold(B)_n^T bold(D)_n^(-2) bold(B)_n$。

  *第三步（配方）*：二次式加线性式的标准配方（矩阵版本：
  $bold(xi)^T bold(A) bold(xi) - 2 bold(xi)^T bold(b) = (bold(xi) - bold(A)^(-1) bold(b))^T bold(A) (bold(xi) - bold(A)^(-1) bold(b)) - bold(b)^T bold(A)^(-1) bold(b)$，
  $bold(A)$ 对称正定时成立）给出

  $ bold(xi)^T (bold(I)_M + bold(K)_n) bold(xi) - 2 bold(xi)^T bold(B)_n^T bold(D)_n^(-2) bold(k)_(1:n) = (bold(xi) - bold(mu)_bold(xi))^T (bold(I)_M + bold(K)_n) (bold(xi) - bold(mu)_bold(xi)) + bold("const.") , $

  $ bold(mu)_bold(xi) = (bold(I)_M + bold(K)_n)^(-1) bold(B)_n^T bold(D)_n^(-2) bold(k)_(1:n) . $

  *第四步（读出后验）*：于是后验密度是
  $e^(-1\/2 (bold(xi) - bold(mu)_bold(xi))^T (bold(I)_M + bold(K)_n) (bold(xi) - bold(mu)_bold(xi)))$
  乘常数，这是均值为 $bold(mu)_bold(xi)$、协方差为
  $(bold(I)_M + bold(K)_n)^(-1)$ 的高斯密度（二次型矩阵的逆即协方差），
  即命题的两个式子。$square$
]

同一条结论也可以走条件分布公式的路：联合分布
$(bold(xi), bold(k))$ 是高斯，条件分布公式给出
$bold(mu)_(bold(xi)|bold(k)) = bold(Sigma)_(bold(xi) bold(k)) bold(Sigma)_(bold(k) bold(k))^(-1) bold(k)$
与 $bold(Sigma)_(bold(xi)|bold(k)) = bold(I) - bold(Sigma)_(bold(xi) bold(k)) bold(Sigma)_(bold(k) bold(k))^(-1) bold(Sigma)_(bold(k) bold(xi))$；
代入 $bold(Sigma)_(bold(xi) bold(k)) = bold(B)_n^T$、
$bold(Sigma)_(bold(k) bold(k)) = bold(B)_n bold(B)_n^T + bold(D)_n^2$，
对 $bold(Sigma)_(bold(k) bold(k))^(-1)$ 用 Woodbury 恒等式
$bold(Sigma)_(bold(k) bold(k))^(-1) bold(B)_n = bold(D)_n^(-2) bold(B)_n (bold(I)_M + bold(K)_n)^(-1)$
（把上面的逆公式两边右乘 $bold(B)_n$ 即得），
两条路算的是同一个后验。完整观测模型（含均值场 $bold(A)$ 与对角噪声以外的部分）
把 $bold(k)$ 换成逆变换后的 $tilde(bold(Z))$，形式不变。

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
这一步是实现要点：一维是最常见的标量 QoI 场景。

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
