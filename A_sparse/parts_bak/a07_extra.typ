#import "@preview/cetz:0.3.4": canvas, draw

#import "../../preamble.typ": *

= 附录补充：术语索引、阅读路线与跨界接口

== 术语索引（按定义位置）

#ttable(
  ([术语], [定义位置], [一句话]),
  (
    ([Christoffel 函数], [第 1 章定义框], [多项式平方和的倒数（或其倒数版），采样与加权的公共主轴]),
    ([三项递推 / Jacobi 矩阵], [第 2 章], [递推系数编码测度；Golub–Welsch 从它读出节点与权]),
    ([Gauss 求积精度], [第 2 章，定理], [$n$ 点法则对 $2n - 1$ 次精确，证明一行带余除法]),
    ([求积权 = Christoffel 函数], [第 2 章，引理], [Christoffel–Darboux 的直接推论]),
    ([基追踪], [第 1 章], [$min norm(c)_1$ s.t. $Psi c = f$，$ell_0$ 的凸松弛]),
    ([RIP], [第 4 章，定义], [矩阵近似保持稀疏向量长度；随机采样高概率成立]),
    ([MIP / MIC], [第 4 章，定义], [最坏列对的归一化内积；逐点可验]),
    ([一致有界系统 / $L$ 因子], [第 6 章，引理], [加权多项式逐点界的常数，决定样本量]),
    ([加权近似 Fekete 点], [第 7 章], [$1\/sqrt(K)$ 进基的贪心点列；加权后贪心全局最优]),
    ([梯度增强], [第 9 章], [导数行拼入测量矩阵；$W\/P$ 预条件保各向同性]),
    ([TL1 罚], [第 10 章], [$rho_a = (a+1)abs(t)\/(a+abs(t))$，$a$ 插值 $ell_0$ 与 $ell_1$]),
    ([Hankel–Cholesky 构基], [第 13 章], [矩 → 对称正定矩阵 → 分解 → 递推系数]),
    ([平衡测度], [第 14 章与谱系章], [$kappa rho$ 的位势论极限；有界域上是 Chebyshev 测度]),
    ([induced 测度], [第 14 章], [经验样本按 $kappa$ 加权的离散测度；配 $1\/kappa$ 权无偏]),
    ([杠杆分数], [谱系章], [随机数值线性代数里同一思想的另一个名字]),
  ),
  widths: (1.15fr, 0.8fr, 2.55fr),
)

== 与方向 B 的接口：数据驱动思想进入深度学习

组内方向 B（神经网络解随机方程）的第一篇（NN-aPC）与本线的
第五、六篇共享同一个出发点：*输入分布未知、只有样本时，
正交基从数据的矩构造*。差异在下游：本线把基交给多项式逼近，
方向 B 把基交给网络的模态结构。两条线对“数据驱动基”的分工
也互补：本线的收敛性由逼近论保证（多项式最优逼近率），
方向 B 的表达力由网络保证（非线性特征），代价是没有收敛率。
一个正在发生的合流：学习式基函数 + 稀疏恢复系数
（用方向 B 的网络算基、用本线的 $ell_1$ 算系数）在文献里
已经出现雏形，本线第六篇的 induced 采样正是这类混合方法
需要的采样部件。

== 阅读路线

#align(center, canvas({
  import draw: *
  import "../../preamble.typ": cnode
  let c1 = rgb("#1a4f9a")
  let c2 = rgb("#b8860b")
  let cg = rgb("#4a7c59")
  content((0, 0), cnode(4.4cm, 1.0cm, c1, [问题与工具（第 1 到 4 章）]))
  content((0, -1.6), cnode(4.4cm, 1.0cm, c1, [采样与选点（第 6 到 7 章）]))
  content((0, -3.2), cnode(4.4cm, 1.0cm, c1, [导数与非凸罚（第 9 到 10 章）]))
  content((0, -4.8), cnode(4.4cm, 1.0cm, c1, [数据驱动（第 13、14 章）]))
  content((0, -6.4), cnode(4.4cm, 1.0cm, cg, [案例研究（第 20 章）]))
  for i in range(4) {
    line((0, -0.55 - 1.6 * i), (0, -1.05 - 1.6 * i), stroke: 1pt + c1, mark: (end: ">", fill: c1))
  }
  line((0, -5.35), (0, -5.85), stroke: 1pt + cg, mark: (end: ">", fill: cg))
  content((4.7, -1.6), cnode(4.0cm, 1.0cm, c2, [只关心理论：第 2、3、5 章], size: 8pt))
  content((4.7, -3.2), cnode(4.0cm, 1.0cm, c2, [只关心上手：第 18 章], size: 8pt))
  content((4.7, -4.8), cnode(4.0cm, 1.0cm, c2, [只关心选型：第 20 章], size: 8pt))
}))

== 参考文献补充（谱系章新增条目）

- N. Wiener. *The homogeneous chaos.* American J. Math., 1938.
- R. Ghanem, P. Spanos. *Stochastic finite elements: a spectral approach.*
  Springer, 1991.
- S. Smolyak. *Quadrature and interpolation formulas for tensor products
  of certain classes of functions.* 1963.
- Gerstner, Griebel. *Numerical integration using sparse grids.* 1998.
- I. Babuška, F. Nobile, R. Tempone. *A stochastic collocation method for
  elliptic PDEs with random input data.* SIAM Review, 2010.
- D. Xiu, G. E. Karniadakis. *The Wiener–Askey polynomial chaos.* 2002.
- X. Wan, G. E. Karniadakis. *Multi-element generalized polynomial chaos.*
  J. Comput. Phys., 2006.
- R. Ahlfeld, B. Belkharraz, et al. (SAMBA). 2016.
- A. Cohen, M. A. Davenport, D. Leviatan. *On the stability and accuracy
  of least squares approximations.* Found. Comput. Math., 2013.
- A. Chkifa, A. Cohen, C. Schwab. *High-dimensional adaptive sparse
  polynomial interpolation.* 2014.
- A. Cohen, G. Migliorati. *Optimal weighted least-squares methods.*
  2017.
- J. Hampton, A. Doostan. *Coherence motivated sampling.* J. Comput.
  Phys., 2015.
- J. Jakeman, A. Narayan, et al. *Generalized compressed sensing
  statistics.* SIAM/ASA JUQ, 2016.
- A. Cohen, R. DeVore, C. Schwab. *Analytic regularity and polynomial
  approximation of parametric stochastic PDEs.* 2010. 稀疏性定理的出处。
- E. Fornberg et al. SPGL1 求解器：M. Friedlander, M. van den Berg.
  *Probing the Pareto frontier for basis pursuit solutions.* 2008.
