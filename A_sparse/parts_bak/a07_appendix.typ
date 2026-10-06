#import "../../preamble.typ": *

= 附录

== A. 记号表

#table(
  columns: (auto, auto),
  inset: 6pt,
  align: (left, left),
  stroke: 0.4pt + luma(190),
  table.header([*记号*], [*含义*]),
  [$d$], [随机维数],
  [$y, bold(y)$], [随机变量 / 随机向量；$Gamma$ 是取值域],
  [$rho^i, rho$], [边际密度与联合密度],
  [$phi^i_n, phi_(bold(n))$], [一元与多元正交多项式基],
  [$Lambda, N$], [指标集与基函数个数],
  [$s$], [稀疏度（显著非零系数个数）],
  [$M$], [样本（配点）数],
  [$bold(Psi), bold(A), bold(D)$], [测量矩阵（$bold(D)$ 为预条件后）],
  [$bold(W), bold(P)$], [行预条件与列归一矩阵（第三篇）],
  [$delta_s$], [RIP 常数],
  [$mu, mu_L, beta_L$], [互不相干常数与逐列上界型 coherence],
  [$K_Lambda, kappa, lambda_n$], [Christoffel 函数的三种记法（各篇略异）],
  [$L(bold(n))$], [采样测度的分布常数（主定理里的乘积界）],
  [$nu_(bold(n)), mu$], [Gauss 网格离散测度与 induced 测度],
  [$sigma_(s, 1)$], [最佳 $s$-项逼近的 $ell_1$ 尾巴],
  [$rho_a, P_a$], [TL1 罚与它的总和],
  [$a_j, b_j$], [正交多项式三项递推系数],
)

== B. 关键公式速查

#boxed(title: "问题格式")[
  - *基追踪*：$arg min norm(c)_1$ s.t. $Psi c = f$；加权版 $sqrt(W) A c = sqrt(W) b$
  - *加权最小二乘*：$arg min norm(W^(1\/2) A c - W^(1\/2) f)_2^2$
  - *梯度增强*：$tilde(Phi) = mat(Phi; Phi_partial)$，
    $arg min norm(c)_1$ s.t. $W tilde(Phi) P c = W tilde(f)$
  - *TL1*：$min P_a (x) = sum ((a+1)abs(x_i))/(a + abs(x_i))$ s.t. $A x = b$
]

#boxed(title: "恢复保证", fill: rgb("#f2fbf5"), stroke: rgb("#8dc9a3"))[
  - *RIP 恢复*：$delta_s < 0.307 => norm(c^sharp - c)_2 <= C sigma_(s,1)(c)\/sqrt(s)$
  - *样本量*：$M >= L(bold(n)) C_1 s log^3(s) log(N)$，
    $L = product_i L_i (n_i)$；Jacobi：$L_i = C$；指数型：$L_i = C n^(2\/3)$
  - *MIP 恢复*：$mu < 1\/(2s - 1)$ 时精确恢复
  - *TL1*：$delta_(2s) < 1/(1 + (a+1)\/a sqrt(2))$，
    误差 $C_0 s^(-1\/2) P_a (x - x_s) + C_1 epsilon$
  - *Cohen–LS*：$kappa(N) <= delta M \/ log M$ 时无权 LS 稳定；
    Christoffel 加权后 $hat(kappa) = N$
]

#boxed(title: "三个显式构造", fill: rgb("#fff9ec"), stroke: rgb("#ecc76a"))[
  - *求积权多项式*：$psi^i_(k, n) = sqrt(n lambda^i_n (x^i)) phi^i_k (x^i)$，
    在离散测度 $nu^i_n$ 下正交归一
  - *加权基*：$Q = "span"{psi_alpha \/ sqrt(K_Lambda)}$，行范数恒为 1
  - *induced 测度*：$mu = sum_j tilde(kappa)_j delta_(z^((j)))$，
    $tilde(kappa)_j prop kappa(z^((j)))$；采样从 $mu$、权重 $1\/kappa$
]

== C. 六篇论文信息

- L. Guo, A. Narayan, T. Zhou, Y. Chen. *Stochastic collocation methods via
  $ell_1$ minimization using randomized quadratures.* SIAM J. Sci. Comput.,
  2017. arXiv:1602.00995.
- L. Guo, A. Narayan, L. Yan, T. Zhou. *Weighted approximate Fekete points:
  sampling for least-squares polynomial approximation.* 2017.
  arXiv:1708.01296.
- L. Guo, A. Narayan, T. Zhou. *A gradient enhanced $ell_1$-minimization for
  sparse approximation of polynomial chaos expansions.* 2018.
  arXiv:1802.08837.
- L. Guo, J. Li, Y. Liu. *Stochastic collocation methods via minimization of
  transformed $ell_1$ penalty.* 2018. arXiv:1805.05416.
- L. Guo, Y. Liu, T. Zhou. *Data-driven polynomial chaos expansions:
  a weighted least-square approximation.* 2018. arXiv:1805.10893.
- L. Guo, A. Narayan, Y. Liu, T. Zhou. *Sparse approximation of data-driven
  polynomial chaos expansions: an induced sampling approach.* 2020.
  arXiv:2008.10121.

=== 背景文献

- D. Doostan, H. Owhadi. *A non-adapted sparse approximation of PDEs with
  stochastic inputs.* J. Comput. Phys., 2011. 压缩感知进入 UQ 的最早工作之一。
- L. Yan, L. Guo, D. Xiu. *Stochastic collocation algorithms using $ell_1$-minimization.*
  SIAM J. Sci. Comput., 2018. 原密度采样路线的代表。
- H. Rauhut, R. Ward. *Interpolation via weighted $ell_1$ minimization.*
  Appl. Comput. Harmon. Anal., 2016. Chebyshev 偏置采样的 RIP 理论基础。
- G. Tang, G. Iyengar (Tang–Zhou). *Compressed sensing with structured
  Gauss quadrature grids.* 2016. 本线第一篇的直接前身（均匀分布情形）。
- A. Cohen, M. A. Davenport, D. Leviatan. *On the stability and accuracy of
  least squares approximations.* 2013. 无权最小二乘 $M tilde N^2$ 的奠基分析。
- A. Cohen, G. Migliorati. *Optimal weighted least-squares methods.*
  2017. induced distribution 概念的来源（密度已知版）。
- J. Hampton, A. Doostan. *Coherence motivated sampling.* 2015.
  coherence 型采样判据的来源。
- J. Jakeman, A. Narayan, et al. *Generalized compressed sensing statistics.*
  2016. CSA（平衡测度采样）的代表。
- O. G. Ernst, et al. *On the convergence of generalized polynomial chaos
  expansions.* 2012. gPC 收敛与矩问题理论。
- X. Wan, G. E. Karniadakis. *Multi-element generalized polynomial chaos.*
  2006. ME-gPC，任意分布 gPC 的早期替代路线。
- L. Zhang, J. Xin. *T$ell_1$: transformed $ell_1$ penalty.* 2016.
  TL1 罚的原始文献。

== D. 自检清单

#boxed(title: "读完后应该能回答")[
  1. 随机配点法为什么天然欠定？稀疏性假设从哪里来？（第 1 章）
  2. RIP 与 MIP 各自衡量什么？为什么梯度增强那篇走 MIP 路线？（第 4 章）
  3. Christoffel 函数的分母大意味着什么？为什么“多采 $K$ 大的点”？
     （第 1 章定义框）
  4. Gauss 网格子采样里，求积权为什么恰好等于 Christoffel 权重？
     离散正交性怎么证？（第 6 章引理）
  5. $L_i$ 在有界与无界分布下的界差在哪里？为什么无界情形的
     $n^(2\/3)$ 无法消除？（第 6 章）
  6. 加权 Fekete 点的“贪心 = 全局最优”证明里，Hadamard 不等式用在何处？
     （第 7 章）
  7. 拼接导数行后矩阵为什么变差？$W\/P$ 预条件靠哪条正交恒等式设计？
     （第 9 章）
  8. TL1 的参数 $a$ 怎么插值 $ell_0$ 与 $ell_1$？自适应版怎么选 $a$？（第 10 章）
  9. Hankel–Cholesky 构基相对直接解矩匹配方程，稳定在哪里？（第 13 章）
  10. induced 测度怎么定义？“从 $mu$ 采样 + 按 $1\/kappa$ 加权”
      为什么等于一次重要性采样？（第 14 章）
]
