#import "../../preamble.typ": *

= 附录

== A. 记号表

#table(
  columns: (auto, auto),
  inset: 6pt,
  align: (left, left),
  stroke: 0.4pt + luma(190),
  table.header([*记号*], [*含义*]),
  [$x in cal(D)$], [物理坐标；$cal(D)$ 是物理区域，维数 $D_x$],
  [$omega, xi$], [随机样本点；$xi$ 多指 KL/参数化后的随机变量],
  [$k(x, omega)$], [随机系数（方程参数）随机场],
  [$u(x, omega)$], [随机解],
  [$overline(u)$], [均值场 $E[u]$],
  [$u_i, Y_i, a_i$], [时变模态展开的空间模态、随机模态、振幅],
  [$⟨dot, dot⟩$], [物理空间内积 $integral cal(D) dot dot d x$],
  [$E[dot]$], [概率空间期望],
  [$cal(N)_x[u]$], [方程的非线性算子（含空间导数）],
  [$D_(|x|)^alpha$], [Riesz 分数阶导数],
  [$(-Delta)^(alpha\/2)$], [分数阶 Laplacian（超奇异积分定义）],
  [$Delta x, w_k^alpha$], [GL 离散的网格步长与递推系数],
  [$theta, phi$], [网络参数],
  [$cal(F)_"KZ"$], [NFF 的逆变换（目标场 $k$ 拉回参考场 $z$；正向是 $cal(F)_"ZK"$）],
  [$M_n$], [KL 截断阶数],
  [$hat(L), hat(L)'$], [MC 算子的两次独立估计],
  [$m$], [MC 样本数],
)

== B. 关键公式速查

#boxed(title: "表示随机性")[
  - *KL 展开*：$k(x, omega) approx overline(k)(x) + sum_(i=1)^M sqrt(lambda_i) phi_i (x) xi_i (omega)$
  - *aPC 正交化*：$psi_alpha = w_alpha^alpha psi_alpha^* - sum_(beta prec alpha) w_beta^alpha psi_beta$，
    正交性对经验测度成立
  - *时变模态展开*：$u_N = overline(u) + sum_(i=1)^N a_i (t) u_i (x, t) Y_i (t; omega)$
  - *flow 密度*：$log p_K (k) = log p_(Z^0)(z^0) - sum_(j=1)^n log |det (partial cal(F)_j)/(partial z^(j-1))|$
  - *NFF 参考场*：$z(x, omega) = A(x) + B(x) xi(omega) + "diag"(C(x)) epsilon(x, omega)$
]

#boxed(title: "PINN 类损失", fill: rgb("#f2fbf5"), stroke: rgb("#8dc9a3"))[
  - *强形式*：$cal(L) = "MSE"_u + "MSE"_k + "MSE"_f$，残差 $cal(N)_x[tilde(u); tilde(k)]$
  - *弱形式三投影*：$E[dot]$、$⟨dot, u_i⟩$、$E[dot Y_i]$
  - *DO 约束*：$⟨partial_t U_i, U_j⟩ = 0$，$E[Y_i] = 0$，$E[Y_i dot(Y)_i] = 0$
  - *BO 约束*：$⟨dot U_i, U_j⟩ + ⟨dot U_j, U_i⟩ = 0$，
    $E[Y_i dot(Y)_j] + E[Y_j dot(Y)_i] = 0$
  - *MC 无偏物理损失*：$hat(L)_"equ" = 1/(m N_u) sum hat(L) dot hat(L)'$，
    $E[hat(L)_"equ"] = L_"equ"$
  - *WGAN-GP*：$cal(L)_D = E[D(G)] - E[D(psi)] + omega E[(norm nabla D_2 - 1)^2]$
]

#boxed(title: "经典模态演化方程（NN 版本绕开它们）", fill: rgb("#fff9ec"), stroke: rgb("#ecc76a"))[
  - *DO*：$sum_i C_(Y_i Y_j) partial_t u_i = product_(V_S^perp) E[cal(N)_x[u] Y_j]$
    （要求协方差矩阵可逆）
  - *BO 闭式*：$M_(i j) = (G_(i j) + G_(j i))\/(-lambda_i + lambda_j)$
    （要求无特征值交叉）
]

== C. 六篇论文信息

- D. Zhang, L. Lu, L. Guo, G. E. Karniadakis. *Quantifying total uncertainty in
  physics-informed neural networks for solving forward and inverse stochastic
  problems.* J. Comput. Phys., 2019. arXiv:1809.08327.
- D. Zhang, L. Guo, G. E. Karniadakis. *Learning in modal space: solving
  time-dependent stochastic PDEs using physics-informed neural networks.*
  SIAM J. Sci. Comput., 2020. arXiv:1905.01205.
- L. Guo, H. Wu, T. Zhou. *Normalizing field flows: solving forward and inverse
  stochastic differential equations using physics-informed diffusion flows.*
  J. Sci. Comput., 2022. arXiv:2108.12956.
- L. Guo, H. Wu, X. Yu, T. Zhou. *Monte Carlo PINNs: deep learning approach for
  forward and inverse problems involving high dimensional fractional partial
  differential equations.* 2022. arXiv:2203.08501.
- L. Ma, R. Li, F. Zeng, L. Guo, G. E. Karniadakis. *Bi-orthogonal fPINN:
  a physics-informed neural network method for solving time-dependent
  stochastic PDEs.* J. Comput. Phys., 2023. arXiv:2303.10913.
- S. Zhou, W. You, L. Guo, X. Meng. *Scalable physics-informed deep generative
  model for solving forward and inverse stochastic differential equations.*
  2025. arXiv:2503.18012.

=== 背景文献

- M. Raissi, P. Perdikaris, G. E. Karniadakis. *Physics-informed neural
  networks.* J. Comput. Phys., 2019. PINN 的原始论文。
- D. Xiu, G. E. Karniadakis. *The Wiener--Askey polynomial chaos for stochastic
  differential equations.* SIAM J. Sci. Comput., 2002. gPC 的奠基工作。
- T. Sapsis, P. Lermusiaux. *Dynamically orthogonal field equations.*
  Physica D, 2009. DO 条件的来源。
- M. Cheng, T. Y. Hou, Z. Zhang. *A dynamically bi-orthogonal method for
  time-dependent stochastic PDEs.* J. Hyp. Diff. Eq., 2013. BO 条件的来源。
- L. Pang, K. Lu, G. E. Karniadakis. *fPINNs: Fractional physics-informed
  neural networks.* SIAM J. Sci. Comput., 2019. GL 离散的来源。
- L. Dinh, J. Sohl-Dickstein, S. Bengio. *Density estimation using Real NVP.*
  ICLR, 2017. RealNVP 耦合块的来源。
- L. Kharazmi, Z. Zhang, G. E. Karniadakis. *vpINNs: Variational
  physics-informed neural networks.* 2020. 变分方程损失的来源。
- Z. Fang, J. Li, et al. *Physics-informed GAN (PI-GAN).* SIAM J. Sci.
  Comput., 2020. 高维随机方程的早期生成式方法。

== D. 自检清单

#boxed(title: "读完后应该能回答")[
  1. 随机场的 KL 展开在做什么？有效维数由什么决定？（第 1 章）
  2. PINN 的损失函数由哪两类项构成？逆问题在这个框架下为什么不需要额外机制？（第 2 章）
  3. aPC 相对 gPC 放弃了什么假设？用什么换的？（第 4 章）
  4. 经典 DO 方法在确定性初值下为什么崩溃？NN-DO 靠什么绕开？（第 6 章）
  5. 弱形式的三条投影分别是什么？为什么它们合起来就足够确定截断解？（第 6 章）
  6. normalizing flow 的对数似然为什么能逐层精确计算？RealNVP 的雅可比
     为什么是三角的？（第 8 章）
  7. NFF 的参考场为什么要加对角噪声？（第 8 章）
  8. 分数阶 Laplacian 的 MC 估计为什么内段与外段要用不同的采样密度？（第 10 章）
  9. 随机算子的平方损失为什么必须用两次独立估计的乘积？（第 10 章）
  10. sPI-GeM 把生成模型的度量估计从什么空间搬到了什么空间？
      代价是什么？（第 13 章）
]

== E. 常见问答

#boxed(title: "六个高频问题")[
  *问：PINN 类方法能达到传统求解器的精度吗？*
  不能指望。模态方法的天花板约 $10^(-5)$，MC 类约 $10^(-3)$，
  传统谱方法可以到机器精度。PINN 系的卖点是逆问题、
  无网格、与数据融合，正向精度竞争不是主战场。

  *问：模态个数 $N$ 怎么选？*
  无理论准则。实践：从 3 到 6 起步，看各阶振幅的能量占比，
  新增模态对末端预测的提升低于噪声水平时停止。
  BO-fPINN 的 6/7/8 模态对比表就是这种收敛性检查。

  *问：MC-PINN 的采样噪声会不会让训练不稳定？*
  损失面确实带噪声，这是它精度停在 $10^(-3)$ 的原因之一。
  缓解手段：增大 $m$（实验显示误差对 $m$ 不敏感，说明默认值已够）、
  Adam 的自适应步长天然耐受噪声、避免二阶优化器。

  *问：NFF 与 sPI-GeM 都用生成模型，选哪个？*
  看空间维数与是否需要密度。空间维数 $<= 2$ 且下游要做
  贝叶斯推断（需要逐点密度）→ NFF；空间维数高、只需要采样
  与统计量 → sPI-GeM。两者的合流（学出来的基 + flow）
  是这条线公认的下一步。

  *问：分数阶阶数 $alpha$ 出现在采样密度里，反问题怎么办？*
  这是 MC-PINN 反问题最微妙的地方：$f_I (r)$ 依赖 $alpha$，
  $alpha$ 又是可训练参数。实践是把当前 $alpha$ 代入采样密度
  （重参数化技巧的思路），梯度通过采样过程用近似传递；
  论文的反问题实验（四参数识别到千分位）验证了这一处理。

  *问：这条线与直接用蒙特卡洛跑传统求解器相比，什么时候划算？*
  MC + 传统求解器：每个样本一次完整求解，样本贵但无训练成本；
  本线：训练贵、单次评估的代价接近于零。划算的分界是“要评估多少次”：
  一次性的单点预测用 MC；需要大量评估（优化、可靠性分析、
  实时数字孪生）时，训练一次的网络摊薄成本。
]

== F. 记号表补充（扩展章新增）

#table(
  columns: (auto, auto),
  inset: 6pt,
  align: (left, left),
  stroke: 0.4pt + luma(190),
  table.header([*记号*], [*含义*]),
  [$bold(B), bold(D)$], [NFF 参考场的低秩因子与对角噪声（snapshot 堆叠形式）],
  [$bold(w)$], [Woodbury 化简后的 $M$ 维中间量],
  [$bold(mu)_xi, bold(Sigma)_xi$], [隐变量 $bold(xi)$ 的后验均值与协方差],
  [$a_n, b_n$], [模态振幅与其平方（平流例子的解析解）],
  [$f_I, f_O, f_(I, t)$], [MC 估计的内段、外段、时间采样密度（Beta 分布）],
  [$xi, r$], [MC 估计的随机方向与随机半径],
  [$lambda_w, lambda_"IC"$ 等], [BO-fPINN 的动态损失权重],
  [$p$], [sPI-GeM 的基函数个数],
  [$cal(U)$], [sPI-GeM 输入函数的摘要表示（传感器 + PCA）],
  [$bold(Sigma)^s$], [第 $s$ 个 snapshot 的低秩加对角协方差],
)

== G. 六方法的超参数与训练配置速查

#ttable(
  ([配置], [NN-aPC], [NN-DO/BO], [NFF], [MC-PINN], [BO-fPINN], [sPI-GeM]),
  (
    ([优化器], [Adam], [Adam + L-BFGS-B], [Adam], [Adam], [Adam + L-BFGS-B], [Adam（GAN 双头）]),
    ([权重方案], [手工分组], [100x 约束 + 0.1x 正则], [三权 $w_"data/equ/bnd"$], [$w_"equ/g/u"$], [动态学习], [GAN 超参 $omega = 10$]),
    ([时域分段], [不适用], [10 段（$T = 10 pi$）], [不适用], [不适用], [5 段（$T = 5$）], [不适用]),
    ([关键超参], [aPC 阶数 $r$], [模态数 $N$], [KL 截断 $M$], [MC 样本 $m$、半径 $r_0$], [GL 阶数、网格 $Delta x$], [基个数 $p$、噪声 $M sigma_u^2$]),
    ([训练量级], [$10^5$ epoch], [$10^5$ epoch], [$10^3$ snapshot], [$10^4$ 到 $10^5$ 步], [迁移后 1/58 时间], [$10^4$ 步]),
  ),
  widths: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
)

读法：训练量级一栏是数量级参考（硬件相关），权重方案一栏
直接对应训练专题的病态讨论；时域分段是时间依赖实验的
通用配置，分段数由“单段训练末段误差开始爬升”触发。

== H. 反向问题的统一视角

六篇里五个做了反问题，做法惊人地一致，值得单独点出：
*把待识别量设为可训练变量，其余全部照旧*。

- NN-aPC：$k$ 的模态网络输出就是识别对象，三类数据（$u, k, f$
  传感器）同时进损失；
- NN-DO/BO：反应系数 $a, b$ 与网络参数联合优化，10 万 epoch 收敛；
- NFF：隐变量 $bold(xi)$ 的后验有闭式（低秩高斯回归），推断不用训练；
- MC-PINN：$alpha, gamma, c, v$ 四个物理参数设为可训练标量，
  采样密度随当前值更新；
- BO-fPINN：反应系数、噪声强度与分数阶阶数同时识别
  （三个初值跨越一倍以上也能收敛到 1% 内）。

统一的解释：反问题在“残差进损失”的框架里不需要新数学，
需要的只是把“哪些量可微”想清楚。NFF 是例外也是对照：
它的反问题走贝叶斯路线（后验分布而非点估计），
两类答案（点估计 + 残差、后验 + 区间）在方向 C 的指标体系里
各有对应的评价方式。

== I. 时间依赖方程的完整实验矩阵

时间依赖的三个实验（Burgers、扩散反应、分数阶反应扩散）
共用一套实验协议，把它们排成矩阵有助于看出各篇的证据链：

#ttable(
  ([检验项], [平流（1905）], [Burgers（1905）], [扩散反应（1905/2303）], [分数阶版本（2303）]),
  (
    ([精确解核对], [有（闭式 DO 解）], [无（MC 对照）], [无（MC/QMC 对照）], [无（QMC+GL 对照）]),
    ([特征值交叉], [不出现], [大量出现], [构造出现], [构造出现]),
    ([确定性初值], [是], [是], [噪声初值也测], [噪声初值也测]),
    ([随机维数], [1], [2], [19（KL）], [19（KL）+ 时间演化强迫]),
    ([反问题], [无], [无], [识别 $a, b$], [识别 $mu, epsilon, alpha$]),
    ([长时间], [无], [柄 $T = 10 pi$], [末端预测], [柄 $T = 5$]),
  ),
  widths: (auto, 1fr, 1fr, 1fr, 1fr),
)

矩阵的读法：每一行是一种“压力测试”，每一列是一个实验。
平流列最弱（一切都有真值所以最可信），Burgers 列补上交叉压力，
扩散反应列补上高维与噪声，分数阶列把全部压力叠在一起。
证据链是递进的：后面的实验假设前面的结构一致性已被核对。

== J. 术语索引补充（扩展章新增术语）

#ttable(
  ([术语], [定义位置], [一句话]),
  (
    ([损失权重病态], [训练专题], [各项梯度范数失衡，被最大项主导；处方是动态重加权]),
    ([谱偏差], [训练专题], [网络先学低频后学高频；特征扩展把高频变成线性权重]),
    ([变分形式 / 弱形式], [训练专题], [导数阶转移到测试函数；ReLU 可用、二阶导不再需要]),
    ([时域分段训练], [训练专题], [长时段切成子域逐段训练；接口由末端预测衔接]),
    ([特征函数解 / 均值衰减], [模态手推], [$E[cos xi t] = e^(-sigma^2 t^2 \/ 2)$：高斯特征函数的三角形式]),
    ([Woodbury 恒等式], [NFF 推导], [低秩加对角的逆与行列式降到 $M$ 维]),
    ([线性高斯后验], [NFF 推导], [观测越多后验越收缩；无观测退回先验]),
    ([稳定 Lévy 过程], [MC 推导], [分数阶 Laplacian 是它的生成元；重尾跳跃]),
    ([对称化], [MC 推导], [$xi$ 与 $-xi$ 平均消掉一阶项，二阶差商除 $r^2$ 有界]),
    ([DCA / 超梯度], [A 线第四篇], [差之差分解的交替迭代；TL1 的求解器]),
  ),
  widths: (auto, auto, 1fr),
)
