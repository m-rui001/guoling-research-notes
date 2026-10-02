#import "../../preamble.typ": *

= 实验画廊补充：BO-fPINN 与 sPI-GeM

== BO-fPINN：二维与噪声鲁棒性

主章节给了特征值交叉与迁移学习两个亮点。这里补二维结果与
噪声初值实验，它们检验的是“组装是否破坏了各组件的鲁棒性”。

#figure(
  align(center, image("../fig_b/x_b7_noise.png", width: 82%)),
  caption: [19 维随机分数阶反应扩散、噪声初值（30 个传感器加
  $sigma = 0.1$ 噪声）下的均值与方差预测（误差棒为参考解的波动）。
  模态分解 + GL 离散的组合没有放大观测噪声。],
)

#figure(
  align(center, image("../fig_b/x_b7_2dvar.png", width: 82%)),
  caption: [二维双 Riesz 随机分数阶反应扩散的方差场
  （$alpha = 1.8$，随机系数 + KL 型随机强迫，QMC+GL 参考解）：
  方差的空间结构与参考解一致，长时间积分按时域分段训练。],
)

== sPI-GeM：从随机场到随机方程

#figure(
  align(center, grid(columns: 2, gutter: 10pt,
    image("../fig_b/x_b8_helm.png"), image("../fig_b/x_b8_darcy.png"))),
  caption: [左：52 维随机 Helmholtz 方程的解均值场
  （$D_zeta = 52$，$D_x = 2$）：生成模型学 52 维系数分布，
  基网络学空间结构。右：2 维 Darcy 流：渗透率 $lambda$ 与解 $u$
  各配一组基网络，$u$ 均值误差 1.01%。],
)

#figure(
  align(center, image("../fig_b/x_b8_invmean.png", width: 78%)),
  caption: [一维随机椭圆反问题：系数网络输出拆成 $u$ 与 $lambda$
  两半（共享潜表示），同时重构解与渗透率场。
  $lambda$ 均值误差 0.56%、标准差误差 2.40%。],
)

== 附录补充：B 线新增文献

- I. Lagaris, A. Likas, D. Fotiadis. *Artificial neural networks for
  solving ordinary and partial differential equations.* IEEE Trans.
  Neural Netw., 1998.
- W. E, B. Yu. *The deep Ritz method.* 2018.
- J. Han, A. Jentzen, W. E. *Solving high-dimensional partial
  differential equations using deep learning.* PNAS, 2018.
- L. Kharazmi, Z. Zhang, G. E. Karniadakis. *vpINNs.* 2020.
- L. Pang, K. Lu, G. E. Karniadakis. *fPINNs.* SIAM J. Sci. Comput.,
  2019.
- Z. Chen, J. Chen. *Universal approximation to nonlinear operators.*
  1993 / 1995.
- L. Lu, P. Jin, G. E. Karniadakis. *DeepONet.* Nature Mach. Intell.,
  2021.
- Z. Li, et al. *Fourier neural operator.* ICLR, 2021.
- P. Yang, et al. *B-PINNs: Bayesian physics-informed neural networks.*
  J. Comput. Phys., 2021.
- S. Wang, Y. Teng, P. Perdikaris. *Understanding and mitigating
  gradient flow pathologies in PINNs.* SIAM J. Sci. Comput., 2021.
- T. Sapsis, P. Lermusiaux. *Dynamically orthogonal field equations.*
  Physica D, 2009.
- M. Choi, et al. *DO–BO equivalence.* 2014.
- L. Dinh, et al. *NICE / RealNVP.* 2015 / 2017.
- D. P. Kingma, P. Dhariwal. *Glow.* 2018.
- I. Goodfellow, et al. *GANs.* 2014; M. Arjovsky, et al. *WGAN.* 2017;
  I. Gulrajani, et al. *WGAN-GP.* 2017.
- K. Sohl-Dickstein, et al. 2015; J. Ho, et al. *DDPM.* 2020.
- P. Yang, et al. *PI-GAN.* SIAM J. Sci. Comput., 2020.
- F. Zhong, et al. *PI-VAE.* 2023.

== 自检清单扩展

#boxed(title: "读完扩展章后追加回答")[
  11. 损失权重病态怎么诊断？动态加权的判据是什么？（训练专题）
  12. 谱偏差为什么对随机方程影响大？特征扩展为什么有效？（训练专题）
  13. 弱形式损失换来什么、放弃什么？哪些方程必须走弱形式？（训练专题）
  14. 平流方程的 DO 分解里，振幅为什么从零开始增长？
      $a_1^2$ 的显式表达式是什么？（模态分解手推）
  15. DO 与 BO 等价性的直觉解释是什么？NN 版本里它表现为
      什么现象？（模态分解手推）
  16. Woodbury 恒等式把 $N times N$ 的求逆降到什么复杂度？
      哪两个矩阵恒等式参与？（NFF 推导）
  17. 后验协方差 $bold(Sigma)_xi$ 在没有观测时退化成什么？
      这对应什么行为？（NFF 推导）
  18. 分数阶 Laplacian 的傅里叶符号是什么？$alpha -> 2$ 与
      $alpha -> 0$ 各退化成什么？（MC 推导）
  19. 内段公式的对称化消掉了什么？为什么被估量除以 $r^2$
      之后有界？（MC 推导）
  20. 六方法对比表里，哪两列互相排斥（不可兼得）？
      选型决策树的第一个分叉是什么？（对比与选型）
]
