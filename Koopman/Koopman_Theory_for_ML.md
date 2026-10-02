# Koopman Operator Theory for ML People
### A technical study note (no ergodic theory required)

**Deliverables of this note.** (i) the 1931 idea and what it costs; (ii) eigenfunctions, eigenvalues, modes, spectral expansion; (iii) spectrum ↔ dynamics dictionary; (iv) DMD / EDMD as regression; (v) Takeishi et al.'s learned observables; (vi) Koopman-invariant subspaces and the global obstruction that bounds all deep Koopman learning.

---

## 1. Koopman (1931): moving a nonlinear problem into a linear one

Take a dynamical system on a state space $M$, either discrete, $x_{t+1}=F(x_t)$, or a flow $\dot x=f(x)$, $x(t)=S^t x_0$. An **observable** is any scalar (or vector) function $g:M\to\mathbb R$ — a measurement. Koopman's move is to stop tracking the point $x$ and start tracking the *function*, via the composition operator

$$(Kg)(x)=g(F(x))=(g\circ F)(x),\qquad (U^t g)(x)=g(S^t x).$$

$K$ is **linear in $g$** no matter how nonlinear $F$ is: $K(a g_1+b g_2)=a\,Kg_1+b\,Kg_2$ — the nonlinearity is absorbed into the features, like a kernel trick with an infinite feature list of *functions*. For continuous time, $U^t=e^{tL}$ with generator $Lg=\nabla g\cdot f$; for Hamiltonian flows preserving a measure $\mu$, $U^t$ is **unitary** on $L^2(M,\mu)$ (Koopman's original setting, where $L$ is self-adjoint).

**What is gained:** superposition, semigroup/matrix-calculus tooling, and — crucially — *spectral* reasoning about a nonlinear system: prediction, mode extraction, stability and linear-optimal control all become statements about eigenvalues of $K$.

**What is paid:** (a) $K$ is generically **infinite-dimensional**; (b) it is **non-normal** for dissipative systems, so eigenvalues alone do not control transient behavior; (c) its spectrum need not consist of eigenvalues at all — for chaotic/mixing dynamics (and even for the pendulum, Koopman & von Neumann 1932) there is **continuous spectrum**; the "eigenvectors" may be generalized eigenfunctions outside $L^2$. Every practical method below is a bet that a low-dimensional piece of this operator can be truncated without too much loss.

## 2. Eigenfunctions, eigenvalues, Koopman modes

A **Koopman eigenfunction** is $\varphi:M\to\mathbb C$ with

$$K\varphi=\lambda\varphi\iff \varphi(F(x))=\lambda\varphi(x);\qquad L\varphi=\mu\varphi \Rightarrow \varphi(S^t x)=e^{\mu t}\varphi(x),\ \lambda=e^{\mu\Delta t}.$$

Products of eigenfunctions are eigenfunctions ($\lambda_i\lambda_j$); level sets of eigenfunctions are invariant sets (isochrons, isostables); eigenfunctions are *coordinates that linearize*. They are defined only up to scale — the same gauge freedom DMD modes inherit.

Choosing the **full-state observable** $g(x)=x$ gives the **Koopman mode decomposition** (Mezić 2005): if $x=\sum_j\varphi_j(x)v_j$ then

$$x(t)=\sum_j e^{\mu_j t}\,\varphi_j(x_0)\,v_j,\qquad\text{discrete: } x_t=\sum_j \lambda_j^{\,t}\,\varphi_j(x_0)\,v_j .$$

The **Koopman modes** $v_j\in\mathbb C^n$ are the spatial patterns ("how the mode looks in state space"), $\varphi_j(x_0)$ is the initial amplitude, and $\lambda_j$ gives frequency $\arg\lambda_j/\Delta t$ and growth $|\lambda_j|$. Modes depend on the observable ($y=Cx$ gives modes $Cv_j$); eigenvalues and eigenfunctions do not. If the discrete spectrum has principal eigenvalues $\lambda_1,\dots,\lambda_n$, it is **lattice-like**, $\{\prod_i\lambda_i^{k_i}\}$ (Mezić 2020) — a useful sanity check on estimated spectra.

## 3. Reading the spectrum: $|\lambda|\gtrless 1$ and the unit circle

- $|\lambda|=1$ (continuous time $\mathrm{Re}\,\mu=0$): **persistent/oscillatory**. For a measure-preserving system the point spectrum lies *on* the unit circle (unitarity), and the eigenvalues form a multiplicative subgroup of it; $\lambda=1$ is the constant eigenfunction, i.e. the invariant measure. Ergodicity $\iff$ $\lambda=1$ is a *simple* eigenvalue; weak mixing $\iff$ $1$ is the **only** eigenvalue. Mixing implies weak mixing but is *not* implied by it (mixing is a decay-of-correlations statement, not merely a spectral-gap statement). So a purely continuous spectrum buys you no finite eigenfunction basis at all — the reason low-dimensional Koopman models degrade on weakly mixing/chaotic dynamics.
- $|\lambda|<1$ ($\mathrm{Re}\,\mu<0$): **decaying transients**. Their level sets are isochrons/isostables; the slowest decaying modes govern late-time behavior.
- $|\lambda|>1$ ($\mathrm{Re}\,\mu>0$): **unstable/growing**. Repelling directions, off-attractor transients, unbounded domains — or numerical artifacts.

**Why near-unit-circle eigenvalues drive prediction:** the "slow subspace" spanned by eigenfunctions with $|\lambda|$ close to $1$ captures long-term behavior once fast transients die; when a **spectral gap** separates slow from fast, a low-dimensional truncation is accurate and the truncation error after $t$ steps scales like $|\lambda_{\text{next}}|^t$. Two warnings from the same literature: for non-normal $K$ (or Jordan blocks $\Rightarrow$ generalized eigenfunctions $t^k e^{\mu t}$) transient growth can dwarf the eigenvalue asymptotics; and finite-dimensional regression estimators produce **spurious eigenvalues** that are eigenvalues of the estimator, not of $K$ (Korda & Mezić 2018 study when these converge).

## 4. DMD and EDMD: Koopman estimation as linear regression

Pick a dictionary $\Psi(x)=[\psi_1(x),\dots,\psi_K(x)]^\top$ (polynomials, RBFs, Fourier modes, random features, a network) and data pairs $(x_m,\;y_m=F(x_m))_{m=1}^M$; collect $\hat X=[\Psi(x_1)\cdots\Psi(x_M)]$, $\hat Y=[\Psi(y_1)\cdots\Psi(y_M)]$ and solve the least-squares problem $\min_K\|\hat Y-K\hat X\|_F^2$:

$$\boxed{K=\hat Y\hat X^{+}=A\,G^{+}},\qquad G=\tfrac1M\hat X\hat X^{*},\quad A=\tfrac1M\hat Y\hat X^{*},$$

writing $\Psi(y)\approx K\Psi(x)$. (Williams et al. 2015 use the row-vector convention $\Psi(x)\,K$ with $\Psi$ a row, giving $K=G^{+}A$; the two differ by transposition, i.e. right vs. left eigenvectors. Use the pseudoinverse because $G$ is often rank-deficient/redundant.) Then: **eigenvalues of $K$ estimate Koopman eigenvalues; right eigenvectors $\xi_j$ give eigenfunctions $\varphi_j(x)=\xi_j^{\top}\Psi(x)$; left eigenvectors give modes $v_j=(w_jB)^{\top}$ for an observable $g(x)=B\Psi(x)$** (for $g(x)=x$, $v_j$ is the Koopman mode). As $M\to\infty$, $K$ converges to the **Galerkin projection** of $\mathcal K$ onto $\mathrm{span}(\Psi)$: the residual is $L^2(\mu)$-orthogonal to the dictionary. With $\Psi(x)=x$ (identity dictionary) $K$ *is* the classical **DMD** matrix $A_{\text{DMD}}=YX^{+}$ — so DMD is the case "linear observables only": it is exact when the state components themselves span an invariant subspace, and provably too restrictive otherwise. DMD yields eigenvalues and modes but not eigenfunctions; EDMD yields all three. EDMD applies unchanged to Markov data, where it estimates eigenfunctions of the Kolmogorov backward operator, the "stochastic Koopman operator" (Williams et al. 2015).

## 5. Learning the observables: Takeishi et al. 2017 (LKIS)

The remaining freedom is $\Psi$. Takeishi, Kawahara & Yairi (NeurIPS 2017) make $\Psi=g_\theta$ a neural encoder and learn $\theta$ so that a linear model *closes* in the learned coordinates. With $Y_0=[g(x_0)\cdots g(x_{m-1})]$, $Y_1=[g(f(x_0))\cdots]$ and $A=Y_1Y_0^{+}$ the minimum-norm least-squares solution, they minimize the **residual sum of squares of the regression** plus a reconstruction term:

$$\mathcal L=\underbrace{\big\|Y_1-AY_0\big\|_F^2}_{\mathcal L_{\text{RSS}}}+\;\gamma\underbrace{\sum_j\|y_j-h(g(\tilde x_j))\|^2}_{\mathcal L_{\text{rec}}}.$$

Their Theorem 1: $\mathcal L_{\text{RSS}}=0$ **iff** $\{g_i\}$ spans a Koopman-invariant subspace — so the loss is a differentiable, data-only surrogate for invariance. The reconstruction term with a decoder $h$ (both MLPs) prevents trivial solutions ($g\equiv$ const). Partial observations are handled by a *learned* linear delay embedding $g(x)\leftarrow$ `Linear(delay-embed(y))`, so encoder, embedder and invariant subspace are learned jointly; afterwards ordinary DMD on the learned coordinates yields eigenvalues, modes $w$, and eigenfunction values $\varphi=z^{H}g(x)$ ($z$ = left eigenvectors of $A$). Practical notes they flag: the RSS loss is non-decomposable over samples (mini-batch SGD has no clean guarantee), and the method is markedly more robust to observation noise than EDMD with a fixed polynomial dictionary. This is the direct ancestor of deep Koopman autoencoders, consistent-Koopman nets, and (with the added device of an auxiliary network parameterizing a **continuous** eigenvalue) Lusch, Kutz & Brunton (2018).

## 6. Koopman-invariant subspaces, and the wall you always hit

A subspace $\mathcal G$ of observables is **Koopman-invariant** if $g\in\mathcal G\Rightarrow Kg\in\mathcal G$. If $\mathcal G$ is finite-dimensional and spanned by $\Psi$, then $K$ restricted to $\mathcal G$ is a $K\times K$ matrix and one obtains an **exact** finite-dimensional linear model, no regression error. Eigenfunctions span 1-dimensional invariant subspaces; any set of eigenfunctions spans one.

Brunton, Brunton, Proctor & Kutz (2016) ask for the useful special case: an invariant subspace **containing the state** ($x_1,\dots,x_n$ among the observables), which is what control (e.g. LQR) requires. Their negative result: this is impossible for **any system with multiple fixed points, or more general attractors**. The argument is short enough to reproduce: if the state is in the span, then $\dot x=f(x)$ must be in the span too, so the flow is exactly linear in these coordinates, $x(t)=Ce^{Kt}$-type dynamics; that is a **topological conjugacy to a finite-dimensional linear system**. But an invertible linear system has a single fixed point, and conjugacy preserves the number of fixed points. Hence finite-dimensional, state-including, Koopman-invariant subspaces exist essentially only for systems with a single isolated fixed point (the paper's constructive examples are Carleman-like: e.g. $\dot x_1=\mu x_1,\ \dot x_2=\lambda(x_2-x_1^2)$ is exactly linear on $\mathrm{span}\{x_1,x_2,x_1^2\}$). Systems with multiple basins still have eigenfunctions, but they are **discontinuous**, supported on one basin, hence not representable by a smooth finite dictionary; EDMD must be run per-basin rather than globally (Williams et al. 2015).

**Consequences for deep Koopman learning.** (1) Learning $\Psi$ that makes dynamics linear *is* learning a conjugacy; topology caps what it can do — expect basin-local or attractor-local accuracy, not a global finite linearization. (2) Linearity/prediction loss alone is satisfiable by degenerate encoders; identifiability requires the reconstruction/invertibility term. (3) $|\lambda|$ near $1$ modes are both the most valuable (long-horizon error) and the hardest to estimate, and continuous-spectrum systems admit no finite eigenfunction basis at all. (4) Estimated spectra require convergence-aware handling: spurious eigenvalues, non-normality/transient growth, and mode-normalization gauge freedom are properties of the *estimator*, not of $\mathcal K$.

---

### Chinese TL;DR
Koopman 1931 把非线性动力学 $x\mapsto F(x)$ 提升为观测函数空间上的**线性**算子 $(Kg)(x)=g(F(x))$，代价是无穷维、非正规性、以及可能有连续谱。谱分解给出 $x_t=\sum_j\lambda_j^t\varphi_j(x_0)v_j$：$|\lambda|=1$ 为持续/振荡模（对保测系统点谱全在单位圆上），$|\lambda|<1$ 为衰减瞬态，$|\lambda|>1$ 为不稳定方向；预测能力主要由单位圆附近的**慢子空间**与谱间隙决定。EDMD 用字典 $\Psi$ 做最小二乘回归 $K=\hat Y\hat X^+ = AG^+$，DMD 是 $\Psi(x)=x$ 的特例；Takeishi et al. 2017 用神经网络学习 $\Psi$，以回归残差平方（+重构项）为损失，RSS=0 当且仅当学到的观测张成 Koopman 不变子空间。Brunton et al. 2016 的根本限制：包含状态变量的有限维不变子空间**只**可能存在于单不动点系统，因为有限维线性系统只有一个不动点，而拓扑共轭保持不动点个数——这正是深度 Koopman 模型难以全局线性化的原因。

### Sources
- B. O. Koopman, *Hamiltonian systems and transformation in Hilbert space*, PNAS 17:315–318, 1931; Koopman & von Neumann, PNAS 18:255–263, 1932 (continuous spectrum).
- I. Mezić, *Spectral properties of dynamical systems, model reduction and decompositions*, Nonlinear Dyn. 41:309–325, 2005; Mezić, *Spectrum of the Koopman operator, spectral expansions in functional spaces, and state-space geometry*, J. Nonlinear Sci. 30:2091–2148, 2020 — https://arxiv.org/abs/1702.07597.
- M. O. Williams, I. G. Kevrekidis, C. W. Rowley, *A data-driven approximation of the Koopman operator: extending DMD*, J. Nonlinear Sci. 25(6):1307–1346, 2015 — https://arxiv.org/abs/1408.4408 (PDF also at cwrowley.princeton.edu/papers/Williams-edmd-2015.pdf).
- N. Takeishi, Y. Kawahara, T. Yairi, *Learning Koopman invariant subspaces for DMD*, NeurIPS 30:1130–1140, 2017 — https://arxiv.org/abs/1710.04340.
- S. L. Brunton, B. W. Brunton, J. L. Proctor, J. N. Kutz, *Koopman invariant subspaces and finite linear representations of nonlinear dynamical systems for control*, PLoS ONE 11(2):e0150171, 2016 — https://arxiv.org/abs/1510.03007. *(Local copy: `../Koopman.pdf`.)*
- B. Lusch, J. N. Kutz, S. L. Brunton, *Deep learning for universal linear embeddings of nonlinear dynamics*, Nat. Commun. 9:4950, 2018 — https://arxiv.org/abs/1712.09707.
- M. Korda, I. Mezić, *On convergence of extended dynamic mode decomposition to the Koopman operator*, J. Nonlinear Sci. 28:687–710, 2018.

**See also:** `koopman_autoencoder_lineage.md` (same folder) — the architecture-by-architecture lineage of deep Koopman autoencoders, which picks up exactly where §5–§6 end.
