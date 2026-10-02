# The Deep Koopman Autoencoder Lineage — Technical Study Note

Scope: background for *Information Shapes Koopman Representation* (Cheng et al., ICLR 2026,
arXiv:2510.13025). Each entry: problem → architecture/objective → limitations that motivate an
information-theoretic regularizer.

---

## 1. Otto & Rowley 2019 — Linearly-Recurrent Autoencoder Networks (LRAN/LRA)
*SIAM J. Appl. Dyn. Syst.* 18(1):558–593; arXiv:1712.01378.

**Problem.** EDMD/KDMD approximate the Koopman operator on a *user-chosen* dictionary. A rich
dictionary overfits: their Example 2.1 shows that adding one superfluous observable to a
2-observable linear system (empirical residual < 6e-15) yields spurious eigenvalues
{1.0413, 0, 0.3515}. EDMD also presumes linear reconstruction of the state from Koopman modes,
i.e. that the full-state observable lies in the learned invariant subspace.

**Architecture.** An *undercomplete* autoencoder: encoder `Ψ(x;θ_enc) = z ∈ R^d` (d ≪ n), a single
learned matrix `K(θ_K)` advancing the latent, `ẑ_{t+τ} = [K^τ]ᵀ z_t`, and a **nonlinear** decoder
`Ψ̃(ẑ;θ_dec) → x̂`. The nonlinear decoder relaxes the "full state is in the subspace" assumption and
takes the place of Koopman modes. **Why linear:** the Koopman operator is linear by definition, so
the finite-rank approximation on the learned dictionary is one shared matrix `K`. This buys
`K = W_R Λ W_L*`, hence Koopman eigenvalues/eigenfunctions/modes exactly as in EDMD, and τ-step
rollout as matrix powers. Note the *bottleneck is the small dimension d*; only the latent **dynamics**
is linear — the encoder/decoder are deep ELU networks.

**Objective (Eq. 25).** Multi-step consistency:
`J = E[1/(1+β) ( Σ_{τ=0}^{T-1} (δ^τ/N₁(δ)) ‖x̂_{t+τ}−x_{t+τ}‖²/(‖x_{t+τ}‖²+ε₁)
+ β Σ_{τ=1}^{T-1} (δ^{τ-1}/N₂(δ)) ‖ẑ_{t+τ}−z_{t+τ}‖²/(‖z_{t+τ}‖²+ε₂) )] + Ω`,
with `N₁(δ)=Σ_{τ=0}^{T-1}δ^τ`, `N₂(δ)=Σ_{τ=1}^{T-1}δ^{τ-1}`. Two terms: decoded multi-step
reconstruction error, and a **latent multi-step consistency** error comparing the rolled-forward
latent `K^τ z_t` against the encoder of the true future state `Ψ(x_{t+τ})`. δ∈(0,1] discounts later
steps (prioritises short-term), β trades the two terms, Ω is optional regularization. **Rationale:**
in non-normal systems (shear flows) low-amplitude but "highly observable" features carry energy into
larger modes several steps later; single-step (snapshot-pair) losses discard them, multi-step training
forces them into the dictionary. Paper generalises to continuous time / irregular sampling with
weights ρ_i(t_k).

**Limitations.** No stability constraint on K → slightly unstable modes (|λ|>1), tolerable
short-term, fatal long-term. d, δ, β, T are unprincipled hyperparameters. Nothing prevents latent
redundancy/collapse: all d modes are trained although only a few are genuine. No uncertainty.

---

## 2. Azencot, Erichson, Lin & Mahoney 2020 — Consistent Koopman Autoencoders (CKA)
ICML 2020; arXiv:2003.02236.

**Problem.** RNNs ignore physics; existing Koopman AEs model only forward dynamics, and their C is
useless for backward prediction. Goal: exploit forward *and* backward dynamics.

**Architecture.** Encoder χ_e, decoder χ_d, and two *bias-free linear layers* C (forward) and
D (backward): `f̂_{k+1} = χ_d∘C∘χ_e(f_k)`, `f̌_{k-1} = χ_d∘D∘χ_e(f_k)`. All connections
bidirectional. Loss: `E = λ_id E_id + λ_fwd E_fwd + λ_bwd E_bwd + λ_con E_con`, with multi-step
(λ_s = 8) forward/backward prediction in *observation* space.

**The consistency loss (Eq. 12).** With `D_{k*}` = top k rows of D and `C_{*k}` = leftmost k columns
of C:
`E_con = Σ_{k=1}^{κ} [ 1/(2k)‖D_{k*}C_{*k} − I_k‖_F² + 1/(2k)‖C_{k*}D_{*k} − I_k‖_F² ]`.
**Why.** `E_bwd` only touches D — it cannot change C. Prop. 1 (continuous case) shows invertibility
⟺ `⟨ξ_i, U_ψ K_φ ξ_j⟩ = δ_ij`, i.e. `DC = I_κ`. Prop. 2 shows that in the *discrete* sampling
setting the correct condition is the nested leading-principal-submatrix sum, hence the graded sum
over k rather than plain `‖DC − I‖_F²` (Eq. 12 approximates the true m-dimensional condition, which
is only recoverable when κ = m). Effect: couples C and D, and *quasi-stability* — eigenvalues are
pushed toward the unit circle, relaxing strict stability. Empirically the regularizer, not the
backward term, drives robustness to noise and to long horizons.

**Limitations.** Assumes the backwards map exists (diffusive systems fail — their own negative
result: adding diffusion to Euler flow degrades CKA). Consistency is enforced per leading k×k block
in latent space only. λ weights are hand-tuned; nothing controls latent content. Treats only
time-reversible/volume-preserving structure.

---

## 3. Pan & Duraisamy 2020 — Physics-informed probabilistic linear embeddings
*SIAM J. Appl. Dyn. Syst.* 19(1):480–509; arXiv:1906.03663.

**Problem.** Unify (i) *continuous-time* Koopman learning, (ii) provable stability, (iii) uncertainty
quantification — all absent from LRAN/deep-DMD (their Table 1 lists these as exclusive to this work).

**Formulation.** Measure-theoretic. Koopman error functional
`J[Φ] = max_ψ min_{f∈F_D} ‖f − Kψ‖²_F` with reconstruction `R[Φ,Ψ] = ‖Ψ∘Φ − I‖²`. Two variants:
*differential* (`F·∇Φ = ΦK`, needs governing equations) and *recurrent*
(`‖Φ e^{tK} − K_t Φ‖²_G`, needs trajectories).

**Guaranteed stability (structural parameterisation).** K is constrained to a tridiagonal
`K_stable` with diagonal `−σ_i² ≤ 0` and skew-symmetric off-diagonals `±ζ_i` (Eq. 31). Every
realisation is stable (eigenvalues in the closed left half-plane), and Theorem 1 states any real
diagonalisable matrix with non-positive-real-part eigenvalues is similar over R to some such
`K_stable` — so **no loss of expressivity**.

**Architecture.** SVD-DMD is embedded in the autoencoder and only the *residual* is learned:
observables = nonlinear NN observables + linear SVD-DMD observables `Φ_svd(z) = zΛV_D` (dominant
POD modes); decoder = nonlinear + linear reconstruction. K is initialised at the DMD estimate, so for
essentially linear dynamics optimisation starts near the desired minimum. Motivation: POD/DMD is a
robust approximation of the Koopman subspace, and pure-nonlinear parameterisations admit spurious
"over-complicated" eigenfunctions (e.g. `φ_j = y_j^{2n_j+1}` ⇒ eigenvalues `(2n_j+1)λ_j`) under the
non-uniqueness of the ideal Φ.

**Probabilistic.** Hierarchical Bayesian; mean-field variational inference with diagonal Gaussians
via ADVI (Edward/TensorFlow); maximise the ELBO. They note prior deep-Koopman weight decay is
equivalent to MAP estimation of the same posterior.

**Limitations.** Assumes eigenvalues have non-positive real parts — excludes genuine instability and
transient growth. Explicitly **does not handle continuous spectra**, which is precisely what chaotic
systems have. Mean-field Gaussian posterior is crude/multimodal-blind. Single-trajectory,
single-attractor regime.

---

## 4. Pan et al. 2023 — PyKoopman (the KAE baseline)
arXiv:2306.12962; github.com/dynamicslab/pykoopman.

**Problem.** Standardisation, not new theory: no common open-source Koopman implementation existed.
sklearn-style `Koopman` estimator = `observables` (Identity, Polynomial, TimeDelay,
RadialBasisFunctions, RandomFourierFeatures, CustomObservables, ConcatObservables) × `regressor`
(DMD, DMDc, EDMD, EDMDc, KDMD, HDMD, HDMDc, **NNDMD**), plus differentiation and a sparsity-based
analytics module. Flagged extensions: bilinearisation, continuous spectrum (resDMD, MP-EDMD, π-DMD).

**The KAE used as a baseline** in Cheng et al. 2026 is `NNDMD` wrapping `DLKoopmanRegressor`
(`pykoopman/regression/_nndmd.py`): FFNN encoder (e.g. tanh, [32,32], output_size d = 6) → Koopman
propagator → FFNN decoder (linear activations) → state. The propagator is a *continuous-time
generator* passed through `matrix_exp(dt·K)`. Three parameterisations: `StandardKoopmanOperator`
(unconstrained full K — default), `HamiltonianKoopmanOperator` (`K = A − Aᵀ`, purely imaginary
spectrum), `DissipativeKoopmanOperator` (`K = diag(−σ²) + (A − Aᵀ)`).

**Objective.** `MaskedMSELoss(look_forward)` multi-step rollout MSE on *decoded states*, plus
reconstruction MSE at x and along the predicted sequence. There is **no explicit latent consistency
term, no forward/backward consistency, no entropy/KL term** — latent structure is only implicit in
the rollout. Practical advice: recurrent (multi-step) loss > one-step; L-BFGS faster than Adam but
can diverge; observable choice should scale with dataset size.

**Limitations as a baseline.** The default operator is unconstrained and unregularized → eigenvalue
drift. d is a hyperparameter; no criterion for how many or which observables to keep (the analytics
module prunes *post hoc*). Continuous spectra remain unsolved by their own admission. Latent is fit
purely by self-supervised reconstruction + rollout, so mode dominance/collapse is uncontrolled — the
target paper reports KAE's learned manifold collapsing to nearly one dimension.

---

## 5. Li, He, Wu, Katabi & Torralba 2020 — Compositional Koopman Operators (CKO)
ICLR 2020; arXiv:1910.08264.

**Problem.** Deep-Koopman methods assume a *fixed-dimensional* state space, so they fail when the
number of objects varies. Naively, K has N²m² parameters, which harms sample efficiency.

**Architecture.** A GNN encoder φ (Interaction-Network style) produces object-centric Koopman
embeddings: relation function `f_R(o_u, o_v, a^r)` → edge effect `e_k`; object function
`f_O(o_i, Σ_{k∈N_i} e_k)` → `g_i^t ∈ R^m`; concatenation `g^t ∈ R^{Nm}`. Dynamics in the embedding is
block-wise linear: `g^{t+1} = K g^t + L u^t` with **`K = σ ⊗ K̂`**, where σ ∈ R^{N×N×h} is a one-hot
index tensor assigning block type and `K̂ ∈ R^{h×m×m}` holds h shared sub-matrices. So *per-object*
("self", diagonal) and *per-relation* ("interaction", off-diagonal) operators are shared across
object pairs of the same type — the count drops from N²m² to `h m²`, independent of N. Decoder ψ
maps back to state.

**Identification & objective.** `min_{K̂,L̂} ‖(σ⊗K̂)g^{1:T-1} + (σ⊗L̂)ũ − g^{2:T}‖₂` — linear least
squares with a closed-form solution, hence cheap online adaptation to unseen physical parameters
(≈8 episodes / 800 samples). Training loss
`L = L_ae + λ₁ L_pred + λ₂ L_metric`, with `L_pred` a T-step rollout decoded to state, and
`L_metric = Σ_ij | ‖g_i−g_j‖ − ‖x_i−x_j‖ |` so that latent-space distances are meaningful and the
control cost can be defined there. **Control:** linear latent dynamics turn the reaching task into a
QP — minimise `Σ_t [1[t=T]‖g^T − g*‖² + λ‖u^t‖²]` s.t. `g^{t+1} = Kg^t + Lu^t`, plus MPC with
periodic feedback. Experiments: Rope (5–9 masses, extrapolating to 10–14), Soft, Swim. Ablations:
Block > None > Diag; m = 32 (best sim at m = 16), λ₁ = 1, λ₂ = 0.3.

**Limitations.** The object/relation vocabulary (`a^o`, `a^r`) and the block-type assignment σ must
be **specified by hand** from the graph structure. Homogeneous, permutation-symmetric interactions
only — no object-specific or history-dependent couplings. K and L come from unconstrained least
squares, so no spectral/stability guarantee. `L_metric` is a heuristic isometry prior. h, m, and the
metric weight are hyperparameters; full observability of all object states is assumed.

---

## 6. Cheng et al. 2025 — Learning Chaos in a Linear Way (PFNN)
ICLR 2025; arXiv:2503.14702.

**Problem.** Dissipative chaos has positive Lyapunov exponents, so stepwise-MSE autoregressive models
accumulate exponentially growing error; long-term *statistics* (the invariant measure on the
attractor) matter more than trajectory accuracy. Prior work aligns invariant measures but needs the
invariant measure, which is generally unknown.

**Key idea.** Split dissipative chaos into two physical phases and enforce each by *structure* with a
simple one-step loss. (i) **Finite-rank approximation:** `G_L = Π_L G Π_L`, `G_ij = ⟨φ_i, Gφ_j⟩`;
Theorem 3.1 gives pointwise convergence; base loss = prediction + γ·reconstruction with bijectivity
`Id = g^de∘g^en`. (ii) **Contraction phase:** attractors imply `‖G_c‖ ≤ 1`; enforce
`⟨G_cφ, G_cφ⟩ ≤ λ²‖φ‖²`, 0<λ≤1, implemented as
`Γ_con = ReLU(σ(Ĝ_cᵀĜ_c − I_d))`, i.e. eigenvalues of `Ĝ_cᵀĜ_c` ≤ 1. The relaxation time — bounded by
the log-Sobolev time, `k ∝ ⌈(1/c_LSI) log(‖z‖₂/ε)⌉` — splits the dataset. (iii) **Measure-invariant
phase:** a measure-preserving T makes the Koopman operator unitary on L² (Koopman–von Neumann ergodic
theorem), with spectrum on the unit circle. Loss (Eq. 8) = forward loss
`‖Ĝ_m g(z_k) − g(z_{k+1})‖₂` + backward loss `‖Ĝ_m* g(z_{k+1}) − g(z_k)‖₂` +
`γ₂(‖Ĝ_m*Ĝ_m − I_d‖_F + ‖Ĝ_mĜ_m* − I_d‖_F)` (consistent constraint) + γΓ_rec. So a bi-directional
unitary operator is learned **without knowing the invariant measure**; learned eigenvalues cluster
densely on the unit circle (continuous spectrum) and the energy spectrum matches truth. Beats LSTM,
Koopman/PyKoopman, FNO and MNO on KLD/MMD/TKE over Lorenz 63/96 (m=9/40/80), KS (m=128) and
Kolmogorov flow; γ₂ ∈ [0.1, 0.5] is the sweet spot.

**Limitations.** Requires the dissipative-chaos-with-global-attractor assumption, and the relaxation
time k must be estimated (via a log-Sobolev constant that is itself inferred). Unitarity is the right
structure *only* on the attractor of a measure-preserving system — Stage I (moduli ≤ 1) and Stage II
(moduli = 1) conflict for forced/controlled/input-driven dynamics. The hard two-way data split and
the two operators add hyperparameters. Continuous spectrum is only *diagnosed*, not parameterised.
Above all, every constraint is **structural/spectral** (eigenvalue geometry) and evaluated
*externally* by KLD/MMD/TKE: nothing in training controls how much information the latent retains,
how modes are distributed, or collapse/redundancy.

---

## Synthesis: the common gap

All six share the same latent objective — reconstruction + multi-step prediction, self-supervised —
so the latent representation is **underdetermined**. Their only regularizers are domain-specific
priors introduced one system class at a time: stability (Pan & Duraisamy), invertibility/consistency
(Azencot), dissipation + unitarity (PFNN), permutation/relation symmetry (CKO), sparsity
(PyKoopman analytics). Every trade-off is mediated by hand-tuned scalars (β, δ, γ, γ₁, γ₂, λ_bwd,
λ_con, λ₁, λ₂, d, m, h, L, k) with no transferable principle for choosing latent dimension or
balancing latent *simplicity/linearity* against *expressiveness*. The empirically observed failure
modes — mode collapse, dominance of a few modes, distorted low-dimensional manifolds — are exactly
the information-bottleneck symptoms (compression vs. predictive relevance; marginal vs. von Neumann
entropy) that Cheng et al. 2026 (arXiv:2510.13025) formalise into an information-theoretic Lagrangian.
