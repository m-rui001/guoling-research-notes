# Information-Bottleneck / Mutual-Information Toolkit for an IB-Koopman Paper
Study note. Notation: `I(A;B)` = mutual information (base chosen per line, stated), `h(·)` = differential entropy, `H(·)` = Shannon entropy, `D_KL` = KL divergence, `ln` natural log unless a base is written.

Standing Markov-chain convention used throughout: `Y — X — Z`, i.e. `p(y,x,z) = p(y,x) q(z|x)`. Everything in items 1–2 and 5 is phrased on this chain; items 3–4 supply *estimators and capacity proxies* for the quantities in it; item 6 supplies the *conversion rate* between information and error.

---

## 1. Tishby–Pereira–Bialek (2000) and Tishby–Zaslavsky (2015): the IB principle

**Problem.** Given a joint law `p(x,y)`, find an encoder (transition kernel) `q(z|x)` that compresses `X` into a representation `Z` while retaining what `Z` needs to know about `Y`. Compression and relevance are both measured in the same currency, mutual information:

	min_{q(z|x)} I(X;Z)  s.t.  I(Z;Y) >= R      (IB-constrained form)
	L_beta(q) = I(X;Z) - beta * I(Z;Y)          (IB Lagrangian, beta >= 0)

with `Y — X — Z`, so by the data processing inequality (item 5) `I(Z;Y) <= I(X;Y)`: **the relevance term is capped by the relevance available at the input**, and the whole game is how much of that you keep per unit of rate. The two forms are equivalent (up to convexity issues) with `beta` the Lagrange multiplier; `beta → 0` gives `I(X;Z) → 0` (constant `Z`, pure compression), `beta → ∞` gives a (approximately) *minimal sufficient statistic* of `X` for `Y`: `I(Z;Y) = I(X;Y)` at the smallest possible `I(X;Z)`.

**The one identity to memorize.** Because `Y — X — Z`,

	I(Z;Y) = I(X;Y) - I(X;Y|Z)   =>   L_beta = I(X;Z) + beta * I(X;Y|Z) - beta * I(X;Y).

So the IB Lagrangian is, up to a constant, `rate + beta × (relevant information loss)` with
`D_rel := I(X;Y|Z)` a genuine, non-negative "distortion" on `Y`. This is why IB is a *rate–distortion* problem: under log-loss distortion, `H(Y|Z) = H(Y) - I(Z;Y) = H(Y|X) + I(X;Y|Z)`, so `D_rel` is exactly the excess log-loss distortion above the irreducible `H(Y|X)`. IB is the rate–distortion function of the *relevance variable* `Y` under logarithmic loss, and it is the natural relaxation of "minimal sufficient statistic" to **stochastic** encoders (stochasticity strictly helps: for correlated Bernoulli `(X,Y)` a randomized `q(z|x)` beats every deterministic one).

**Solution structure.** Stationary points satisfy the self-consistent equations

	q(z|x) = p(z)/Z(x,beta) * exp( -beta * D_KL[ p(y|x) || p(y|z) ] )
	p(y|z) = sum_x p(y|x) q(z|x) p(x) / p(z),     p(z) = sum_x p(x) q(z|x),

solved by a Blahut–Arimoto-type alternating iteration with **deterministic annealing in `beta`** (continue in `beta` to avoid poor local optima). Note what the exponent contains: the *KL between conditional label distributions*. IB clusters **functionals `p(y|x)`**, not points `x` — a lifted coordinate that is useless for `Y` costs rate and is pruned. That is precisely the "relevance-geometric" intuition one wants for choosing Koopman observables.

**Two dual readings, and what `beta` means.** The Lagrangian has a second, engineering reading: at an optimum the IB curve's slope is `dI(Z;Y)/dI(X;Z) = 1/beta`, so `1/beta` is the *marginal exchange rate* — nats of relevance bought per extra nat of rate. Equivalently the problem is `min I(X;Z) + beta D_rel` with `D_rel = I(X;Y|Z)`, a rate–distortion problem whose distortion is a KL between label distributions (item 6). Small `beta` sits on the steep, high-yield part of the curve; large `beta` is where the curve flattens and relevance per bit collapses — the Gaussian tail `e^{-2R}` below is exactly that regime. This is the quantitative form of the "`beta` is the compression knob" folklore and it is the reason to report whole `(I(X;Z), I(Z;Y))` pairs (a *curve*), not a single objective value.

The self-consistent equations are the stationarity conditions of `L_beta` w.r.t. `q(z|x)` under the normalization constraint `Σ_z q(z|x) = 1`: `q(z|x)` is a Gibbs/exponential tilt of the prior `p(z)` by the *negative* relevant cost `D_KL[p(y|x)||p(y|z)]`, which is the same algebraic shape as a Blahut–Arimoto channel update or an EM step — hence the alternating algorithm (update `q(z|x)` and `p(y|z)` in turn) and the standard practice of annealing `beta` from small to large values to avoid poor local optima. Semantically, `Z` is a *soft partition of the conditional label distributions*: two inputs with the same `p(y|x)` are interchangeable for the task and are merged first, which is exactly the sense in which IB clusters "functions, not points". The useful relaxed yardsticks are `ε`-sufficiency (`I(Z;Y) >= I(X;Y) - ε`) and `ε`-relevance (`I(X;Y|Z) <= ε`); all generalization arguments in this literature are stated in those terms.

**Gaussian IB (Chechik–Globerson–Tishby–Weiss 2005).** For jointly Gaussian `(X,Y)` the optimum is `T = A X + noise` (the noise is essential), `A` spanning CCA directions = the eigenvectors of `Σ_{X|Y} Σ_X^{-1}` with eigenvalues `λ_i`; the projection rank jumps discretely at `beta_i* = 1/(1-λ_i)`, so **`1/beta` acts as a soft rank selector**, a CCA-with-tunable-rank complexity knob. Scalar case (`ρ² = corr²(X,Y)`, `R = I(X;Z)`):

	I(Z;Y) = -1/2 * ln( 1 - ρ² (1 - e^{-2R}) ),   slope at R=0: 1-λ_1 = ρ²,
	I(X;Y|Z) = 1/2 ln( 1 + ρ²/(1-ρ²) * e^{-2R} )  ~  (ρ²/(2(1-ρ²))) e^{-2R}.

The IB curve is concave with slope `dI(Y;T)/dI(X;T) = 1/beta` ("diminishing returns"), and both limits are visible: linear relevance gain per nat of rate at small `R`, and **exponential decay of the missing relevance `I(X;Y|Z)` at large `R`** — the rate/error exchange rate of item 6.

**Tishby–Zaslavsky 2015 (the DL reading).** Treat the layers as a Markov chain `X — T_1 — … — T_L — Ŷ`, and plot each layer in the **information plane** `(I(X;T_l), I(T_l;Y))`. Claims: (i) `I(T_l;Y)` is non-decreasing in depth while `I(X;T_l)` is non-increasing (DPI, item 5) — "increasingly compressed sufficient statistics"; (ii) training has a fitting/drift phase (both MIs rise) followed by a compression/diffusion phase (`I(X;T)` falls, `I(T;Y)` flat) driven by SGD noise; (iii) compression is the mechanism behind generalization, with the network's effective `beta` inversely related to the injected noise level; (iv) depth buys better `relevance per rate` (the IB curve is the optimal achievable frontier). Caveats you must state if you lean on this: for *deterministic* nets with continuous inputs `I(X;T)` is degenerate (`H(X)` for discrete data, `∞` for continuous), so measured trajectories were quantization artifacts; the compression phase is *not* universal (absent with ReLU; it tracks saturating nonlinearities — Saxe et al. 2018); MI is invariant to bijective reparameterization, so it cannot by itself measure a representation's usable dimensionality (hence item 4); the "noisy DNN" reformulation adds explicit noise to make the quantities parameter-dependent.

**Why the toolbox needs it.** IB supplies (a) the objective — predict the future, pay for bits; (b) a diagnostic two-dimensional plane to report; (c) the identification "relevance term = `I(Z;Y)`" with `Y` = the future/next state in a dynamics problem; and (d) the warning that if you compute `I(X;Z)` for a *deterministic* lifting map you get `H(X)` or `∞` unless you either make the encoder stochastic (→ item 2) or use a capacity proxy such as effective rank (→ item 4).

**Which IB facts carry the weight in a paper.** Three, and they should be named explicitly: (a) the identity `I(Z;Y) = I(X;Y) - I(X;Y|Z)` — the loss-accounting statement that makes "what the lift forgot" a measurable, non-negative quantity; (b) concavity and the `1/beta` slope — which justifies annealing/sweeping `beta` and reporting the whole frontier rather than one operating point; (c) the Gaussian/CCA solution — which supplies the *closed-form baseline* that any nonlinear lift must beat at matched rate, because in the Gaussian world the optimal `beta`-tunable representation is precisely a CCA projection of tunable rank, and lower-rank projections are the dashed suboptimal segments of the IB curve. A nonlinear Koopman lift that cannot beat that curve at equal `I(X;Z)` is not buying anything beyond second-order structure.

---

## 2. Alemi et al. 2017: Deep Variational IB — making the tradeoff trainable

IB assumes known `p(x,y)` and solves a nonconvex problem with a hand-rolled iteration. Deep learning needs stochastic-gradient-compatible bounds. VIB gives two one-sided bounds pointing in *opposite* directions, which is the crux of the whole method.

**Rate side (upper bound on `I(X;Z)`).** For any fixed prior `r(z)`,

	I(X;Z) = E_{p(x)} D_KL[ q(z|x) || p(z) ]  <=  E_{p(x)} D_KL[ q(z|x) || r(z) ],
	because E_{p(x)}D_KL[q(z|x)||r(z)] = I(X;Z) + D_KL[ p(z) || r(z) ] >= I(X;Z).

Equality iff `r(z) = p(z) = ∫ q(z|x)p(x)dx`. **Minimizing the KL-to-prior term therefore shrinks an *upper* bound on the true rate** — you never over-estimate the rate, but you may under-estimate it (the term is not the rate itself). With `q(z|x) = N(μ(x), diag σ²(x))` and `r(z) = N(0,I)` (reparameterization trick for gradients),

	E D_KL = 1/2 Σ_j ( μ_j² + σ_j² - ln σ_j² - 1 ).

**Relevance side (lower bound on `I(Z;Y)`).** With a variational decoder `q(y|z)`,

	I(Z;Y) = E log q(y|z) + E log[ p(y|z)/q(y|z) ] + H(Y) >= E_{q(z|x)p(x,y)} log q(y|z) + H(Y).

Both `H(Y)` and the KL slack are constants/negligible for optimization, so **maximizing the log-likelihood of `y` from `z` maximizes a *lower* bound on relevance.**

**Objective.** With `beta_VIB` on the rate term (note the convention clash: `beta_VIB = 1/beta_IB` when comparing to `min I(X;Z) - beta I(Z;Y)`),

	L_VIB = E[ -log q(y|z) ] + beta_VIB * E[ D_KL(q(z|x) || r(z)) ],
	       >= H(Y|Z) + beta_VIB * I(X;Z) - H(Y)     (in expectation, per-sample)

i.e. a variational *upper* bound on the IB Lagrangian. Implementation: one encoder net emitting `(μ, σ)`, one decoder/classifier net, reparameterized sampling, minibatch estimate (note the well-known per-datum vs. batch-aggregate KL distinction in Alemi et al. 2018). **The encoder must be stochastic**: with a continuous deterministic encoder the KL to a continuous prior diverges, which is the practical face of the "vacuous MI" problem of item 1.

**Derivation check, and the form actually implemented.** Both bounds are one Gibbs step plus one Jensen step. *Relevance:* `H(Y|Z) <= E[-log q(y|z)]` (Gibbs), hence `I(Z;Y) = H(Y) - H(Y|Z) >= H(Y) + E log q(y|z)`, and the slack is exactly `E_z D_KL[ p(y|z) || q(y|z) ] >= 0` — the average decoder mismatch, so a stronger decoder *tightens* the bound and the slack is itself a usable diagnostic. *Rate:* the excess over the true rate is `D_KL[ p(z) || r(z) ]`, so a better prior (aggregate posterior, mixture, VampPrior) tightens the bound while leaving the objective shape unchanged. The recipe that appears in papers: encoder `q_theta(z|x) = N(mu_theta(x), diag sigma_theta^2(x))`, sample `z = mu + sigma ⊙ ε`, `ε ~ N(0,I)`, decoder `q_phi(y|z)`, and

	L = (1/N) Σ_i [ -log q_phi(y_i | z_i) + beta * D_KL( q_theta(z|x_i) || r(z) ) ],

with `beta` typically annealed `0 → beta_target` over training, and with a decision to be stated about the KL: per-datum `D_KL(q(z|x_i)||r(z))` (the standard VIB choice) versus batch-aggregate forms (the ELBO convention), which differ by a term and give different `beta`-to-rate mappings. "Free bits"/capacity-limited variants clip each dimension's KL at `lambda`, which weakens the rate bound slightly but protects against posterior collapse.

**Evidence that the term does work.** VIB's original experiments are why this machinery is borrowed so widely: adding the KL-rate term produced representations that generalized better *and* were markedly more robust to adversarial perturbation than the same architecture without it, with the effect monotone in `beta` over a range (too large `beta` → the encoder stops transmitting and the model underfits). The practical lesson for an IB-Koopman paper is that compression is a data-dependent regularizer whose strength is chosen by `beta`, not by architectural width.

**How this differs from a plain autoencoder, and the extremes.** A β-VAE-style objective reconstructs `X` from `Z` and is the case `Y = X`, so its rate term does exactly the same job (an upper bound on `I(X;Z)`) but its relevance term is *reconstruction* rather than task relevance — which is why a generative objective can spend bits on nuisance that a discriminative/IB objective would prune. VIB's discriminative `Y` is what makes the two terms fight over *task-relevant* information. The two extremes are instructive: `beta_VIB → 0` gives a high-rate, essentially lossless channel (rate bound unrestricted, nothing compressed); `beta_VIB → ∞` drives the KL to zero, `q(z|x) → r(z)`, no bits flow, and the relevance bound collapses to `H(Y)` (underfitting via posterior collapse). The useful operating régime is the elbow of the `(R, D)` tradeoff, where the *marginal* drop in distortion per nat of rate stops paying — the same statement as `1/beta` being the IB curve's slope (item 1).

**Diagnostics and failure modes.** Alemi et al. 2018 ("Fixing a broken ELBO") recast training as rate–distortion: rate `R ≈ E D_KL(q(z|x)||r(z))` (the upper bound above), distortion `D = E[-log q(x|z)]` (or the label analogue), and plot trajectories in the `(R, D)` plane, revealing three regimes: (i) `beta` too large → rate collapses, the channel is unused (posterior collapse); (ii) `beta` too small → high rate, uninformative bits about `Y` accumulated; (iii) the "good" regime tracing the `R(D)` frontier. Useful facts for a paper: the reachable `(R,D)` frontier is the rate–distortion curve of item 6; at a *fixed* rate, a model may still fail to *use* the bits (`I(X;Z)` high, `I(Z;Y)` low) — "the elbow" is where you want to sit.

**Why the toolbox needs it.** VIB is the only reason the `beta` in an IB-Koopman objective is a gradient-trainable knob at all; it tells you which direction each logged quantity errs (`KL` is an upper bound on rate, `log q(y|z)` a lower bound on relevance, so both are conservative in the *right* direction for the tradeoff); and it tells you exactly what to make stochastic (the lifting to `z`) so that "information" is a well-defined, parameter-dependent quantity rather than `H(X)`.

**Numbers to report.** (i) the rate *upper bound* `E D_KL` — labelled as an upper bound, not as the rate; (ii) the relevance *lower bound* `E log q(y|z)` together with the decoder/classifier that produced it; (iii) the `(R, D)` pairs for the whole `beta` sweep, including the trivial endpoints; (iv) a proxy for `I(X;Y)` so that the deficit `I(X;Y) - I(Z;Y) = I(X;Y|Z)` is explicit and interpretable in nats. That set makes the paper's claims falsifiable without ever computing an exact mutual information.

---

## 3. InfoNCE (van den Oord et al. 2018; Wu et al. 2020) and MI estimation by contrast

**MI as an optimal log-ratio (Donsker–Varadhan).** `I(X;Y) = sup_f E_{p(x,y)}[f(x,y)] - ln E_{p(x)p(y)}[e^{f(x,y)}]`, attained at `f* = ln[p(x|c)/p(x)]` (up to a function of `c`). Contrastive learning is the multi-sample, empirically estimable version of this statement.

**The estimator.** Draw a positive pair `x_1 ~ p(x|c)`, `c ~ p(c)`, and `K-1` negatives `x_{2:K} ~ p(x)` i.i.d. With critic `f` define the softmax classifier over candidate positions and its cross-entropy loss

	L_NCE = E[ -f(x_1,c) + ln Σ_{j=1}^{K} e^{f(x_j,c)} ],     Î_NCE := ln K - L_NCE.

**Why `ln K - L` is a lower bound.** (i) Gibbs' inequality: the cross-entropy of *any* predictor of the positive index `I ∈ {1..K}` is at least the conditional entropy, `L_NCE >= H(I | data)`. (ii) `H(I) = ln K` and `I(I;data) <= I(X;C)` (the index is a function of the sampled pair and negatives). Hence

	L_NCE >= ln K - I(X;C)   <=>   I(X;C) >= ln K - L_NCE = Î_NCE,   with Î_NCE <= ln K.

The optimal critic `f(x,c) = ln[p(x|c)/p(x)] + const(c)` makes the bound tight: `L_NCE = ln K - I(X;C)` exactly. Two consequences that matter: **(a) the estimator saturates at `ln K`** — with a batch of `K` you cannot measure (and, via gradient saturation, cannot chase) MI much above `ln K`; **(b) it is a *bound*, not an estimate**: it is only as good as the critic's ability to represent the true density ratio, so reported values are floors. In the ML-visible convention (Wu et al. 2020, Eq. 5), the same object is written as the unnormalized form with an empirical partition estimate:

	Î_NCE = E_{p(x,y_1)}[ f(x,y_1) - E_{p(y_{2:K})} ln ( (1/K) Σ_j e^{f(x,y_j)} ) ]  <= I(X;Y_1),

with dual encoders `f(x,y) = g_θ(x)^T g_φ(y)`; the constant shifts by `ln K` between conventions, and the K-including-positive-vs-leave-one-out partition choice is exactly the source of the `ln K` bias noted in this literature.

**CPC: temporal-neighborhood positives.** Encoder `z_t = g_enc(x_t)`, autoregressive context `c_t = g_ar(z_{<=t})`, critic `f_k(x_{t+k}, c_t) = z_{t+k}^T W_k c_t` (a *different* `W_k` per offset), loss = InfoNCE over `(z_{t+k}, c_t)` positives and negatives drawn from other time steps/sequences in the batch. Maximizing `Î_NCE` maximizes a lower bound on `I(z_{t+k}; c_t)`. Reasons for this design: (i) the quantity is *predictive* mutual information — the exact analogue of the IB relevance term with `Y` = future observation, so it is the natural unsupervised stand-in for the task MI; (ii) the context forces the positive to be predictable *from the past* rather than merely equal to its own augmentation; (iii) `k > 1` suppresses instantaneous nuisance (sensor noise, high-frequency content) that is present in `I(z_t; z_{t+1})` but not in `I(z_{t+k}; c_t)` — CPC's empirical finding is that the usefulness grows with `k` up to a point; (iv) `k` gives you a *monotone profile*: by DPI along `z_t — x_t — x_{t+k}`, `I(z_{t+k}; c_t)` is non-increasing in `k` for a fixed encoder geometry, so the curve `k ↦ Î_NCE(k)` is a readable "predictability spectrum" — directly reusable as a diagnostic for how many nats of the future a Koopman lift has captured. In a dynamics setting the temporal-positive construction is also what makes the *linear* structure meaningful: contrastive training pushes the lift toward coordinates whose future is linearly decodable from the past by the critic, i.e. toward an approximately invariant/epoch-free observable space.

**Conventions, numerics, and what actually moves the bound.** (i) *Which partition:* with the positive included in the softmax denominator the loss is `L = E[-f(x_1,c) + ln Σ_{j=1}^K e^{f(x_j,c)}]` and `I >= ln K - L`; with a leave-one-out denominator (only negatives) the constant becomes `ln(K-1)` and the estimator is marginally less biased. State the convention. (ii) *Temperature:* `f = z^T W c / τ` rescales logits; `τ` does not change the math of the bound, it changes critic conditioning (and thus how close the bound sits to `I`). (iii) *Normalization:* `ℓ2`-normalizing `z` and `c` bounds the logits, stabilizes `logsumexp` and prevents a few large scores from dominating — but it also caps achievable `f` and hence the estimable MI, so "at least `ln K - L` nats, up to `ln K`" is the only defensible phrasing. (iv) *Effective K:* memory queues / momentum encoders raise the effective number of negatives without enlarging the batch, raising the `ln K` ceiling. (v) *Variance dominates bias in practice:* the dominant error is the heavy-tailed variance of `ln (1/K) Σ e^f`, cured by more negatives rather than a bigger critic; Wu et al. also show that replacing the original Instance-Discrimination softmax with hardcoded normalization constants by plain InfoNCE *improves* transfer accuracy, i.e. the MI reading is also the numerically better implementation.

**Where DPI caps the temporal MI.** In CPC's chain `z_t — x_t — x_{t+k}`, item 5 gives `I(z_{t+k}; c_t) <= I(x_{t+k}; c_t)`, and for a fixed encoder geometry `k ↦ I(z_{t+k}; c_t)` is non-increasing. Contrastive training therefore maximizes a lower bound on a quantity that is itself capped by the state's own predictability, and the measured `Î_NCE(k)` is a lower bound of a lower bound. The defensible claim is "at least this many nats of the future are *linearly* (or neurally) accessible from the lift", never "the lift is sufficient".

**Wu et al. 2020 (what the estimator is really doing).** They (i) state the bound above with a bilinear witness; (ii) show that Instance Discrimination, CMC, SimCLR and Local Aggregation are all *the same* bound with different "views" (augmentation sets) and negative distributions — IR ≡ InfoNCE + `ln K + const`; (iii) prove that negatives may be drawn from a *restricted* conditional `q_T(y_{2:K}|y_1)` and still yield a valid lower bound (VINCE), which justifies "hard negative" schemes (ball/ring/cave sampling) and yields a general family; (iv) prove a permutation-invariance result: without augmentations (views), the objective only constrains the *arrangement* of embeddings, not their semantic grouping — the information that gets retained is decided by the **view/positive definition**, not by the objective. This is the sharpest available argument that "the choice of positives is the choice of what information the representation keeps".

**Caveats to quote.** (1) Estimation bias/variance: variational MI estimators are loose in high dimension and the bounds blow up in variance with weak critics (Poole et al. 2019, which unifies BA/DV/NWJ/InfoNCE/JSD and shows e.g. the JSD bound caps at `ln 2`); any estimator from finite samples has fundamental limits (McAllester & Stratos 2020); the practical looseness of InfoNCE with neural critics is quantified in Song & Ermon 2020. (2) MI alone does not explain transfer performance (inductive bias/augmentation-aware arguments exist, e.g. Wang et al. 2022) — so treat InfoNCE values as *relative* diagnostics under a fixed architecture/positive scheme, never as absolute information budgets.

**Why the toolbox needs it.** InfoNCE is the only tractable, single-sample-estimable, gradient-friendly handle on `I(Z;Y)`. In an IB-Koopman paper it plays two roles: it is the *training signal for the relevance term* when labels/futures are high-dimensional, and it is the *measurement instrument* for "how much of the future the lifting carries".

**Reporting protocol.** Freeze the positive scheme (temporal offsets `k`, augmentations/views), the negative pool, the temperature, the normalization, and the batch size `K`; then report `Î_NCE(k)` curves for `k = 1, 2, 4, 8, …` as *relative* diagnostics under that single protocol, alongside (a) a linear-probe or one-step prediction error, which separates "accessible" from "available", and (b) a decoder-based relevance bound (item 2), so that at least one number is a genuine bound on task MI rather than on the contrastive surrogate. Finally, note that MI maximization alone does not explain transfer performance — the augmentation-reliance results of Wu et al. and the theoretical critique of the "MI-only" account show that inductive bias carries part of the load — so InfoNCE should be presented as one aligned signal among several, not as the explanation of success.

---

## 4. von Neumann entropy, effective rank / effective dimension, Roy–Vetterli (2007)

**Von Neumann entropy.** For a density matrix `ρ ⪰ 0`, `Tr ρ = 1`, with spectral decomposition `ρ = Σ_i λ_i |v_i><v_i|` (`λ_i >= 0`, `Σ λ_i = 1` — the spectrum *is* a probability vector),

	S(ρ) = -Tr(ρ ln ρ) = -Σ_i λ_i ln λ_i = H(λ).

Facts: `S(ρ) ∈ [0, ln d]` with `S = 0` iff `ρ` is pure (rank 1) and `S = ln d` iff `ρ = I/d` (maximally mixed); unitary invariance (basis-independent); concavity in `ρ`; additivity `S(ρ_A ⊗ ρ_B) = S(ρ_A) + S(ρ_B)`; `S(ρ) <= ln rank(ρ)`; and, writing `S(ρ) = ln d - D(ρ || I/d)`, it measures "distance from maximally mixed", i.e. **how many directions the state actually occupies**. In nats it is the entropy of the eigenvalue distribution — which is exactly the object people compute for a normalized feature covariance/Gram matrix in representation learning, with `λ_i/Tr Σ` as the eigenvalues.

**Effective dimension.** `d_eff(ρ) := exp(S(ρ)) ∈ [1, d]`: the number of "equally weighted directions" with the same entropy, a smooth, differentiable surrogate for rank. Companion proxies: participation ratio `1/Tr(ρ²) = 1/Σλ_i²`, which satisfies `ln(1/Tr ρ²) <= S(ρ) <= ln rank ρ` (Rényi monoticity/Schur-concavity), and stable rank `‖A‖_F²/‖A‖_2²` for non-normalized matrices.

**Roy–Vetterli 2007 (general matrix version).** For `A` with singular values `σ_1 >= … >= σ_r > 0`, set `p_i = σ_i / ‖σ‖_1` and define

	erank(A) = exp( -Σ_i p_i ln p_i ) = exp( H(p) ).

Properties (this is why it is the right extension of "rank"): `1 <= erank(A) <= rank(A) <= min(m,n)`; `erank` is continuous and **scale-invariant** (`erank(cA) = erank(A)`); `erank(A) = rank(A)` iff all nonzero singular values are equal, and `erank → 1` as the spectrum degenerates onto one direction; it interpolates "rank 1 for a single dominant direction … rank `r` for a flat spectrum". A clean quantitative reading (to use when you need a bound rather than intuition): by Schur-concavity/Rényi monotonicity,

	erank(A) >= 1 / Σ_i p_i² = ‖A‖_*² / ‖A‖_F²       (nuclear-to-Frobenius ratio, equality for flat spectra),

i.e. effective rank is at least the "flatness" of the spectrum, with the top-`k` energy fraction bounded through the usual Eckart–Young decomposition — so `erank` is precisely a *concentration* functional of the spectrum. Roy & Vetterli's motivation is numerical/algebraic: integer rank is discontinuous under perturbation, so rank-revealing algorithms need a real-valued notion of "how many singular values are meaningful"; `erank` gives that with an operational reading in terms of the number of directions required to reach a prescribed relative accuracy.

**Operational reading and relatives.** `exp(H(p))` is the size of the *typical set* of the eigenvalue distribution — the number of directions that carry non-negligible mass, as opposed to the number that are merely nonzero — so `erank` counts "directions that matter". It lives in the Rényi family `H_alpha(p) = (1/(1-alpha)) ln Σ_i p_i^alpha`; monotonicity `H_0 >= H_1 >= H_2 >= …` gives the sandwich

	ln(1/Σ_i p_i²) = H_2 <= S = H_1 <= ln rank(A),

which yields the computable lower bound `erank(A) >= 1/Σ_i p_i²` quoted above (equality for flat spectra). Related surrogates in use: stable rank `‖A‖_F²/‖A‖_2²`, numerical rank at a threshold, the "effective dimension" of a Fisher/feature matrix in generalization bounds, and the participation ratio `1/Tr(ρ²)`. Two implementation notes: eigen-decompose the symmetric part of `A^T A` (guard against spurious negative eigenvalues), and always state the log base — `exp(H)` in nats versus `2^{H}` in bits differ by a factor `1/ln 2 ≈ 1.44` in the exponent, a classic source of confusion.

**Why this is the right capacity proxy for MI.** For a Gaussian latent `Z`, `h(Z) = 1/2 ln((2πe)^d det Σ_z) = 1/2 Σ_i ln(2πe λ_i)`; by AM–GM, at fixed trace (fixed total power) this is *maximized by a flat spectrum*, i.e. by the largest `erank`. So `erank` is the scale-free coordinate that governs how much volume — hence how many distinguishable codes, hence how much rate — a given amount of variance can support. Compression wants the spectrum squeezed (`erank` small), expressivity wants it fat. That is the same tradeoff as `I(X;Z)`, expressed in a quantity that remains finite, differentiable, and meaningful for a *deterministic* lifting map — precisely the case in which items 1 and 5 leave MI degenerate.

**Why ML wants it, and why an IB paper needs it.** (i) *Capacity proxy for a rate whose exact value is degenerate*: for a deterministic lifting map `φ`, `I(X;Z)` is `H(X)` or `∞` (item 1), but the spectrum of the lifted feature covariance still tells you the *effective number of independent lifted coordinates* — the natural, differentiable stand-in for `R`. (ii) *Rate vs. capacity bookkeeping*: for Gaussian features, `h(Z) = 1/2 ln((2πe)^d det Σ_z)` is maximized, at fixed total power `Tr Σ_z`, by a flat spectrum, i.e. by maximal effective rank (AM–GM/Hadamard), so `erank` controls the achievable log-det (log-volume) at fixed "energy" — compression wants a *low* `erank` (concentrated spectrum), expressivity wants a *high* one, and the IB tradeoff becomes a tradeoff in a scale-invariant, computable quantity. (iii) *Diagnostics for rank collapse and for identity of the latent*: `erank → 1` means the model has collapsed to a single direction (constant feature); `erank` of the Graham/eigen-spectrum is a standard tool for measuring collapse in depth and for counting *active modes* in a latent linear dynamical representation (how many Koopman modes actually carry variance). (iv) *Caveat to state explicitly*: `erank` is scale-invariant, so it is a measure of spectral *shape*, not of magnitude/rate; report it alongside `Tr Σ_z` (or `ln det`) and never as a substitute for a bound on `I(X;Z)`.

**Relation to the other sections.** `S(ρ)` and `erank` are the "information-theoretic but geometry-aware" measure that fills the gap left by DPI's invariance argument (item 5): MI cannot see the difference between an uninformative `z` and a rich one under bijective reparameterization, while `erank`/`S` can. For quantum channels, the analogue of DPI is monotonicity of relative entropy: since `D(Λρ ‖ I/d) <= D(ρ ‖ I/d)` for any trace-preserving map, any **unital** channel (`Λ(I/d) = I/d`) satisfies `S(Λρ) >= S(ρ)`, i.e. mixing never reduces the effective dimension, whereas non-unital compression/decoherence can. That "second law"-like asymmetry is a useful sanity check when interpreting `erank` trajectories during training.

**Concrete usage in a Koopman paper.** (i) Report the eigenvalue spectrum of the lifted feature covariance (or the singular values of the empirical Koopman operator) with `erank` as its scalar summary: this is the "how many effective observable coordinates / active modes are in play" number and it is the honest complement to any claimed compression rate. (ii) Track it jointly with the relevance loss: `erank → 1` while relevance saturates means over-compression (the reverse of item 2's failure mode (ii)); `erank` pinned at the ambient dimension while future-MI is flat means rate is being spent on nuisance directions. (iii) Never present `erank` as a rate; pair it with `ln det Σ_z` or `Tr Σ_z`. (iv) When the point of the paper is that the lift is *low-dimensional*, the precise claim reads "the lifted representation occupies `erank` effective directions while retaining `Î_NCE` nats about the future", which is a joint statement about shape and content — the two things MI alone cannot state together.

---

## 5. The data processing inequality and information-loss bounds in encoder–decoder pipelines

**Statement and proof.** If `X — Y — Z` (i.e. `X ⫫ Z | Y`), then `I(X;Z) <= I(X;Y)` and `I(X;Z) <= I(Y;Z)`. Proof of the first: expand two ways,
`I(X;Y,Z) = I(X;Z) + I(X;Y|Z) = I(X;Y) + I(X;Z|Y) = I(X;Y)`, using `X ⫫ Z | Y` for the last term. Hence

	I(X;Z) = I(X;Y) - I(X;Y|Z)   =>   I(X;Y|Z) = "information lost about X by going through Z" >= 0.

**Equality conditions (the useful part).** `I(X;Z) = I(X;Y)` iff `I(X;Y|Z) = 0` iff `Z` is a *sufficient statistic* of `X` for `Y` (iff `X — Z — Y` is also Markov). So "no information about `Y` is lost" is a checkable Markov condition, and the loss is quantified by the conditional MI. Replacing `X` by the target `Y` in the chain `Y — X — Z` gives the IB-relevant statement: `I(Z;Y) = I(X;Y) - I(X;Y|Z)`, i.e. **the relevance deficit is exactly the conditional MI**, and it is non-negative — this is what makes "information-loss bound" a well-posed phrase in an encoder–decoder pipeline.

**Scope and strengthenings.** (i) DPI is not special to mutual information: it holds for all `f`-divergences (and for Rényi divergences of order `alpha > 1/2`), which is why loss bounds of the same shape can be stated for `χ²`, Hellinger or total variation when those are the estimable quantities. (ii) The *strong* DPI holds for `f`-divergences, so the monotonicity applies to every internal pair of the pipeline, not only to the endpoint. (iii) Composition is free: `X — Z — W` gives `I(X;W) <= I(X;Z) <= I(X;Y)`, so a chain of encoder → projection → quantizer → noise can be analyzed pairwise, and *adding independent noise* at any stage is just one more Markov step that can only reduce relevance. (iv) The quantum counterpart (Uhlmann/Lindblad monotonicity of relative entropy under completely positive, trace-preserving maps, plus the operator Schwarz inequality for such maps) is the precise form of "no physical operation creates correlations"; for unital channels it implies `S(Lambda ρ) >= S(ρ)`, i.e. the von Neumann-entropy analogue of a second law used in item 4. (v) The chain rule `I(X;Y,Z) = I(X;Z) + I(X;Y|Z) = I(X;Y) + I(X;Z|Y)` is the algebraic workhorse that converts an IB objective into "rate plus relevant distortion".

**Consequences for pipelines (all of these get used in practice).**
1. *Non-increasing along depth/stacking:* for `Y — X — T_1 — … — T_L`, both `I(X;T_l)` and `I(T_l;Y)` decrease with `l`. You cannot add information about `X` or `Y` by adding layers or by post-processing; downstream performance ceilings are set upstream. In particular `I(future ; decoder_output) <= I(future ; z) <= I(future ; x)` — the decoder cannot manufacture what the encoder dropped.
2. *Invariance to bijections:* `I(T;Y) = I(ψ(T);Y)` for invertible `ψ` — the layer-wise "compression" measured by MI in deterministic nets is vacuous (item 1), and MI cannot tell a whitened/scaled/batched representation from the original. Consequence: MI can bound what is *available*, never what is *linearly accessible*; if your claim is about linear accessibility (Koopman!), you need an extra argument (linear probe/CCA/eigenvalue error), not MI alone.
3. *Quantitative error bounds.* Fano: `H(Y|Z) <= h_b(P_e) + P_e ln(|Y|-1)`, so `P_e >= (I(X;Y) - I(Z;Y) - 1)/ln|Y|`-type floors — any decoder's error is lower-bounded by the *retained* relevance. In continuous/MSE form, the entropy-power (Bretagnolle–Huber) inequality gives
   `E‖X - g(Z)‖² >= (1/2πe) exp(2 h(X|Z)) >= (1/2πe) exp(2 h(X) - 2 I(X;Z))`,
   an explicit exponential "information-loss bound" for encoder–decoder pipelines, and the natural bridge to item 6.
4. *Generalization (why IB people care):* `I(Z;Y) >= I(X;Y) - ε` ⟺ `I(X;Y|Z) <= ε`, and combining with Fano gives risk bounds `risk <= empirical risk + sqrt((I(X;Z) + ln(1/δ))/(2n))`-type statements (Shamir–Sabato–Tishby 2008; the "information / PAC-Bayes" line, Xu & Raginsky 2017 for the algorithm-output version `gen gap <= sqrt(I(W;S)/2n)`). So compression (`I(X;Z)`) is a *capacity* term and the relevance deficit an *approximation* term — exactly the two error sources of a learning bound.
5. *Non-increase of MI is the reason "views" matter:* if a view/augmentation drops label-relevant information, then `I(view; Y) <= I(x; Y)` and no contrastive critic can recover it (this is DPI applied to the view map); the contrastive objective is bounded by the view's own information content — the formal counterpart of Wu et al.'s "views are the crux".

**Nuisance/informativeness tradeoff (the design corollary).** Apply DPI to a task variable `Y`, a nuisance variable `N`, and the representation `Z`: any reduction of `I(Z;N)` achieved by *discarding* information about `X` also caps `I(Z;Y)` unless the discarded subspace is specifically the `N`-only subspace. This is the formal content of the standard invariance/informativeness tension (Achille & Soatto's information-dropout and "emergence of invariance and disentanglement", and the "sufficiency" view of representation learning): invariance to nuisance is paid for in relevance, and the cheapest payment is to align the discarded directions with the nuisance directions — which is exactly what a contrastive objective with well-chosen positives/views, or an IB objective with a well-chosen `Y`, induces. In a dynamics pipeline, if the lift is required to be invariant to an operating-condition nuisance, the predicted cost is exactly `I(X;Y|Z)` and it should be *reported*, not assumed to be free.

6. *Discrete capacity ceilings (immediate corollary).* If `Z` is a finite code (`|Z| = K`), then `I(Z;Y) <= H(Z) <= log K`: a `K`-way quantized latent can carry at most `log K` bits about anything downstream. Similarly a rank-`r` linear lift in `R^d` cannot retain more than the information in its `r`-dimensional image. Both statements are DPI plus the trivial entropy bound, they are checkable in one line, and they are the most common "your architecture is information-starved" findings in practice — which is why the effective rank of item 4 is the natural continuous companion to the integer `K` here.

**Why the toolbox needs it.** DPI is the *impossibility theorem* that keeps an IB-Koopman paper honest: it converts "we compress the lift" into "here is the exact quantity of relevance we forfeit", provides the exact identity used to write the loss, converts to a hard floor on prediction error (items 5.3/6), and explains why the paper must separately address (i) stochasticity (so MI is finite), (ii) linear accessibility (which DPI-invariant MI cannot see), and (iii) the fact that no decoder or contrastive head can exceed the MI retained by the encoder.

---

## 6. Rate–distortion basics, and the `D ~ exp(-I)` law

**Definition.** For a source `X`, reproduction `X̂`, distortion `d(x,x̂)`:
`R(D) = inf_{q(x̂|x) : E[d] <= D} I(X;X̂)`, `D(R) = R^{-1}(D)` (Shannon's source-coding theorem: `R(D)` is the infimum of achievable rates; the converse says no code does better — this is the operational content of item 5's DPI bound).

**Gaussian template (memorize this one).** `X ~ N(0, σ²)`, squared error, base-2 logs:
`R(D) = 1/2 log(σ²/D)`, `D <= σ²`; `D(R) = σ² 2^{-2R} = σ² e^{-2R}` (nats). In general (Shannon lower bound + Zador high-rate / Shannon–Bucklew),
`D(R) >= (1/2πe) e^{2h(X) - 2R}` (nats, MSE), with equality for a Gaussian source, and for a `d`-dimensional source `D(R) ≈ c_d · ‖p‖_{d/(d+2)} · e^{-2R/d}` — **exponential decay with the exponent inversely proportional to the source dimension** (shoulder/plateau effects appear when the source concentrates near a lower-dimensional set, which is exactly the "effective dimension" of item 4 entering the exponent). Two readings, both useful:
- *Cost per unit of accuracy:* to reduce error by a factor `e²`, you must buy one more nat of rate; to reduce it by `10×` you need `R ≈ 0.5 ln 10 = 1.15` nats more. **Rate is logarithmic in precision** — improvement is cheap in bits and expensive in "distortion decades".
- *Budget statement:* if the pipeline can pass at most `R` nats (by construction, by the bottleneck, or by `I(X;Z)` measured as in items 2–3), then the best possible decoder MSE is `≈ σ² e^{-2R}` — a *hard* floor independent of architecture or training.

**Where `D(R) = σ² e^{-2R}` comes from, and its finite-rate refinements.** The Gaussian case is solved by a Lagrangian over `q(x̂|x)` using `I(X;X̂) = h(X̂) - h(X̂|X)` with a second-moment constraint; the result is the memorized formula. For a vector source with independent components of variances `λ_i`, **reverse water-filling** gives `R(D) = Σ_i 1/2 ln(λ_i / D_i)` with `D_i = min(λ_i, θ)` where `θ` is chosen so that `Σ D_i = D`: components below the water level are *never coded at all*. That is the multivariate face of "the dimension enters the exponent", and it is the mechanism behind shoulders/plateaus in measured `R(D)` curves for sources that concentrate near a lower-dimensional set. Achievability is by random codebooks/typicality; the converse is the Shannon lower bound `R(D) >= h(X) - 1/2 ln(2πe D)`, from which `D >= (1/2πe) e^{2h(X) - 2R}` follows directly (entropy-power inequality). At high rate, Zador's formula `D ≈ c_d · ‖p‖_{d/(d+2)} · e^{-2R/d}` fixes the exponent at `2/d` for a `d`-dimensional MSE source, and the familiar `(πe/6) ≈ 1.53 dB` gap of scalar uniform quantization says that practical quantizers sit a *constant* factor above the law but cannot change the exponent — a useful reality check when comparing compressor families. Convexity of `R(D)` and concavity of `D(R)` (true for MSE and Gaussian sources) are what make the Lagrangian/`beta`-sweep machinery of items 1–2 coherent.

**The `D ~ exp(-I)` bound, stated three ways (use the one that matches the paper's distortion).**
1. **MSE / continuous, any pipeline:** since `I(X;Z) <= R` for any `z`-based pipeline by DPI,
   `E‖X - g(Z)‖² >= (1/2πe) e^{2h(X) - 2 I(X;Z)}`, i.e. `D = Θ(e^{-2I})`. Error decays exponentially in the *mutual information the representation carries*; equivalently, information about the target must grow logarithmically in `1/D`. For a jointly Gaussian `(state, lift)` pair this is exact: `MSE ≃ σ² e^{-2I}`.
2. **Gaussian IB / relevant-information loss (the sharpest version for an IB paper).** With the scalar Gaussian IB formula of item 1,
   `I(X;Y|Z) = 1/2 ln(1 + (ρ²/(1-ρ²)) e^{-2R}) ~ (ρ²/(2(1-ρ²))) e^{-2 I(X;Z)}` as `R → ∞`.
   So **the information lost about the relevance variable decays exponentially in the rate** — the exact "`D ~ exp(-I)`" statement inside the IB framework itself, with the prefactor `ρ²/(2(1-ρ²))` set by how much relevance was available (`ρ²`) in the first place. Corollary: to keep, say, `ε` nats of relevance, the rate must be `R ≈ 1/2 ln(ρ²/(2(1-ρ²)ε)) = O(ln(1/ε))` — logarithmic in the target error.
3. **Log-loss distortion (the IB-native case), for contrast.** Under logarithmic loss, `D = H(Y|Z) = H(Y) - I(Z;Y) = H(Y|X) + I(X;Y|Z)` — the distortion is *linear* (affine) in relevance and floored at `H(Y|X)`. So when a paper reports `H(Y|T)`-type distortion, expect a straight line, not an exponential; when it reports MSE/generalization error against measured MI, expect `e^{-2I}`. Mixing the two is a common and detectable error.

**Rate–distortion as the diagnostic plane.** Alemi et al.'s `(R, D)` plane (item 2), the IB curve `(I(X;Z), I(Z;Y))` (item 1), and the Zador/Gaussian predictions above are three views of the same frontier: only points on/above `D(R)` are achievable; a *learned* representation sitting well above the frontier is preserving bits that do not help (item 2's failure mode) and one sitting on it is at the optimal accuracy/complexity exchange. For an IB-Koopman paper, `D(R) = σ²e^{-2R}` also supplies the practical justification for the chosen `β`: the marginal value of the last nat of rate is `|dD/dR| = 2σ²e^{-2R}` — small once `R` exceeds `≈ 1/2 ln(2σ²/ε_target)`, so beyond that point the IB term should win.

**Fitting `beta` from this law.** Three concrete consequences for the paper. (1) The marginal distortion reduction per nat is `2D`, so taking distortion from `D_0` to `ε` costs `R = 1/2 ln(D_0/ε)` — one nat per factor `e² ≈ 7.4`; `beta` should therefore be annealed on a logarithmic scale in the residual error, not a linear one. (2) Because the exponent involves the source's effective dimension (item 4), spending rate beyond that dimension buys nothing — expect plateaus. (3) Any claimed improvement in prediction error can be restated as an information gain: `ΔR = 1/2 ln(ε_baseline / ε_method)`, so the claim becomes "the IB-Koopman lift retains `ΔR` nats more about the future than the linear/CCA baseline, which is why its error is `e^{-2ΔR}` smaller" — architecture-independent and falsifiable.

**Why the toolbox needs it.** It is the only bridge from "nats kept" to "error achieved": it converts the IB/InfoNCE quantities of items 1–3 into a reported, auditable claim about prediction quality, tells you the exchange rate (logarithmic rate for exponential accuracy), warns you about the dimension in the exponent (item 4), and gives a principled stopping rule for the compression/relevance tradeoff.

---

## Synthesis: which tool does what in an IB-Koopman pipeline

| Need | Tool |
|---|---|
| Objective: predict the future, pay for bits | IB Lagrangian `I(X;Z) - beta I(Z;Y)` (1) |
| Trainable `beta`, two-sided error bookkeeping | VIB: KL-to-prior upper-bounds `I(X;Z)`, `log q(y|z)` lower-bounds `I(Z;Y)` (2) |
| Unsupervised relevance signal for high-dim futures + measurement instrument | InfoNCE with temporal-neighborhood positives; `Î = ln K - L <= I`, ceiling `ln K` (3) |
| Finite, geometry-aware capacity/rate proxy when MI degenerates; "how many latent modes" | `S(ρ) = -Tr ρ ln ρ`, `erank(A) = exp H(σ/‖σ‖_1)` (4) |
| Why the loss has the form it has; hard floors on any decoder | DPI: `I(Z;Y) = I(X;Y) - I(X;Y|Z)`, Fano, `(1/2πe)e^{2h(X)-2I(X;Z)}` (5) |
| Nats → error conversion; `beta` tuning; sanity-check of reported curves | `D(R) = σ²e^{-2R}`; Gaussian IB `I(X;Y|Z) ~ (ρ²/2(1-ρ²))e^{-2I(X;Z)}` (6) |

Three traps worth writing into the paper: (a) a *deterministic* encoder makes `I(X;Z)` vacuous — go stochastic or use `erank`; (b) InfoNCE values are *lower bounds capped at `ln K`* — compare only within a fixed batch/positive scheme; (c) DPI invariance means MI cannot certify *linear* accessibility, which is precisely what a Koopman claim needs.

### Key sources
- Tishby, Pereira, Bialek, *The Information Bottleneck Method*, arXiv:physics/0004057 (2000); Tishby & Zaslavsky, *Deep Learning and the Information Bottleneck Principle*, ITW 2015, arXiv:1503.02406; Goldfeld & Polyanskiy, *The Information Bottleneck Problem and Its Applications in ML*, arXiv:2004.14941 (tutorial, incl. Gaussian IB and the caveats).
- Chechik, Globerson, Tishby, Weiss, *Information Bottleneck for Gaussian Variables*, JMLR 6:165–188, 2005.
- Alemi, Fischer, Dillon, Murphy, *Deep Variational Information Bottleneck*, ICLR 2017, arXiv:1612.00410; Alemi et al., *Fixing a Broken ELBO*, ICML 2018.
- van den Oord, Li, Vinyals, *Representation Learning with Contrastive Predictive Coding*, arXiv:1807.03748; Wu, Zhuang, Mossé, Yamins, Goodman, *On Mutual Information in Contrastive Learning for Visual Representations*, arXiv:2005.13149; Poole et al., *On Variational Bounds of Mutual Information*, ICML 2019; McAllester & Stratos, *Formal Limitations on the Measurement of Mutual Information*, AISTATS 2020; Song & Ermon, *Understanding the Limitations of Variational MI Estimators*, ICLR 2020.
- Roy & Vetterli, *The Effective Rank: A Measure of Effective Dimensionality*, EUSIPCO 2007.
- Cover & Thomas, *Elements of Information Theory* (DPI, Fano, rate–distortion); Xu & Raginsky, *Information-Theoretic Analysis of Generalization*, NeurIPS 2017; Shamir, Sabato, Tishby, *Learning and Generalization with the Information Bottleneck*, 2008; Saxe et al., *On the Information Bottleneck Theory of Deep Learning*, ICLR 2018.
