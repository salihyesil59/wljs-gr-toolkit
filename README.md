# General Relativity notebooks

Wolfram Language notebooks for symbolic general relativity and modified gravity, written for
[WLJS Notebook](https://jerryi.github.io/wljs-docs/). Give them a metric or an action and they
hand back curvature tensors, field equations, Friedmann equations and observables — non-zero,
independent components only — and then keep going: cosmological distances, structure growth,
the effective gravitational coupling, the teleparallel and symmetric-teleparallel families, a
Hamiltonian mode count, the price of the approximation all of it rests on, what a gravitational
wave does on the way here, why the laboratory has not already ruled the whole thing out, and
which of the three families a standard siren could actually catch.

The through-line is that every result is checked against a limit where the answer is already
known. Each notebook ends with verification cells, and those cells are not decoration: they
caught a silent iterator collision that had been returning $R = 0$ for FLRW, an
Euler–Lagrange operator that dropped the term carrying $k^2$, an $f(\mathcal{T})$ growth
solution that had walked through a pole, and a degree-of-freedom count that gave general
relativity a propagating scalar it does not have. A method that loses general relativity
cannot be trusted with anything else.

| | in | out |
|---|---|---|
| **GR-01** | a metric | Christoffel, Riemann, Ricci, Einstein, Kretschmann |
| **GR-02** | an action | field equations, Friedmann equations, $E(z)$ |
| **GR-03** | $E(z)$ | distances, times, horizons, the BAO ruler |
| **GR-04** | $E(z)$ and $G_{\rm eff}$ | growth, $f\sigma_8$ |
| **GR-05** | perturbed field equations | $G_{\rm eff}$, slip, lensing |
| **GR-06** | a tetrad | torsion, the TEGR identity, $f(\mathcal{T})$ cosmology |
| **GR-07** | a metric and a flat connection | non-metricity, $f(Q)$ cosmology |
| **GR-08** | quadratic actions | kinetic matrices, wave speeds, strong coupling |
| **GR-09** | a quadratic action | constraints, first and second class, a mode count |
| **GR-10** | the same perturbed action | the exact linear system, and what quasi-static costs |
| **GR-11** | a transverse traceless mode | the wave equation, $c_{\rm GW}$, and the siren distance |
| **GR-12** | a dense body in a background | the thin shell, PPN $\gamma$, and the bound on $|f_{R0}|$ |
| **GR-13** | GR-08's quadratic actions | the tensor coefficient, and which families sirens can see |

Conventions: signature $(-,+,+,+)$, geometrized units $G = c = 1$, $\kappa = 8\pi G/c^4$, and
$R^{\rho}{}_{\sigma\mu\nu}$ with the first index up. The connection is Levi-Civita in GR-01
through GR-05; GR-06 replaces it with a flat metric-compatible connection carrying torsion, and
GR-07 and GR-09 with a flat torsion-free connection carrying non-metricity. Those three build
their connection explicitly rather than inheriting it.

## The notebooks

### `GR-01-Curvature-Tensors`

Metric tensor in, curvature out:

- Christoffel symbols $\Gamma^{\lambda}{}_{\mu\nu}$ (printed once per independent symbol,
  since the lower pair is symmetric)
- Riemann tensor, shown as $R_{\rho\sigma\mu\nu}$ with one representative per symmetry orbit
- Ricci tensor, Ricci scalar, Einstein tensor
- Kretschmann scalar, plus vacuum and flatness checks

A metric library sits in a comment block ready to copy: Minkowski (Cartesian and spherical),
Schwarzschild, Reissner–Nordström, Schwarzschild–de Sitter, Kerr in Boyer–Lindquist
coordinates, FLRW flat and with spatial curvature, de Sitter static patch, anti-de Sitter in
the Poincaré patch, and the round 2-sphere.

### `GR-02-Field-Equations-From-Action`

Action in, field equations and cosmology out. The action need not be Einstein–Hilbert:
$f(R)$, $f(R,T)$, $f(\mathcal{T})$, $f(Q)$ and Gauss–Bonnet $f(R,\mathcal{G})$ are all
handled, through two engines.

|  | **Part I** covariant | **Part II** FLRW minisuperspace |
|---|---|---|
| varies | the metric $g_{\mu\nu}$ | the lapse $N(t)$ and scale factor $a(t)$ |
| valid for | any metric ansatz | FLRW only |
| handles | $R$, $f(R)$, $f(R,T)$ | $f(R)$, $f(\mathcal{T})$, $f(Q)$, $f(R,\mathcal{G})$ |
| returns | $\mathcal{E}_{\mu\nu} = 0$ | Friedmann I and II, and $E(z) = H(z)/H_0$ |

There are two engines because $f(R)$ and $f(R,T)$ live on a Levi-Civita connection and can be
varied covariantly, while $f(\mathcal{T})$ and $f(Q)$ have a tetrad or a flat non-metric
connection as their fundamental variable. Both routes derive rather than quote, and the
notebook checks that they agree on the overlap.

The dimensionless Hubble parameter is read off Friedmann I by writing $H = H_0 E(z)$ with
$1+z = 1/a$ and converting time derivatives with $d/dt = -(1+z)H_0E\,d/dz$. It comes out in
closed form for GR, $\Lambda$CDM, $f(\mathcal{T})$ and $f(Q)$, and as an ODE for $f(R)$ and
$f(R,\mathcal{G})$, which carry $\dot H$ and $\ddot H$; the engine reports which case you are
in and hands the ODE to `NDSolve` when it has to.

### `GR-03-Cosmological-Distances`

$E(z)$ in, observables out — the step that makes a model falsifiable. Everything here is an
integral of $E(z)$:

- comoving, transverse comoving, angular diameter and luminosity distances, and the distance
  modulus $\mu(z)$
- lookback time, age at redshift $z$, age of the universe
- particle and event horizons, comoving volume and its element
- the BAO sound horizon $r_s$

Distances come in Mpc and times in Gyr. A cosmology is just an $E(z)$ plus $h$ and $\Omega_k$,
so models sit side by side: the library covers flat, open and closed $\Lambda$CDM,
Einstein–de Sitter, de Sitter, $w$CDM and CPL, and shows how to carry a model over from GR-02
— including one whose $E(z)$ is only defined implicitly, which `ENumeric` solves by
root-finding at each redshift.

Reproduces, from scratch: age $13.796$ Gyr for Planck parameters (matching the closed form
$t_0 = \tfrac{2}{3H_0\sqrt{\Omega_\Lambda}}\operatorname{arcsinh}\sqrt{\Omega_\Lambda/\Omega_m}$),
the $D_A$ turnover near $z \approx 1.6$, Etherington reciprocity $D_L = (1+z)^2 D_A$, a
comoving distance to last scattering of $13.86$ Gpc, and $r_s \approx 147$ Mpc.

### `GR-04-Structure-Growth`

Distances cannot separate a modified-gravity model tuned to mimic the $\Lambda$CDM expansion
history. Growth can. This notebook solves the linear growth equation

$$\delta'' + \left[\frac{E'}{E} - \frac{1}{1+z}\right]\delta' = \frac{3\Omega_{m0}(1+z)}{2E^2}\frac{G_{\rm eff}}{G}\delta$$

and returns the growth factor $D(z)$, the growth rate $f(z) = d\ln\delta/d\ln a$, and the
redshift-survey observable $f\sigma_8(z)$, plus the growth index $\gamma$ read off rather
than assumed.

The background $E(z)$ and $G_{\rm eff}/G$ are chosen independently, which is the point: hold
the background fixed and vary only the gravity, and every distance in GR-03 stays identical
while $f\sigma_8$ moves by several percent. The library covers GR, the phenomenological
$\mu(z)$ used in DESI and Euclid forecasts, $f(\mathcal{T})$ and $f(Q)$, and the $f(R)$
large-scale, small-scale and general scale-dependent forms.

Reproduces $f(0) = 0.527$ against $\Omega_m^{0.55} = 0.530$, $\gamma(0) = 0.554$,
$f\sigma_8(0) = 0.427$, the exact $\Lambda$CDM growing-mode quadrature to two parts in $10^4$,
and $f = 1$ exactly in Einstein–de Sitter.

### `GR-05-Effective-Gravitational-Coupling`

GR-04 takes $G_{\rm eff}/G$ as input. This notebook derives it: it perturbs the field
equations in Newtonian gauge, expands in a bookkeeping parameter, applies the quasi-static
sub-horizon approximation, and solves for the metric potentials.

For $f(R)$, with $m = \dfrac{k^2}{a^2}\dfrac{f_{RR}}{f_R}$, it returns

$$\frac{G_{\rm eff}}{G} = \frac{1}{f_R}\frac{1+4m}{1+3m}, \qquad \eta = \frac{\Phi}{\Psi} = \frac{1+2m}{1+4m}, \qquad \Sigma = \frac{1}{f_R}.$$

The third is the one worth noticing, and it was not put in by hand: the scalaron cancels out
of the lensing combination, so light is deflected as in general relativity up to a constant
rescaling, with **no scale dependence at all**. Growth feels $m$; lensing does not.

Checks: the $\epsilon^0$ order reproduces the GR-02 Friedmann constraint for both $f = R$ and
$f = R - 2\Lambda$ — an independent cross-check between covariant perturbation theory here and
minisuperspace variation there — and the limits $m \to 0$, $m \to \infty$ give $1/f_R$ and
$4/(3f_R)$ with $\eta \to 1$ and $1/2$. The derived $G_{\rm eff}(k,z)$ is then fed back into
the GR-04 growth equation to show the scale-dependent growth that is the $f(R)$ fingerprint.

Scope: this is metric perturbation theory, so it covers $f(R)$ and any other theory whose
variable is $g_{\mu\nu}$. $f(\mathcal{T})$ and $f(Q)$ need tetrad and connection perturbations
and are **not** covered — those entries in GR-04's library stay quoted.

### `GR-06-Teleparallel-Geometry`

The first notebook in the series whose fundamental variable is not the metric. Teleparallel
gravity trades curvature for **torsion**: a flat, metric-compatible connection built from a
tetrad $e^A{}_\mu$, with

$$T^\lambda{}_{\mu\nu} = e_A{}^\lambda\left(\partial_\mu e^A{}_\nu - \partial_\nu e^A{}_\mu\right), \qquad \mathcal{T} = S_\lambda{}^{\mu\nu}T^\lambda{}_{\mu\nu}.$$

Everything is built from its definition — torsion, contortion, superpotential, torsion scalar
— and the torsion scalar is computed twice by routes sharing nothing but the torsion tensor,
because index placement in the superpotential is where such a calculation quietly goes wrong.

What it establishes:

- the TEGR identity $R = -\mathcal{T} + B$, verified on three unrelated geometries, which is
  why teleparallel gravity is equivalent to general relativity at linear order and stops being
  equivalent the moment $f$ is non-linear;
- the **good-tetrad problem as arithmetic**: two tetrads related by a local Lorentz rotation
  build the same flat Minkowski metric and have different torsion scalars, $-2/r^2$ and $0$;
- the $f(\mathcal{T})$ Friedmann equations, $\kappa\rho = \tfrac{1}{2}(f - 2\mathcal{T}f_\mathcal{T})$,
  derived by varying the tetrad and agreeing with what GR-02 got from the metric side.

It then derives the effective gravitational coupling, by varying **all sixteen** components of
a perturbed tetrad in a second order action:

$$\frac{G_{\rm eff}}{G} = -\frac{1}{f_{\mathcal{T}}}, \qquad \eta = \frac{\Phi}{\Psi} = 1,$$

scale independent and with no gravitational slip — sharply unlike $f(R)$, where GR-05 finds a
factor of $4/3$ running between large and small scales and a slip that goes to $1/2$. In the
sign convention of GR-02 and GR-04 this reads $1/f_{\texttt{Ts}}$, which is exactly what GR-04's
library had been quoting. It is quoted no longer.

Two wrong turns are kept in the notebook because each is worth more than the result. Section 8
shows a **diagonal** tetrad losing general relativity: $\delta\mathcal{T}$ carries no spatial
gradient and the restricted variation returns zero even in the TEGR limit, where the answer has
to be the Poisson equation. Section 9 shows a subtler trap — the obvious Wolfram spelling of
the field-theory Euler–Lagrange operator silently discards the $\partial_x$ term, which is
precisely the one carrying $k^2$, so a Lagrangian full of correct physics yields equations with
no Poisson term at all. A toy Lagrangian whose answer can be written by hand catches it in one
line.

The TEGR limit is the acceptance test throughout, and it is demanding: the $00$ equation must
carry $k^2$, the transverse equation must force $\Phi = \Psi$, the two momentum equations must
coincide, and both Lorentz modes must disappear.

### `GR-07-Symmetric-Teleparallel`

The third corner of the geometric trinity, and the last gap in the series. Symmetric
teleparallel gravity drops both curvature and torsion; what is left is **non-metricity**,
$Q_{\lambda\mu\nu} = \nabla_\lambda g_{\mu\nu}$, and its scalar $Q$, which differs from $R$ by
a total derivative — so STEGR is general relativity again, and stops being so as soon as $f$
is non-linear.

The awkward part, and the reason GR-05 and GR-06 both deferred it, is that the connection is a
**second field**. Being flat and torsion-free means it can be written from four Stückelberg
functions $\xi^\rho$, and setting $\xi^\rho = x^\rho$ is the coincident gauge that makes $Q$
come out as $-6H^2$. But a diffeomorphism shifts $\delta\xi^\mu$, so the coincident gauge
competes for the same freedom Newtonian gauge wants. This notebook spends it on Newtonian
gauge, keeps $\delta\xi^\mu$ as physical fields, and varies fourteen functions — ten metric,
four connection. The result:

$$\frac{G_{\rm eff}}{G} = \frac{1}{f_Q}, \qquad \eta = \frac{\Phi}{\Psi} = 1,$$

scale independent. Matter couples to the metric and not to the connection, so the connection
equations carry no source at all; they are pure constraint, and they are what closes the
system. Leave out the $xx$ equation and a connection mode stays undetermined and contaminates
the answer — a result that still depends on one has not closed.

With this the three families line up:

| | $G_{\rm eff}/G$ | slip | scale dependent |
|---|---|---|---|
| $f(R)$ | $\frac{1}{f_R}\frac{1+4m}{1+3m}$ | $\frac{1+2m}{1+4m}$ | **yes**, running to $4/3$ |
| $f(\mathcal{T})$ | $-1/f_{\mathcal{T}}$ | $1$ | no |
| $f(Q)$ | $1/f_Q$ | $1$ | no |

$f(R)$ carries a scalaron with a Compton wavelength, so it has a scale to compare $k$ against.
Neither teleparallel family does, and their modification is a pure rescaling.

### `GR-08-Strong-Coupling`

Every earlier notebook ends with the same caveat — a correct linear calculation is necessary,
not sufficient. This one makes that quantitative. It reuses the quadratic actions GR-06 and
GR-07 built and asks a different question of them: not what the equations say, but whether the
modes in them propagate at all.

The diagnostic is the kinetic matrix $K_{ij} = \partial^2 L^{(2)}/\partial\dot q_i\partial\dot q_j$
for one Fourier mode. A zero eigenvalue is a direction with no kinetic term; the ratio of the
gradient coefficient to the kinetic one is $c^2$.

- **Both theories propagate gravitational waves at exactly the speed of light**, $c_{\rm GW}^2 = 1$
  identically rather than approximately, both polarisations, every $f$. Neither is touched by
  the GW170817 constraint.
- **In $f(\mathcal{T})$ the extra Lorentz modes carry no kinetic term at all** — the boost and
  rotation rows of $K$ vanish identically, exactly like the lapse, which is a Lagrange
  multiplier. Around flat FLRW they do not propagate, which is strong coupling computed rather
  than asserted, and the surviving scalar entry is $6a^3(f_{\mathcal T} + 2\mathcal{T}f_{\mathcal{TT}})$.
- **In $f(Q)$ the bare kinetic matrix cannot settle the question**, and the notebook says so.
  It has full rank — but it has full rank at STEGR too, where general relativity has no
  propagating scalars, because $\sqrt{-g}Q$ keeps the time derivatives of the lapse that the
  Einstein–Hilbert boundary term removes. Settling it needs a Hamiltonian constraint analysis
  the notebook does not do — GR-09 does it.

Both families share a degeneracy locus, $f_X + 2Xf_{XX} = 0$ with $X = \mathcal{T}$ or $Q$. For
a concrete model that is a redshift: quadratic $f = X + \alpha X^2$ with $\alpha H_0^2 = -0.02$
degenerates at $z = 0.31$, in the middle of the range surveys measure. A model whose degeneracy
sits inside your fitting range is one whose linear predictions there should not be trusted.

Method note: real perturbations with $\cos kx$ and $\sin kx$ and an average over one wavelength,
not a complex plane wave — a quadratic action built from $e^{ikx}$ picks up $e^{2ikx}$ and the
phases stop cancelling.

### `GR-09-Hamiltonian-Constraints`

GR-08 left one question open and named the tool needed to close it. This notebook is that tool:
an Ostrogradsky reduction followed by the full Dirac algorithm — velocity Hessian, primary
constraints, canonical Hamiltonian, consistency conditions, secondary constraints, and the
split into first and second class from the rank of the matrix of Poisson brackets.

**General $f(Q)$ propagates one scalar mode where general relativity propagates none.** The
mechanism is in the algebra rather than in the counting: the number of constraints is the same
twelve in both cases, but two of them move from first class to second once $f_{QQ}
eq 0$. A
first-class constraint costs two phase-space dimensions and a second-class one costs one, so
moving two across the line frees exactly one degree of freedom. General relativity's scalar
sector is pure gauge; $f_{QQ}$ breaks one of those gauge symmetries and the mode that was gauge
becomes physical.

Two implementation traps are worth repeating, because both produced confident wrong answers
before an acceptance test caught them.

- **Primary constraints come from the null space of the velocity Hessian, not from momenta that
  happen to be velocity free.** One cross term $\dot q_1\dot q_2$ makes every momentum depend on
  a velocity while the Hessian stays degenerate, and a scan over momenta then misses the
  constraint entirely.
- **A constraint is new only when it is linearly independent of the ones already found.** The
  chain regenerates old constraints rescaled by background factors; an equality test never
  recognises them and the algorithm runs forever.

And one modelling trap, which is the reason this could not be done with the notebooks that
already existed: **fixing Newtonian gauge in the action before varying deletes the momentum
constraint**, since that constraint is what varying the shift produces. Do the count that way
and general relativity comes out with one scalar mode instead of none — off by exactly one,
because exactly one constraint was thrown away. The scalar sector here is therefore kept
complete, lapse and shift and curvature and anisotropy, with no gauge fixed, and the metric part
is written in ADM form so that lapse and shift carry no time derivatives by construction.

Three acceptance tests frame the result, and all three must return zero: general relativity in
the metric ADM sector, STEGR as $f = Q - 2\Lambda$, and the linear branch $f'' = 0$. The count
is also reproduced with concrete numbers, so it does not rest on a symbolic rank.

**$f(R)$ gets the same treatment, by a second route.** It propagates **one** scalar mode, the
scalaron, where general relativity propagates none — and here the constraint structure is
*identical* between the two, four first class and eight second class in both. Only the number of
coordinates differs, nine against eight, the extra one being the scalaron. So $f(R)$ adds a
field, while $f(Q)$ moves a constraint across the first/second-class line; two different ways of
freeing exactly one mode.

$f(R)$ cannot be reached the ADM way, and the reason is the third trap in this notebook. In ADM
variables ${}^{(4)}R = K_{ij}K^{ij} - K^2 + {}^{(3)}R$ **plus total derivatives**, and those are
droppable in general relativity, where they multiply a constant, but not in $f(R)$, where they
multiply $f'$. So nothing is integrated by parts by hand: the full four-metric is rebuilt out of
the ADM pieces, ${}^{(4)}R$ is computed from it directly, and the second time derivatives go to
`OstrogradskyReduce` intact. The eight second-class constraints are then a check on the route
itself — they are what removes the Ostrogradsky variables again, the algorithm discovering on its
own that general relativity is not fourth order. A mistake in the four-metric would show up as a
control that no longer returns zero.

Scope: de Sitter background, scalar sector, one Fourier mode, linearised theory. Constraint
structure can be background dependent, so a background with $\dot H 
eq 0$ is not covered.

### `GR-10-Beyond-Quasi-Static`

Every $G_{\rm eff}$ in the series rests on two approximations GR-05 states plainly and then
never tests. This notebook drops both and solves the linear scalar system exactly.

Two things had to be added first. Matter gets a **velocity** — GR-05 leaves it comoving, which
is free under the quasi-static approximation and fatal without it, because the continuity
equation needs somewhere for the density to flow to. And the trace of the $ij$ equations turns
out to be **fourth order in $\Phi$**: that is the scalaron hiding inside the metric variables.
Promoting $\delta R$ to a field of its own drops the order, leaving five first-order equations
in $(\Phi, \delta R, \dot{\delta R}, \delta, v)$. The step that makes it work is an identity,

$$\mathcal{E}_{xx} - \mathcal{E}_{yy} = -k^2\left[f_R(\Phi - \Psi) - f_{RR}\delta R\right],$$

so the gravitational slip is sourced by the scalaron and by nothing else.

**The headline is that the obvious comparison is the wrong one.** Setting the exact $\mu$
against GR-05's formula gives a $15\%$ disagreement at $k = 3aH$, which looks like the
quasi-static approximation failing badly. It is not. $\mu \equiv -2k^2\Psi/(a^2\kappa\rho\delta)$
is not $1$ in general relativity either — the exact $00$ equation carries
$3H(\dot\Phi + H\Psi)$ — and almost all of that $15\%$ is a relativistic correction $\Lambda$CDM
has too, which GR-04's growth equation drops in exactly the same way. Run general relativity
through the same solver and divide it out, and what is left is the part that is really about
$f(R)$:

| | $\mu_{f(R)}/\mu_{\rm QS} - 1$ | the same, GR divided out | as a fraction of the modification |
|---|---|---|---|
| $k = 3aH$ | $-14.7\%$ | $3.6\times10^{-6}$ | $1.2\%$ |
| $k = 10aH$ | $-1.5\%$ | $-9.8\times10^{-6}$ | $0.44\%$ |
| $k = 50aH$ | $-0.06\%$ | $6.8\times10^{-6}$ | $0.015\%$ |

at $|f_{R0}| = 10^{-4}$ on a $\Lambda$CDM background. The third column is the honest one: the
modification itself is only a few parts in $10^4$ there, so an error that is negligible against
$\mu$ need not be negligible against the thing being measured. Repeated at $|f_{R0}| = 10^{-2}$,
where the modification is $O(1)$ and the ratio is unambiguous, it reads $3.7\%$ of the
modification at $k = 3aH$, $2.0\%$ at $k = 10aH$, $0.74\%$ at $k = 20aH$ and $0.12\%$ at
$k = 50aH$.

So GR-05's $\mu$, used as a **ratio to general relativity** — which is how GR-04 and CosmoFit
use it — holds to a few percent of the modification right down to $k \approx 3aH$, and to about
one percent there for the realistic $|f_{R0}| = 10^{-4}$: far further than "sub-horizon"
suggests. What quasi-static does *not* capture is the near-horizon
correction to the Poisson equation, and that is a property of the growth equation rather than
of $G_{\rm eff}$.

Two traps are kept, because both produce confident wrong answers. **$f_R$ is a function of
time**, through the background curvature; abbreviate it to an inert symbol and then
differentiate, and every chain-rule term is silently dropped — the notebook runs the reduction
both ways and shows the trace equations differ. And the sign in the slip relation is easy to
invert, which leaves a system whose derivative structure looks perfectly healthy and whose
numbers are wrong.

Scope: metric perturbation theory, so $f(R)$ and nothing teleparallel; pressureless matter, no
radiation; one Fourier mode; the $\Lambda$CDM background is exact only to $O(|f_{R0}|)$, which
the notebook measures rather than assumes.

### `GR-11-Gravitational-Wave-Propagation`

GR-08 asked how fast the tensor modes travel and found $c_{\rm GW}^2 = 1$ for the two
teleparallel families. Speed is half the story: a wave is also **damped** on the way here, and
if it is damped by more than the expansion alone then a standard siren reports the wrong
distance. This notebook does the tensor sector of $f(R)$, which GR-08 does not touch at all.

The transverse traceless mode is the easy corner of perturbation theory, and for a good reason:
it does not change the Ricci scalar, so $\delta R = 0$ and the scalaron that made GR-10 fourth
order is simply not excited. Through the same engine GR-05 and GR-10 use, the result is

$$\ddot h + \left(3H + \frac{\dot f_R}{f_R}\right)\dot h + \frac{k^2}{a^2}h = 0,$$

so $c_{\rm GW}^2 = 1$ identically for $f(R)$ too — the gradient and acceleration terms carry the
same $f_R$ — while the friction is $(3 + \alpha_M)H$ with $\alpha_M = d\ln f_R/d\ln a$. One trap
on the way: the $34$ equation at first order still contains (background $\mathcal{E}_{yy}$)
$\times\,h$, which looks exactly like matter sourcing a tensor mode until the background
equations are imposed and it disappears.

Reducing in conformal time shows the wave carries $a\sqrt{f_R}h$ unchanged rather than $ah$ —
$h = u/(a\sqrt{f_R})$ removes the friction exactly, leaving $u'' + (k^2 - z''/z)u = 0$ — and the
notebook then confirms that numerically by following a mode's envelope: multiplied by $a$ it
drifts by $5\%$, multiplied by $a\sqrt{f_R}$ by $0.3\%$. Hence

$$\frac{d_L^{\rm GW}(z)}{d_L^{\rm EM}(z)} = \sqrt{\frac{f_R(0)}{f_R(z)}},$$

which falls with redshift and **saturates**, because there is nothing left to accumulate once
$f_R(z)$ has reached one. The whole effect is bounded at every redshift at once by
$1 - \sqrt{1 - |f_{R0}|} \approx |f_{R0}|/2$.

**That is a clean negative result.** Three observables in the series now track the same single
function: $G_{\rm eff}/G = (1/f_R)(1+4m)/(1+3m)$ for growth, $\Sigma = 1/f_R$ for lensing, and
now the siren ratio. At $z = 0$ lensing gives $\Sigma - 1 \approx |f_{R0}|$ and sirens give
$|f_{R0}|/2$ — lensing wins by a factor of two on exactly the same parameter — and the Solar
System pins that parameter below $10^{-6}$, which GR-12 derives. So standard sirens cannot
constrain viable $f(R)$, and not by bad luck: the $f_R \approx 1$ that lets the theory survive a laboratory is
what forces the wave to travel as it does in general relativity.

It does **not** follow for the other two families. GR-08's tensor kinetic coefficients are
$f_{\mathcal{T}}$ and $f_Q$, the same coefficients that appear in $G_{\rm eff}/G$ there, and
those are not tied to unity the way $f_R$ is. Whether sirens can see them is settled in
GR-13, which loads GR-08's quadratic actions rather than rebuilding the machinery.

Scope: metric perturbation theory, so $f(R)$ only. The propagation effect alone — a full siren
prediction also has to ask whether the source's own emission is modified, which is a question
about screening near the binary rather than about the wave.

### `GR-12-Screening-And-The-Solar-System`

Every other result in the series is derived. One was not: $|f_{R0}| \lesssim 10^{-6}$, the Solar
System bound, which GR-05, GR-11, this README and the companion fitting library all leaned on
and all took on trust. GR-11 made that worse by using it twice. This notebook derives it.

The problem is already in GR-05. Its slip $\eta = (1+2m)/(1+4m)$ goes to $\tfrac{1}{2}$ on small
scales, and in a static weak field that slip *is* the post-Newtonian $\gamma$. Cassini measures
$\gamma - 1 = (2.1 \pm 2.3)\times10^{-5}$, so $f(R)$ with a light scalaron is not marginally
disfavoured — it is out by four orders of magnitude, and needs a mechanism rather than a fit.

The mechanism is the whole content of one relation. Linearising the trace equation gives a
Poisson equation for the scalaron with $-\kappa/3$ where gravity has $+\kappa/2$, so sourced by
the same matter $\delta f_R = -\tfrac{2}{3}\Phi_N$: **the field cannot move further than
$\tfrac{2}{3}|\Phi_N|$ however dense the body gets.** If the background value it would have to
climb from is larger than that, it never reaches its interior minimum and the body is bare.

Working that out for a uniform sphere with a pinned core gives a closed form, all of it checked
rather than quoted — with $\varepsilon = |f_{R,\rm bg}|/\Phi_N$,

$$w = \frac{(x-x_s)^2(x+2x_s)}{3\varepsilon x}, \qquad \frac{x_s}{r_b} = \sqrt{1-\varepsilon},
\qquad \frac{A}{A_{\rm linear}} = 1 - (1-\varepsilon)^{3/2},$$

so a screened core exists only for $\varepsilon < 1$, and deep in the screened regime the fifth
force is suppressed by $\tfrac{3}{2}\varepsilon$ — linear in how small the background field is.
Requiring the **Galaxy** to have a thin shell, which is what actually screens the Solar System
sitting inside it, is then $|f_{R0}| < \Phi_{\rm gal} \approx 10^{-6}$. The bound looks like a
potential because it is one: the depth of the well the scalaron has to climb out of.

The pinned core is an idealisation, and section 6 measures it instead of waving at it. What was
dropped is a fraction $1/(D\sqrt{w})$ of the source with $D$ the density contrast, so it matters
only in a boundary layer where $w \lesssim 1/D^2$ — reaching $\sqrt{\varepsilon}/D$ out from the
core, which for a real body is a part in $10^{27}$ of the shell. That is worth doing
algebraically: the interior minimum is an *equilibrium*, so a numerical solution started exactly
on it never leaves, and started just off it, where it leaves depends on rounding.

One caveat the notebook keeps in view, because it cuts the other way: screening is **local**. It
suppresses the fifth force near dense bodies and does nothing in the voids and filaments where
GR-04 and GR-05 work. A model can be invisible in the Solar System and still move $f\sigma_8$ by
percent, which is why the cosmological fits are worth doing rather than foreclosed by the
laboratory. The bound constrains $|f_{R0}|$, not the shape of $f$.

### `GR-13-Sirens-And-The-Other-Two`

GR-11 answered the siren question for $f(R)$ and ended by saying flatly that the answer does not
carry over. This closes that, and the answer is the opposite.

It loads GR-08 rather than rebuilding it. GR-08 already constructed the second-order Lagrangians
for both teleparallel families and checked them against their general-relativity limits; it then
asked those objects about the propagation *speed*. This asks about the **size of the kinetic
coefficient**, whose running is the friction, and a second copy of that index algebra would be a
second thing to keep correct.

$$\text{$f(\mathcal{T})$}: \;\; \frac{\partial^2L}{\partial\dot h^2} = -\tfrac{1}{2}f_\mathcal{T}a^3,
\qquad
\text{$f(Q)$}: \;\; \frac{\partial^2L}{\partial\dot h^2} = +\tfrac{1}{2}f_Q a^3,$$

with $c_{\rm GW}^2 = 1$ falling out again as a by-product — GR-08's own result, reached from the
coefficient instead of the ratio. Normalised so general relativity gives one, $M^2 =
-f_\mathcal{T}$ and $M^2 = f_Q$, the two sign conventions cancelling exactly as they should. So
all three families say GR-11's sentence with a different noun:
$d_L^{\rm GW}/d_L^{\rm EM} = \sqrt{M^2(0)/M^2(z)}$.

**What differs is how large the noun may be, and there they part company completely.**

| | largest $1 - d_L^{\rm GW}/d_L^{\rm EM}$ |
|---|---|
| $f(R)$, with $\|f_{R0}\| < 10^{-6}$ | $5\times10^{-7}$ |
| $f(\mathcal{T})$ or $f(Q)$, $b = 0.1$ | $4\times10^{-2}$ |
| $f(\mathcal{T})$ or $f(Q)$, $b = 0.2$ | $1\times10^{-1}$ |

Five orders of magnitude, and not because the teleparallel families are wilder. The chain that
catches $f(R)$ starts with GR-05's slip running to $\eta = 1/2$ on small scales, which Cassini
excludes, which forces the screening of GR-12, which forces $|f_{R0}| < 10^{-6}$ — the very
number that then bounds the siren effect. GR-06 and GR-07 find $\eta = 1$ *identically* for the
other two: that chain never starts.

**One constraint the notebook raises against its own result, and cannot settle.** Nothing local
pins the *value* of $f_\mathcal{T}$; something does pin its *rate*. GR-06 gives
$G_{\rm eff}/G = -1/f_\mathcal{T} = 1/M^2$, so a running $M^2$ is a running gravitational
constant, and lunar laser ranging allows
$|\alpha_M| = |d\ln M^2/d\ln a| \lesssim 2\times10^{-3}$. The model above has
$\alpha_M = 0.16$ at $b = 0.2$ — seventy times that, and fifteen times it even at $b = 0.05$.

Whether that excludes it is a question the series cannot answer yet, and the reason is exactly
GR-12: lunar ranging measures $G$ *inside the Solar System*, so what enters is the local
$G_{\rm eff}$, and for $f(R)$ the local and cosmological values differ by many orders of
magnitude precisely because the chameleon pins the field. Whether $f(\mathcal{T})$ or $f(Q)$ do
anything analogous needs the static spherically symmetric solution in those theories — the
teleparallel counterpart of GR-12, which GR-06 and GR-07 do not build. So the ten percent stands
only if the cosmological running survives into the local system, and if it does, the same running
is in trouble with a much older experiment. Either way both are measuring the same function.

Two smaller results. The background of the fitted power law is rebuilt here from GR-06's own
Friedmann equation, independently of GR-02's `HubbleFunction`, and agrees with it to twelve
digits — which is also a check on `cosmofit-export.wls`, since those are the numbers it ships.
And $f(\mathcal{T})$ and $f(Q)$ power laws give the *same* $M^2$ as well as the same background,
so sirens cannot separate those two families either, only both from general relativity.

Finally, one value of the exponent is worth knowing about: at $b = \tfrac12$ the correction
carries $(1-2b)$ in $f - 2\mathcal{T}f_\mathcal{T}$ and leaves the Friedmann constraint entirely,
so $E(0) = 1$ has no solution and the closure fails. What makes it worth saying is that $M^2$
still depends on $A_0$ there. A fit near that exponent would have a parameter that does nothing
to the expansion history and everything to the waves.

## Conventions worth knowing before you trust the output

- **Torsion and non-metricity scalars.** In flat FLRW the notebook uses
  $\mathcal{T} = Q = -6H^2$, chosen so that the linear theory is ordinary GR: $f = \mathcal{T}$
  is TEGR and $f = Q$ is STEGR, exactly as $f = R$ is GR. Papers that define
  $\mathcal{T} = +6H^2$ write TEGR as $f = -\mathcal{T}$; flip the sign in `$flrwTorsion` and
  read your $f$ unchanged. The notebook verifies the GR limit rather than asserting it.
- **The matter Lagrangian in $f(R,T)$.** Both $\mathcal{L}_m = -p$ and $\mathcal{L}_m = -\rho$
  appear in the literature. They give the same $T_{\mu\nu}$ but different $\Theta_{\mu\nu}$,
  so they are genuinely different theories. `FieldEquations` refuses to guess and asks for it
  explicitly whenever the action depends on $T$.
- **Curvature backend.** Wolfram has no built-in Riemann/Ricci/Einstein functions, so the
  notebooks use the Function Repository resources `MetricTensor`, `ChristoffelSymbols`,
  `RiemannTensor`, `RicciTensor` and `EinsteinTensor`. These download once and are cached, so
  the first evaluation of the setup cell needs a network connection and later ones do not.
  The Gauss–Bonnet term is assembled explicitly as
  $R^2 - 4R_{\mu\nu}R^{\mu\nu} + R_{\rho\sigma\mu\nu}R^{\rho\sigma\mu\nu}$ rather than taken
  from the resource's `"EulerScalar"` property, which is a differently normalised object.
- **Nothing is quoted.** GR-04 takes $G_{\rm eff}/G$ as input, and every entry its library
  offers is derived inside the series: $f(R)$ in GR-05, $f(\mathcal{T})$ in GR-06, $f(Q)$ in
  GR-07. What remains are the quasi-static and sub-horizon approximations, shared by all three.
  For $f(R)$ they are no longer untested either: GR-10 solves the same system exactly and puts
  a number on them. The last quoted quantity anywhere in the series was the Solar System bound
  $|f_{R0}| \lesssim 10^{-6}$, and GR-12 derives that too. What is quoted now is measurement —
  Cassini's $\gamma$, the Galaxy's potential — which is the only kind of input that should be.
- **Every $f(R)$ deviation in this series is one function.** $G_{\rm eff}$, the lensing
  $\Sigma$, the gravitational-wave friction and the siren distance are all built from $f_R$ and
  its running, so the bound on $|f_{R0}|$ caps all of them at once. GR-05 derives the first two,
  GR-11 the last two, and GR-12 the bound — which was the last quoted number in the series.
- **The $G_{\rm eff}$ results are metric-sector statements.** GR-08 shows that in
  $f(\mathcal{T})$ the extra Lorentz modes have identically zero kinetic terms around flat
  FLRW, so linear theory does not describe them, and that for $f(Q)$ the status of the
  connection modes is open at that level of analysis. This does not touch the algebra of GR-06
  and GR-07; it bounds how far to trust it. $f(R)$ is unaffected — its scalaron has a healthy
  kinetic term.
- **Two sign conventions for the torsion scalar, both used here on purpose.** GR-06 computes
  the standard superpotential contraction, $\mathcal{T} = +6H^2$ in flat FLRW, so TEGR is
  $f = -\mathcal{T}$. GR-02 uses `Ts` $= -6H^2$ instead, chosen so a linear $f$ is GR, matching
  the convention in much of the $f(T)$ cosmology literature. They are related by
  `Ts` $= -\mathcal{T}$, and GR-06 checks that both give the same Friedmann equations.
- **GR-05's approximations are real approximations, and GR-10 says how real.** Used as a ratio
  to general relativity — the way GR-04 consumes it — $\mu$ is good to about a percent of the
  modification down to $k pprox 3aH$. Used as an absolute Poisson equation it is $15\%$ off
  there, but so is $\Lambda$CDM, for the same reason and by nearly the same amount. The
  sub-horizon truncation is an explicit, readable term filter rather than something hidden, so
  you can change it; GR-10 does exactly that, by not applying it. $f(\mathcal{T})$ and $f(Q)$
  are **not** covered — GR-10 is metric perturbation theory, like GR-05.
- **The CMB acoustic scale in GR-03 is indicative only.** $r_s$ is accurate, but $\theta_*$
  comes out around $100\theta_* = 1.06$ against Planck's $1.0411$, because $z_d$ and $z_*$ are
  put in by hand and $\Omega_r$ is a single number rather than a proper photon plus neutrino
  background. Precision CMB work needs a Boltzmann code.

## Building

Each notebook is authored as an annotated `.wl` source and compiled to the WLJS `.wln`
format by `wl2wln.wls`:

```
wolframscript -file wl2wln.wls GR-01-Curvature-Tensors.wl
```

The markup is two comment forms — `(*::md::` … `::*)` for a Markdown cell and `(*::code::*)`
for a Wolfram input cell. Everything outside a marker is ignored, so the `.wl` stays a valid
source file you can evaluate directly with `Get`, which is how the physics is tested before
it ever becomes a notebook.

Both the `.wl` sources and the built `.wln` notebooks are committed, so you can open a
notebook straight away or rebuild it from source. The committed notebooks also carry their
saved evaluation outputs, which is why the results are visible on GitHub without running
anything — and why a rebuild is not free: `wl2wln.wls` regenerates the notebook from source
and discards every one of them.

For a small change, `patchwln.wls` edits a built notebook in place instead:

```
wolframscript -file patchwln.wls GR-04-Structure-Growth.wln patch.txt
```

where `patch.txt` holds the old and new text between `%%%OLD%%%` and `%%%NEW%%%` marker
lines. A Markdown passage lives in a `.wln` exactly twice, once in the hidden `.md` input cell
and once in its rendered output cell, so the tool refuses to write unless it finds precisely
that many matches — a partial replacement would leave the source and the rendered copy
disagreeing. It also rejects replacement text containing a blank line followed by `%`, which
is how the reader recognises a cell separator, and it parses the result with WLJS's own reader
before and after so you can see the output cells survived. Edit the matching `.wl` too, or the
next rebuild reverts the change.

Neither tool looks at the physics. `check.wls` does — it evaluates the sources and reports every
verification cell's verdict, without opening a notebook. See [Checks](#checks).

## Requirements

- Wolfram Engine or Mathematica (developed against 15.0)
- [WLJS Notebook](https://jerryi.github.io/wljs-docs/) to open the `.wln` files
- A network connection the first time, for the Function Repository resources

One environment note, found the hard way. `D[expr, Derivative[0,1][f][t,x]]` &mdash;
differentiating with respect to a derivative of a two-variable function &mdash; does not give
the same answer in WLJS as it does under `wolframscript`, on the same Wolfram Engine 15.0. In
WLJS it comes back zero; under `wolframscript` it is correct. Nothing in these notebooks or in
the WLJS kernel sources accounts for the difference, and it is why every Euler-Lagrange
operator here works by substituting plain symbols for the derivatives before differentiating,
and why no check asserts what the naive spelling returns. GR-06 section 9 shows both side by
side.

## Checks

Every notebook carries its own verification cells, and the committed `.wln` files ship with
those cells already evaluated, so the results below are visible on GitHub without running
anything:

- Schwarzschild: $R = 0$, $G_{\mu\nu} = 0$, $K = 48M^2/r^6$, $R_{trtr} = -2M/r^3$
- Kerr: $R_{\mu\nu} = 0$
- Round 2-sphere: $R = 2/R_0^2$
- FLRW: $G_{tt} = 3\dot a^2/a^2$, and $3H^2/\kappa = \rho$, $p = -(3H^2 + 2\dot H)/\kappa$
- TEGR ($f = \mathcal{T}$) and STEGR ($f = Q$) reproduce GR exactly
- General relativity, STEGR and the linear branch $f'' = 0$ each give **zero** propagating
  scalar modes in the Dirac count; general $f(Q)$ gives **one**, and so does general $f(R)$ —
  reached by a second route, from the four-metric rather than from ADM, with general relativity
  again giving zero as the control
- Gauss–Bonnet is topological in four dimensions: $f = R + \alpha\mathcal{G}$ gives plain GR
- $f(R) = R + \alpha R^2$ agrees between the covariant and minisuperspace routes, and matches
  the textbook $3FH^2 = \kappa\rho + \tfrac{1}{2}(FR - f) - 3H\dot F$
- $\Lambda$CDM: $E^2 = \Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_k(1+z)^2 + \Omega_\Lambda$,
  $q_0 = -1 + \tfrac{3}{2}\Omega_m$, acceleration from $z = 0.671$ at $\Omega_m = 0.3$

### Running them headless

Those cells only report when someone opens a notebook and looks. `check.wls` runs them instead:

```
wolframscript -file check.wls
```

It reads each `.wl` source with the same cell reader `wl2wln.wls` uses, evaluates the code cells
in order, and collects every association carrying an `"ok"` key — wherever it sits, inside a
`Dataset` or not. Failures are printed with the cell they came from, and the exit code is 0 only
when every check passed and no cell errored, so it can be wired to anything that reads exit
codes.

```
GR-06-Teleparallel-Geometry.wl
  15 cells, 45.7 s, 28 checks, all pass
```

The whole series is **214 checks**, and all of them pass. It takes anywhere from
seven to twenty minutes depending on how busy the machine is, most of that inside GR-09 and
GR-10; GR-12 runs in under a second. Flags: `--verbose` lists
every check rather than only the failures, `--parse-only` reads the sources and counts cells
without evaluating anything, and `--timeout=N` caps the seconds any one expression may take
(`0` removes the cap). Named files run instead of the whole series.

Two things are worth knowing before reading its output.

**GR-01 reports no checks, and that is true rather than a bug.** Its verification cells print
tensors for a reader to inspect — Schwarzschild's vanishing Ricci, the 2-sphere's $R = 2/R_0^2$ —
and never reduce them to a boolean, so there is nothing for a runner to collect. Those claims are
checked by eye.

**The label key is not uniform.** Most notebooks write `"check"`, GR-05 and GR-08 write
`"statement"`, and eleven rows write `"acceptance test"`. All three are read, rather than edited
into agreement: changing a notebook means rebuilding its `.wln`, and that discards the saved
outputs which are the reason the results are visible on GitHub at all. A row spelled some fourth
way is still reported, under whatever string it carries, so a new spelling shows up as a check
with an odd label instead of disappearing.

`check-self-test.wl` is the negative control, and it is meant to fail:

```
wolframscript -file check.wls check-self-test.wl
```

Every notebook in the series passes, so a green run on its own does not distinguish a working
checker from one that reports nothing. That file holds one of each way a check can go wrong —
`ok -> False`, an `"ok"` that never evaluated to a boolean at all, a row outside a `Dataset` under
the other spelling, and a cell that aborts — and the runner has to see all of them, and has to
keep going past the cell that died. Its name does not begin with `GR-`, so the default run leaves
it alone.

### Against a second implementation

The notebooks have a numerical counterpart: [CosmoFit](https://github.com/salihyesil59/CosmoFit),
a Python cosmological parameter-estimation library that fits these same models to BAO, supernova,
cosmic-chronometer and CMB data. It derives its Friedmann constraints independently, in sympy,
and solves them with its own root finder — no code is shared with these notebooks.

`cosmofit-reference.wls` prints $E(z)$ at 25 digits by loading GR-02 and driving its
`FriedmannEquations` and `HubbleFunction` directly, rather than re-deriving anything:

```
wolframscript -file cosmofit-reference.wls
```

Those numbers are pinned on the other side in `tests/test_notebook_agreement.py`. Two derivations
written years and languages apart currently agree to **3 parts in $10^{16}$** — machine precision
— for $\Lambda$CDM and for the $f(\mathcal{T})$ power law at $n = -0.5,\, 0.2,\, 0.7$, both
through CosmoFit's hand-written models and through its action compiler.

`cosmofit-mu-reference.wls` does the same for the **perturbations**, driving GR-06's own
`geffSubHorizon`. That comparison matters more than it sounds, because the two codebases
disagree about the sign of the torsion scalar deliberately: GR-06 uses $\mathcal{T} = +6H^2$ with
TEGR at $f = -\mathcal{T}$ and derives $G_{
m eff}/G = -1/f_{\mathcal{T}}$, while CosmoFit uses
$\mathcal{T} = -6H^2$ and $+1/f_{\mathcal{T}}$. Two minus signs that cancel — so the same theory
must give the same number, and would not if either were ever flipped alone. They agree to
**1.4 parts in $10^{16}$**.

The two also agree about where the model is *sick*: at $n = 0.7$ this expression gives
$G_{
m eff}/G = -4.44$, then $+25.6$ — negative, then through a pole — and CosmoFit refuses that
region rather than returning the numbers.

### The bridge the other way

Those two scripts compare the codebases on models someone has already written twice.
`cosmofit-export.wls` goes the other way — derive something new here, and get the Python that
fits it:

```
wolframscript -file cosmofit-export.wls "Ts + A0 (-Ts)^bb" A0 bb=0.2
```

It emits a runnable `Action(...).build(...)` **with this repository's own $E(z)$ carried along as
an assertion**, so the generated script fails loudly if the two derivations ever disagree. Both
built-in examples currently come out at $0.0\mathrm{e}0$.

Carrying the numbers is not decoration; it is what made the translation correct. Torsion needs
no translation — GR-02 and CosmoFit both put $\mathcal{T} = -6H^2$ — but non-metricity is not a
rename. GR-02 uses $Q = -6H^2$ and CosmoFit $Q = +6E^2$, and since **both** call a linear $f$ the
GR limit, their two STEGR Lagrangians differ by an overall sign that the $E(0)=1$ closure hides
for linear $f$ and for nothing else. What survives is

$$f_{\rm CosmoFit}(Q) = -f_{\rm GR02}(-Q),$$

sending $Q_s + A_0(-Q_s)^b$ to $Q - A_0Q^b$. That was measured, not reasoned: the literal
$Q_s \to Q$ produces a complex closure and $Q_s \to -Q$ alone is rejected as having no real
solution. Both were caught by the assertion rather than by inspection.

The script also refuses rather than guessing. It will not export a fourth-order $f(R)$, where
GR-02 returns an ODE and there is no closed form to assert against; and it will not export a
model whose closure pushes $E(z)$ off the real line. That second case is not hypothetical:
$f = Q + \alpha Q^2$ looks entirely reasonable, and once the closure fixes $\alpha$ its
discriminant goes negative at $z = 0.06$. Root-finding through it returns a confident $0.845$
that barely moves with redshift.

Worth knowing: GR-02 hands the $f(\mathcal{T})$ and $f(Q)$ power laws the *identical* background
constraint, since both scalars are $-6H^2$ on flat FLRW. That is correct rather than an oversight,
and it means background data alone cannot separate those two families — a growth or perturbation
observable is needed, which is what GR-04 and GR-05 are for.

## License

MIT — see [LICENSE](LICENSE). Use it, change it, publish with it; just keep the copyright
notice. If it ends up being useful in published work, a citation is appreciated but not
required.
