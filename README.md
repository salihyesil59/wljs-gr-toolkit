# General Relativity notebooks

Wolfram Language notebooks for symbolic general relativity, written for
[WLJS Notebook](https://jerryi.github.io/wljs-docs/). Give them a metric or an action and
they hand back the curvature tensors, the field equations, the Friedmann equations and the
dimensionless Hubble parameter — non-zero, independent components only.

Conventions throughout: signature $(-,+,+,+)$, geometrized units $G = c = 1$,
$\kappa = 8\pi G/c^4$, Levi-Civita connection, and $R^{\rho}{}_{\sigma\mu\nu}$ with the first
index up.

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

Scope: de Sitter background, scalar sector, one Fourier mode, linearised theory. Constraint
structure can be background dependent, so a background with $\dot H 
eq 0$ is not covered.

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
  GR-07. What remains are the quasi-static and sub-horizon approximations, shared by all three,
  which fail on scales approaching the horizon.
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
- **GR-05's approximations are real approximations.** Quasi-static and sub-horizon both fail
  on scales approaching the horizon. The sub-horizon truncation is implemented as an explicit,
  readable term filter rather than hidden, so you can change it.
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

## Requirements

- Wolfram Engine or Mathematica (developed against 15.0)
- [WLJS Notebook](https://jerryi.github.io/wljs-docs/) to open the `.wln` files
- A network connection the first time, for the Function Repository resources

## Checks

Both notebooks carry their own verification cells, and the results below are reproduced by
running them:

- Schwarzschild: $R = 0$, $G_{\mu\nu} = 0$, $K = 48M^2/r^6$, $R_{trtr} = -2M/r^3$
- Kerr: $R_{\mu\nu} = 0$
- Round 2-sphere: $R = 2/R_0^2$
- FLRW: $G_{tt} = 3\dot a^2/a^2$, and $3H^2/\kappa = \rho$, $p = -(3H^2 + 2\dot H)/\kappa$
- TEGR ($f = \mathcal{T}$) and STEGR ($f = Q$) reproduce GR exactly
- General relativity, STEGR and the linear branch $f'' = 0$ each give **zero** propagating
  scalar modes in the Dirac count; general $f(Q)$ gives **one**
- Gauss–Bonnet is topological in four dimensions: $f = R + \alpha\mathcal{G}$ gives plain GR
- $f(R) = R + \alpha R^2$ agrees between the covariant and minisuperspace routes, and matches
  the textbook $3FH^2 = \kappa\rho + \tfrac{1}{2}(FR - f) - 3H\dot F$
- $\Lambda$CDM: $E^2 = \Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_k(1+z)^2 + \Omega_\Lambda$,
  $q_0 = -1 + \tfrac{3}{2}\Omega_m$, acceleration from $z = 0.671$ at $\Omega_m = 0.3$

## License

MIT — see [LICENSE](LICENSE). Use it, change it, publish with it; just keep the copyright
notice. If it ends up being useful in published work, a citation is appreciated but not
required.
