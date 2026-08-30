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
- **$G_{\rm eff}$ in GR-04 is input, not output.** The growth equation is derived in the
  notebook; the $G_{\rm eff}/G$ expressions in its library are the standard quasi-static,
  sub-horizon results quoted from the literature. Deriving them means perturbing each
  theory's field equations, which the series does not yet do. Check the convention against
  your source.
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
notebook straight away or rebuild it from source.

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
- Gauss–Bonnet is topological in four dimensions: $f = R + \alpha\mathcal{G}$ gives plain GR
- $f(R) = R + \alpha R^2$ agrees between the covariant and minisuperspace routes, and matches
  the textbook $3FH^2 = \kappa\rho + \tfrac{1}{2}(FR - f) - 3H\dot F$
- $\Lambda$CDM: $E^2 = \Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_k(1+z)^2 + \Omega_\Lambda$,
  $q_0 = -1 + \tfrac{3}{2}\Omega_m$, acceleration from $z = 0.671$ at $\Omega_m = 0.3$

## License

MIT — see [LICENSE](LICENSE). Use it, change it, publish with it; just keep the copyright
notice. If it ends up being useful in published work, a citation is appreciated but not
required.
