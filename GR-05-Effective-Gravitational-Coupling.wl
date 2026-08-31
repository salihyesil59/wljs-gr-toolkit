(* ::Package:: *)

(* ============================================================================
   GR-05 : The Effective Gravitational Coupling, Derived
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-05-Effective-Gravitational-Coupling.wl

   ============================================================================ *)

(*::md::
# The Effective Gravitational Coupling

GR-04 solved the growth equation but took $G_{\text{eff}}/G$ from the literature. This
notebook closes that gap: it perturbs the field equations, applies the quasi-static
sub-horizon approximation, and derives $G_{\text{eff}}$ instead of quoting it.

Work in Newtonian (longitudinal) gauge,

$$
ds^2 = -\left(1+2\Psi\right)dt^2 + a(t)^2\left(1-2\Phi\right)\delta_{ij}dx^i dx^j ,
$$

and define the three functions that any metric theory has to supply:

$$
\frac{k^2}{a^2}\Psi = -4\pi G_{\text{eff}}\,\rho_m\delta_m , \qquad
\eta = \frac{\Phi}{\Psi} , \qquad
\Sigma = \frac{G_{\text{eff}}}{2G}\left(1+\eta\right) .
$$

$G_{\text{eff}}$ is what a galaxy falls with, so it drives growth. $\Sigma$ is what a photon
is deflected by, since lensing responds to the Weyl combination $\Phi+\Psi$. The slip $\eta$
measures the gap between them, and it is exactly $1$ in general relativity without anisotropic
stress. A theory that changes $G_{\text{eff}}$ but not $\Sigma$ shifts galaxy clustering while
leaving weak lensing alone, which is why the two are quoted separately.

## How the derivation runs

1. Write the metric with a bookkeeping parameter $\epsilon$ in front of $\Psi$ and $\Phi$, and
   put the perturbation in a single plane wave $e^{ikx}$.
2. Build the full non-linear field equations on that metric with $f$ left as an unspecified
   function, so the answer holds for every $f(R)$ at once.
3. Expand in $\epsilon$. Order $\epsilon^0$ returns the background Friedmann equations, which
   is a free check that the setup is right. Order $\epsilon^1$ is the perturbation system.
4. Apply the two approximations, stated explicitly below.
5. Solve the resulting algebra for $\Psi$ and $\Phi$ in terms of $\delta_m$, and read off the
   three functions.

## The two approximations, stated plainly

Everything above step 4 is exact. These two are not, and they are where the physics input sits:

- **Quasi-static.** Time derivatives of $\Psi$, $\Phi$ and $\delta$ are dropped. The
  perturbations are assumed to evolve on the Hubble timescale, so $\dot\Psi \sim H\Psi$, which
  is negligible next to the spatial gradients kept below.
- **Sub-horizon**, $k/a \gg H$. A metric perturbation survives only where it carries at least
  $k^2$ and no factor of $\dot a$ or $\ddot a$; anything with a Hubble factor in place of a
  gradient is smaller by $(aH/k)^2$. Matter source terms carry no $k$ at all and are kept.

Both are standard and both fail on scales approaching the horizon. The notebook implements the
second as an explicit term filter rather than burying it, so you can read the rule and change
it.

## Scope

This is metric perturbation theory, so it covers $f(R)$ and, with the field equations swapped,
any other theory whose fundamental variable is $g_{\mu\nu}$. It does **not** cover
$f(\mathcal{T})$ or $f(Q)$: their fundamental variables are a tetrad and a flat non-metric
connection, and perturbing those needs machinery this notebook does not have. Those entries in
GR-04's library remain quoted rather than derived. Section 10 says exactly what would have to
change.

## The series

| | |
|---|---|
| **GR-01** | a metric in, curvature tensors out |
| **GR-02** | an action in, field equations, Friedmann and $E(z)$ out |
| **GR-03** | $E(z)$ in, distances, times and the BAO ruler out |
| **GR-04** | $E(z)$ and $G_{\text{eff}}$ in, growth and $f\sigma_8$ out |
| **GR-05** | perturbed field equations in, $G_{\text{eff}}$, slip and lensing out |

This is the last link and the only one that feeds backwards: GR-04 treats $G_{\text{eff}}$ as
something you look up, and this is where the $f(R)$ entry it looks up is manufactured.
::*)

(*::md::
## 1. Setup

Curvature comes from the same Function Repository resources as the rest of the series.
Simplification is deliberately kept off the hot path: these expressions are large, and
`Simplify` is applied only after the perturbative order has been extracted.
::*)

(*::code::*)
ClearAll["Global`*"];

$grMetric = ResourceFunction["MetricTensor"];
$grRicci  = ResourceFunction["RicciTensor"];
$grChr    = ResourceFunction["ChristoffelSymbols"];

$coords = {t, x, y, z};
$wave   = Exp[I k x];        (* one plane wave along x is enough to expose the k^2 structure *)

(*::md::
## 2. The perturbed metric

$\epsilon$ is a bookkeeping parameter, not a physical quantity: it counts perturbative order
and is set to zero at the end of each `Series`. Putting the wave along a single axis costs no
generality for scalar perturbations and keeps every expression a function of $t$ and $x$ only,
which is what makes the whole derivation run in seconds instead of hours.
::*)

(*::code::*)
ClearAll[PerturbedMetric];

PerturbedMetric[] := DiagonalMatrix[{
   -(1 + 2 eps psi[t] $wave),
   a[t]^2 (1 - 2 eps phi[t] $wave),
   a[t]^2 (1 - 2 eps phi[t] $wave),
   a[t]^2 (1 - 2 eps phi[t] $wave)}];

gPert = PerturbedMetric[];
MatrixForm[gPert]

(*::md::
## 3. Curvature and the field equations

The metric is handed to the resources once and everything else is built from what comes back.
Leaving $f$ as the undefined symbol `F` is the trick that keeps the derivation general:
Wolfram differentiates it symbolically, so $f_R$, $f_{RR}$ and $f_{RRR}$ appear on their own as
`F'`, `F''` and `F'''` evaluated on the background scalar curvature.
::*)

(*::code::*)
ClearAll[covariantHessian, boxOperator];

mt  = $grMetric[gPert, $coords];
gdn = mt["MatrixRepresentation"];
gup = mt["InverseMetricTensor"]["MatrixRepresentation"];
chr = $grChr[mt]["TensorRepresentation"];
ric = $grRicci[mt]["MatrixRepresentation"];
rsc = Sum[gup[[i, j]] ric[[i, j]], {i, 4}, {j, 4}];

covariantHessian[f_] := Table[
   D[f, $coords[[m]], $coords[[q]]]
     - Sum[chr[[l, m, q]] D[f, $coords[[l]]], {l, 4}], {m, 4}, {q, 4}];

boxOperator[f_] := Sum[gup[[m, q]] covariantHessian[f][[m, q]], {m, 4}, {q, 4}];

(* pressureless matter, comoving to first order *)
uUp    = {1 - eps psi[t] $wave, 0, 0, 0};
uDn    = Table[Sum[gdn[[i, j]] uUp[[j]], {j, 4}], {i, 4}];
rhoTot = rho[t] (1 + eps del[t] $wave);
tmn    = Table[rhoTot uDn[[i]] uDn[[j]], {i, 4}, {j, 4}];

(* f(R) field equations, with f generic *)
fieldEq = With[{fv = F[rsc], fr = F'[rsc]},
  Table[fr ric[[i, j]] - (1/2) fv gdn[[i, j]] + gdn[[i, j]] boxOperator[fr]
      - covariantHessian[fr][[i, j]] - kap tmn[[i, j]], {i, 4}, {j, 4}]];

Dimensions[fieldEq]

(*::md::
## 4. Expanding in $\epsilon$

Order $\epsilon^0$ is the background. It should be exactly the $f(R)$ Friedmann constraint
GR-02 derives by minisuperspace variation, which is a genuine cross-check between two
independent routes through this series.
::*)

(*::code::*)
AbsoluteTiming[
 order0 = Table[Normal[Series[fieldEq[[i, j]], {eps, 0, 0}]], {i, 4}, {j, 4}];
 order1 = Table[Coefficient[Normal[Series[fieldEq[[i, j]], {eps, 0, 1}]], eps], {i, 4}, {j, 4}];
 backgroundR = Simplify[Normal[Series[rsc, {eps, 0, 0}]]];
]

(*::code::*)
Column[{
  Row[{Style["background R = ", Bold], backgroundR}],
  Row[{Style["background 00 equation:  ", Bold], Simplify[order0[[1, 1]]], " = 0"}]}]

(*::md::
## 5. The two approximations

`quasiStatic` kills time derivatives of the perturbation amplitudes. `subHorizon` is the term
filter: it keeps a term if the term has no metric perturbation in it at all (those are the
matter sources), or if it carries $k^2$ or more with no Hubble factor. Naming the background
$f$-derivatives first is what lets the filter see a $\dot a$ for what it is, rather than
tripping over the ones buried inside the argument of `F'`.
::*)

(*::code::*)
ClearAll[quasiStatic, subHorizon, reduceEquation];

abbreviations = {
  Derivative[3][F][backgroundR] -> fRRR,
  Derivative[2][F][backgroundR] -> fRR,
  Derivative[1][F][backgroundR] -> fR,
  F[backgroundR] -> f0};

quasiStatic = {Derivative[_][psi][t] -> 0, Derivative[_][phi][t] -> 0,
               Derivative[_][del][t] -> 0};

subHorizon[e_] := Total[Select[List @@ Expand[e],
   (FreeQ[#, psi | phi] ||
     (Exponent[#, k] >= 2 && FreeQ[#, Derivative[_][a][t]])) &]];

reduceEquation[e_] := subHorizon[Expand[Simplify[(e /. quasiStatic)/$wave] /. abbreviations]];

eq00       = reduceEquation[order1[[1, 1]]];
eqTraceless = reduceEquation[order1[[2, 2]] - order1[[3, 3]]];

Column[{
  Row[{Style["00 :          ", Bold], eq00, " = 0"}],
  Row[{Style["traceless ij: ", Bold], eqTraceless, " = 0"}]}]

(*::md::
## 6. The result

Two equations, two unknowns. Solving them gives $\Psi$ and $\Phi$ in terms of $\delta_m$, and
the three response functions follow by definition.
::*)

(*::code::*)
solution = Simplify[First @ Solve[{eq00 == 0, eqTraceless == 0}, {psi[t], phi[t]}]];

psiSol = Simplify[psi[t] /. solution];
phiSol = Simplify[phi[t] /. solution];

geffExpr  = Simplify[-2 k^2 psiSol/(a[t]^2 kap rho[t] del[t])];
slipExpr  = Simplify[phiSol/psiSol];
sigmaExpr = Simplify[geffExpr (1 + slipExpr)/2];

Column[{
  Row[{Style["G_eff / G = ", Bold], geffExpr}],
  Row[{Style["eta       = ", Bold], slipExpr}],
  Row[{Style["Sigma     = ", Bold], sigmaExpr}]}]

(*::md::
## 7. In the usual variables

Writing $m = \dfrac{k^2}{a^2}\dfrac{f_{RR}}{f_R}$, the same three expressions read

$$
\frac{G_{\text{eff}}}{G} = \frac{1}{f_R}\,\frac{1+4m}{1+3m}, \qquad
\eta = \frac{1+2m}{1+4m}, \qquad
\Sigma = \frac{1}{f_R}.
$$

The last one is worth pausing on, and it was not put in by hand: the scalaron cancels out of
the lensing combination entirely, so light in $f(R)$ gravity is deflected exactly as in general
relativity up to the constant rescaling $1/f_R$, with **no scale dependence at all**. Growth
feels $m$; lensing does not. That is a sharp, falsifiable signature, and it is why a
disagreement between clustering and lensing amplitudes is read as a hint of modified gravity.

Between the two limits $G_{\text{eff}}$ moves by a factor of $4/3$: on scales far outside the
scalaron's Compton wavelength ($m \to 0$) gravity is standard, and far inside it ($m \to \infty$)
it is a third stronger.
::*)

(*::code::*)
Module[{m, target, targetSlip},
 m = k^2 fRR/(a[t]^2 fR);
 target = (1/fR) (1 + 4 m)/(1 + 3 m);
 targetSlip = (1 + 2 m)/(1 + 4 m);
 Dataset @ {
  <|"statement" -> "G_eff/G = (1/fR)(1+4m)/(1+3m)",
    "ok" -> (Simplify[geffExpr - target] === 0)|>,
  <|"statement" -> "eta = (1+2m)/(1+4m)",
    "ok" -> (Simplify[slipExpr - targetSlip] === 0)|>,
  <|"statement" -> "Sigma = 1/fR, with no k dependence",
    "ok" -> (Simplify[sigmaExpr - 1/fR] === 0 && FreeQ[Simplify[sigmaExpr], k])|>,
  <|"statement" -> "general relativity (fR = 1, fRR = 0) gives G_eff = G and eta = 1",
    "ok" -> (Simplify[geffExpr /. {fR -> 1, fRR -> 0}] === 1 &&
             Simplify[slipExpr /. {fR -> 1, fRR -> 0}] === 1)|>,
  <|"statement" -> "large scales, k -> 0: G_eff -> G/fR and eta -> 1",
    "ok" -> (Simplify[Limit[geffExpr, k -> 0] - 1/fR] === 0 &&
             Simplify[Limit[slipExpr, k -> 0] - 1] === 0)|>,
  <|"statement" -> "small scales, k -> infinity: G_eff -> 4G/(3 fR) and eta -> 1/2",
    "ok" -> (Simplify[Limit[geffExpr, k -> Infinity] - 4/(3 fR)] === 0 &&
             Simplify[Limit[slipExpr, k -> Infinity] - 1/2] === 0)|>}]

(*::md::
## 8. A concrete model

Take Starobinsky's $f(R) = R + \alpha R^2$, so $f_R = 1 + 2\alpha R$ and $f_{RR} = 2\alpha$.
On a $\Lambda$CDM-like background the scalar curvature follows from $R = 6(\dot H + 2H^2)$,
which in redshift is

$$
R(z) = 6H_0^2\left[2E^2 - \tfrac{1}{2}(1+z)\frac{d(E^2)}{dz}\right] = 3H_0^2\left[\Omega_m(1+z)^3 + 4\Omega_\Lambda\right].
$$

Units are $H_0 = 1$ from here on, so $\tilde\alpha = \alpha H_0^2$ and $\tilde k = k/H_0$ are
both dimensionless. The plot shows the scale dependence that is the whole point: different
wavenumbers feel different gravity, which never happens in general relativity and never
happens in $f(\mathcal{T})$ or $f(Q)$ either.

**The coupling below is deliberately oversized.** $\tilde\alpha = 0.05$ gives
$f_R(0) \approx 1.92$, so $G_{\text{eff}}$ comes out roughly half of Newton's constant and
$f\sigma_8$ drops to about half its $\Lambda$CDM value — nothing like that survives contact
with data. It is chosen so the scale dependence is visible on a linear axis rather than being
a line thickness. Turn $\tilde\alpha$ down to $10^{-3}$ or below for a model you would
actually defend, and the curves converge onto each other, which is itself the point: viable
$f(R)$ has to hide.
::*)

(*::code::*)
ClearAll[om, ol, alpha, Rbg, fRof, fRRof, geffOf];

om = 0.315; ol = 1 - om;
Rbg[zz_]   := 3 (om (1 + zz)^3 + 4 ol);          (* H0 = 1 *)
alpha      = 0.05;
fRof[zz_]  := 1 + 2 alpha Rbg[zz];
fRRof      = 2 alpha;
geffOf[zz_, kv_] := Module[{m},
   m = kv^2 (1 + zz)^2 fRRof/fRof[zz];           (* k^2/a^2 with a = 1/(1+z) *)
   (1/fRof[zz]) (1 + 4 m)/(1 + 3 m)];

Column[{
  Row[{"R(0) = ", Rbg[0], " H0^2,     fR(0) = ", fRof[0]}],
  Plot[Evaluate[Table[geffOf[z, kv], {kv, {0.3, 1, 3, 10}}]], {z, 0, 3},
   PlotLegends -> {"k = 0.3", "k = 1", "k = 3", "k = 10"},
   AxesLabel -> {"z", "G_eff / G"}, PlotRange -> All, ImageSize -> 440,
   PlotLabel -> "f(R) = R + 0.05 R^2, in units H0 = 1"]}]

(*::md::
## 9. Feeding it back into growth

This is the loop closing. GR-04's growth equation is reused verbatim, but $G_{\text{eff}}$ is
now the derived expression rather than a library entry. Because it depends on $k$, so does the
growth: small scales grow faster. In general relativity every curve here would lie on top of
the others.
::*)

(*::code::*)
Module[{eBg, growthFor, s8 = 0.811, zi = 100, curves},
 eBg = Sqrt[om (1 + z)^3 + ol];

 growthFor[geffFun_] := Module[{del2, sol},
   sol = NDSolve[{
      del2''[z] + (D[eBg, z]/eBg - 1/(1 + z)) del2'[z]
        - (3/2) om (1 + z) geffFun[z] del2[z]/eBg^2 == 0,
      del2[zi] == 1/(1 + zi), del2'[zi] == -1/(1 + zi)^2},
     del2, {z, 0, zi}, MaxSteps -> 10^6];
   del2 /. First[sol]];

 curves = Table[
   With[{kv = kv}, {kv, growthFor[Function[zz, geffOf[zz, kv]]]}],
   {kv, {0.3, 1, 3, 10}}];

 Column[{
   Grid[Prepend[
     {#[[1]], -s8 (1 + 0.) #[[2]]'[0]/#[[2]][0]} & /@ curves,
     {Style["k / H0", Bold], Style["f sigma8 (0)", Bold]}],
    Alignment -> Left, Spacings -> {3, 0.7}],
   Plot[Evaluate[(-s8 (1 + z) #[[2]]'[z]/#[[2]][0]) & /@ curves], {z, 0, 3},
    PlotLegends -> (("k = " <> ToString[#[[1]]]) & /@ curves),
    AxesLabel -> {"z", "f\[Sigma]8"}, PlotRange -> All, ImageSize -> 440,
    PlotLabel -> "scale-dependent growth, the f(R) fingerprint"]}]]

(*::md::
## 10. What would have to change for $f(\mathcal{T})$ or $f(Q)$

Nothing in sections 2 to 6 assumed anything about $f(R)$ except the form of the field
equations fed into `fieldEq`. Swap those for another **metric** theory and the rest of the
pipeline runs unchanged.

Torsion and non-metricity break that, and not on a technicality. In both the object being
varied is something other than the metric.

**Teleparallel, $f(\mathcal{T})$.** The variable is the tetrad $e^A{}_\mu$: sixteen components
against the metric's ten. The extra six are exactly the dimension of the local Lorentz group,
and in TEGR they are pure gauge, which is why the linear theory can be general relativity in
disguise at all. For a non-linear $f$ that invariance is lost in the pure-tetrad formulation,
so the six stop being gauge and start being physics &mdash; two different tetrads building the
*same* metric then satisfy different field equations. This is the good-tetrad problem, and it
means a perturbation calculation has to commit to a tetrad, not merely to a metric. The
covariant formulation buys the invariance back by carrying a spin connection alongside, at the
price of that connection having equations of its own.

**Symmetric teleparallel, $f(Q)$.** The variable is a flat, torsion-free connection, and the
coincident gauge that makes $Q$ come out as the tidy $-6H^2$ used in GR-02 is a statement about
that connection rather than about coordinates in the usual sense. It does not stack for free
on top of Newtonian gauge: fixing the connection spends freedom that the metric gauge choice
also wants. A consistent calculation either keeps the connection perturbations explicitly or
gives up the convenient form of $Q$.

Both families also carry known perturbative pathologies. Extra modes propagate, and their
kinetic terms can degenerate on cosmological backgrounds, so linear theory is not automatically
the right description of them &mdash; the strong-coupling problem. A correct linear calculation
is necessary but not obviously sufficient.

None of this is unreachable. It is a different derivation over different variables, and it is
the natural starting point for a GR-06. What matters here is that the boundary is stated: the
$f(\mathcal{T})$ and $f(Q)$ entries in GR-04's library are quoted from the literature, and
nothing in this notebook licenses them.
::*)

(*::code::*)
Dataset @ <|
  "derived here"      -> "GR, and any metric theory: f(R) shown in full",
  "quoted in GR-04"   -> "f(T), f(Q): need tetrad and connection perturbations",
  "approximations"    -> "quasi-static and sub-horizon, both fail near the horizon",
  "cross-check"       -> "order eps^0 reproduces the GR-02 Friedmann constraint"|>
