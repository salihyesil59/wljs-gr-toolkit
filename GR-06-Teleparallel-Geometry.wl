(* ::Package:: *)

(* ============================================================================
   GR-06 : Teleparallel Geometry and f(T) Cosmology
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-06-Teleparallel-Geometry.wl

   ============================================================================ *)

(*::md::
# Teleparallel Geometry and $f(\mathcal{T})$ Cosmology

GR-05 derived the effective gravitational coupling for $f(R)$ and stopped at the boundary of
metric theories, saying that $f(\mathcal{T})$ needs a tetrad and that perturbing $g_{\mu\nu}$
alone cannot reach it. This notebook crosses into that territory and builds the geometry:
tetrad, torsion, contortion, superpotential, torsion scalar, and from those the $f(\mathcal{T})$
Friedmann equations, derived on the tetrad side and checked against the metric-side answer
GR-02 got by a completely different route.

It stops short of $G_{\text{eff}}$, and section 8 shows exactly where and why. That section is
not an apology; it is the most informative part of the notebook.

## Teleparallel gravity in one paragraph

General relativity puts gravity in the curvature of the Levi-Civita connection, which is
metric-compatible and torsion-free. Teleparallel gravity makes the opposite choice: a
connection that is flat and metric-compatible but has **torsion**. Curvature is identically
zero, and what encodes gravity is the failure of parallel-transported frames to close. The
dynamical object is the tetrad $e^A{}_\mu$, a set of four vector fields with
$g_{\mu\nu} = \eta_{AB}\,e^A{}_\mu e^B{}_\nu$, and the Weitzenböck connection built from it,
$\Gamma^\lambda{}_{\nu\mu} = e_A{}^\lambda\,\partial_\mu e^A{}_\nu$.

$$
T^\lambda{}_{\mu\nu} = e_A{}^\lambda\left(\partial_\mu e^A{}_\nu - \partial_\nu e^A{}_\mu\right),
\qquad
\mathcal{T} = S_\lambda{}^{\mu\nu}\,T^\lambda{}_{\mu\nu}
$$

The remarkable fact is that $\mathcal{T}$ and the Levi-Civita Ricci scalar differ only by a
total derivative,

$$
R = -\mathcal{T} + B, \qquad B = \frac{2}{e}\,\partial_\mu\!\left(e\,T^{\nu\mu}{}_{\nu}\right),
$$

so the two Lagrangians give the same field equations. That is TEGR, the teleparallel
equivalent of general relativity. The equivalence breaks the moment the Lagrangian stops being
linear: $f(R)$ and $f(\mathcal{T})$ are genuinely different theories, because a boundary term
inside a non-linear function is no longer a boundary term.

The notebook verifies that identity rather than assuming it, on three separate geometries.

## Conventions, and the sign that has to be watched

Here $\mathcal{T}$ is the standard superpotential contraction, which gives $\mathcal{T} = +6H^2$
in flat FLRW. GR-02 uses the opposite sign, `Ts` $= -6H^2$, chosen there so that a *linear*
$f$ reproduces general relativity. The two are related by $\texttt{Ts} = -\mathcal{T}$, so TEGR
is $F(\mathcal{T}) = -\mathcal{T}$ here and $f(\texttt{Ts}) = \texttt{Ts}$ there. Section 7
checks that the two conventions give the same Friedmann equations, which is the only thing
that matters.
::*)

(*::md::
## The series

| | |
|---|---|
| **GR-01** | a metric in, curvature tensors out |
| **GR-02** | an action in, field equations, Friedmann and $E(z)$ out |
| **GR-03** | $E(z)$ in, distances, times and the BAO ruler out |
| **GR-04** | $E(z)$ and $G_{\text{eff}}$ in, growth and $f\sigma_8$ out |
| **GR-05** | perturbed field equations in, $G_{\text{eff}}$, slip and lensing out |
| **GR-06** | a tetrad in, torsion, the TEGR identity and $f(\mathcal{T})$ cosmology out |

This one branches off the main line. Everything before it lives on the metric; this is the
first notebook whose fundamental variable is something else.
::*)

(*::md::
## 1. Setup

Only the curvature backend is borrowed, and only so that the teleparallel results can be
checked against Levi-Civita ones computed independently.
::*)

(*::code::*)
ClearAll["Global`*"];

$grMetric = ResourceFunction["MetricTensor"];
$grRicci  = ResourceFunction["RicciTensor"];

$eta = DiagonalMatrix[{-1, 1, 1, 1}];

(*::md::
## 2. The teleparallel toolkit

Every object below is built from its definition, with no formula quoted from a paper:

$$
K^{\mu\nu}{}_{\lambda} = -\tfrac{1}{2}\left(T^{\mu\nu}{}_{\lambda} - T^{\nu\mu}{}_{\lambda} - T_{\lambda}{}^{\mu\nu}\right),
\qquad
S_{\lambda}{}^{\mu\nu} = \tfrac{1}{2}\left(K^{\mu\nu}{}_{\lambda} + \delta^{\mu}_{\lambda}T^{\alpha\nu}{}_{\alpha} - \delta^{\nu}_{\lambda}T^{\alpha\mu}{}_{\alpha}\right)
$$

Index placement in these two is where the whole calculation can quietly go wrong, so the
function returns the torsion scalar computed **twice**, once as the superpotential contraction
$S_\lambda{}^{\mu\nu}T^\lambda{}_{\mu\nu}$ and once through the equivalent three-term form

$$
\mathcal{T} = \tfrac{1}{4}T^{\lambda\mu\nu}T_{\lambda\mu\nu} + \tfrac{1}{2}T^{\lambda\mu\nu}T_{\nu\mu\lambda} - v^\mu v_\mu,
\qquad v_\mu = T^{\lambda}{}_{\lambda\mu}.
$$

The two share no intermediate quantity beyond the torsion tensor itself, so agreement between
them is a real check on the index bookkeeping. Getting this wrong the first time produced
$\tfrac{21}{2}H^2$ instead of $6H^2$, which is exactly the sort of error that a single
calculation cannot catch.
::*)

(*::code::*)
ClearAll[Teleparallel, RicciScalarOf, BoundaryTerm];

Teleparallel[tet_?MatrixQ, xs_List] :=
 Module[{n = Length[xs], iv, eup, g, ig, det, tor, tDn, tUp12, tUp23, tUpAll,
         kon, tr1, sup, tS, vec, vup, tS2, capA, capB, l, m, i, p, q, al},

  iv  = Inverse[tet];
  eup = Table[iv[[m, capA]], {capA, n}, {m, n}];                 (* e_A^mu *)
  g   = Table[Sum[$eta[[capA, capB]] tet[[capA, i]] tet[[capB, m]],
     {capA, n}, {capB, n}], {i, n}, {m, n}];
  ig  = Inverse[g];
  det = Det[tet];

  (* T^l_{m i} *)
  tor = Table[
    Sum[eup[[capA, l]] (D[tet[[capA, i]], xs[[m]]] - D[tet[[capA, m]], xs[[i]]]), {capA, n}],
    {l, n}, {m, n}, {i, n}];

  tDn    = Table[Sum[g[[l, q]] tor[[q, m, i]], {q, n}], {l, n}, {m, n}, {i, n}];
  tUp12  = Table[Sum[ig[[m, p]] ig[[i, q]] tDn[[p, q, l]], {p, n}, {q, n}],
     {m, n}, {i, n}, {l, n}];
  tUp23  = Table[Sum[ig[[m, p]] ig[[i, q]] tDn[[l, p, q]], {p, n}, {q, n}],
     {l, n}, {m, n}, {i, n}];
  tUpAll = Table[Sum[ig[[l, al]] tUp23[[al, m, i]], {al, n}], {l, n}, {m, n}, {i, n}];

  kon = Table[-(1/2) (tUp12[[m, i, l]] - tUp12[[i, m, l]] - tUp23[[l, m, i]]),
    {m, n}, {i, n}, {l, n}];
  tr1 = Table[Sum[tUp12[[al, i, al]], {al, n}], {i, n}];
  sup = Table[(1/2) (kon[[m, i, l]]
      + KroneckerDelta[m, l] tr1[[i]] - KroneckerDelta[i, l] tr1[[m]]),
    {l, n}, {m, n}, {i, n}];

  tS = Simplify[Sum[sup[[l, m, i]] tor[[l, m, i]], {l, n}, {m, n}, {i, n}]];

  vec = Table[Sum[tor[[l, l, m]], {l, n}], {m, n}];
  vup = Table[Sum[ig[[m, p]] vec[[p]], {p, n}], {m, n}];
  tS2 = Simplify[
    (1/4) Sum[tUpAll[[l, m, i]] tDn[[l, m, i]], {l, n}, {m, n}, {i, n}]
     + (1/2) Sum[tUpAll[[l, m, i]] tDn[[i, m, l]], {l, n}, {m, n}, {i, n}]
     - Sum[vup[[m]] vec[[m]], {m, n}]];

  <|"Tetrad" -> tet, "TetradInverse" -> eup, "Metric" -> Simplify[g], "MetricInverse" -> ig,
    "Determinant" -> det, "Torsion" -> tor, "Superpotential" -> sup,
    "TorsionVector" -> vec, "TorsionScalar" -> tS, "TorsionScalarAlt" -> tS2,
    "Coordinates" -> xs|>];

(* independent Levi-Civita comparison, same backend as GR-01 *)
RicciScalarOf[g_?MatrixQ, coords_List] :=
 Module[{mt = $grMetric[g, coords]},
  Simplify[Tr[mt["InverseMetricTensor"]["MatrixRepresentation"] .
     $grRicci[mt]["MatrixRepresentation"]]]];

(* B = (2/e) d_mu (e v^mu) *)
BoundaryTerm[tp_Association] :=
 Module[{xs = tp["Coordinates"], n, e = tp["Determinant"], vup, m, p},
  n = Length[xs];
  vup = Table[Sum[tp["MetricInverse"][[m, p]] tp["TorsionVector"][[p]], {p, n}], {m, n}];
  Simplify[(2/e) Sum[D[e vup[[m]], xs[[m]]], {m, n}]]];

(*::md::
## 3. Flat FLRW

The diagonal tetrad $e^A{}_\mu = \mathrm{diag}(1, a, a, a)$ reproduces the FLRW metric, and
its torsion scalar is $\mathcal{T} = 6H^2$ &mdash; built purely from first derivatives, where
the Ricci scalar needs second ones. That difference is the whole of teleparallel gravity in one
observation.
::*)

(*::code::*)
xs = {t, x, y, z};

flrw = Teleparallel[DiagonalMatrix[{1, a[t], a[t], a[t]}], xs];

Column[{
  Row[{"metric from the tetrad : ", MatrixForm[flrw["Metric"]]}],
  Row[{"torsion scalar         : ", flrw["TorsionScalar"]}],
  Row[{"same, three-term form  : ", flrw["TorsionScalarAlt"]}],
  Row[{"Levi-Civita R          : ", RicciScalarOf[flrw["Metric"], xs]}],
  Row[{"boundary term B        : ", BoundaryTerm[flrw]}]}]

(*::md::
## 4. The TEGR identity, checked three ways

$R = -\mathcal{T} + B$ is what makes teleparallel gravity equivalent to general relativity at
the linear level, and it is an identity, so it has to hold for **every** tetrad, not just
convenient ones. The three below are deliberately different in character: a time-dependent
cosmology, a static vacuum with a horizon, and an exponentially expanding one.
::*)

(*::code::*)
Module[{ys = {t, r, th, ph}, cases},
 cases = {
   {"flat FLRW", Teleparallel[DiagonalMatrix[{1, a[t], a[t], a[t]}], xs], xs},
   {"Schwarzschild tetrad",
    Teleparallel[DiagonalMatrix[{Sqrt[1 - 2 M/r], 1/Sqrt[1 - 2 M/r], r, r Sin[th]}], ys], ys},
   {"de Sitter", Teleparallel[DiagonalMatrix[{1, Exp[hh t], Exp[hh t], Exp[hh t]}], xs], xs}};

 Dataset[
  Function[c,
   Module[{tp = c[[2]], cs = c[[3]], rr},
    rr = RicciScalarOf[tp["Metric"], cs];
    <|"tetrad" -> c[[1]],
      "T" -> tp["TorsionScalar"],
      "R" -> rr,
      "both T forms agree" -> (Simplify[tp["TorsionScalar"] - tp["TorsionScalarAlt"]] === 0),
      "R = -T + B" -> (Simplify[rr - (-tp["TorsionScalar"] + BoundaryTerm[tp])] === 0)|>]] /@ cases]]

(*::md::
## 5. The good-tetrad problem, computed

GR-05 asserted that in $f(\mathcal{T})$ the extra tetrad components stop being gauge. Here is
that statement as arithmetic rather than prose.

Both tetrads below build the **same** metric &mdash; flat Minkowski space in spherical
coordinates &mdash; and they are related by a local Lorentz rotation. In general relativity
nothing could possibly distinguish them. They have different torsion scalars.

For TEGR that difference is harmless: it sits entirely in the boundary term $B$, and the field
equations never see it. Put the same two tetrads inside a non-linear $f$ and they become
different theories of the same spacetime. Choosing the tetrad is now physics, and the
literature's "good tetrad" is the one whose field equations do not force $f_{\mathcal{TT}} = 0$
and collapse the theory back to TEGR.
::*)

(*::code::*)
Module[{ys = {t, r, th, ph}, diag, rot, tpD, tpR},
 diag = DiagonalMatrix[{1, 1, r, r Sin[th]}];
 rot  = {{1, 0, 0, 0},
         {0, Sin[th] Cos[ph], r Cos[th] Cos[ph], -r Sin[th] Sin[ph]},
         {0, Sin[th] Sin[ph], r Cos[th] Sin[ph],  r Sin[th] Cos[ph]},
         {0, Cos[th],        -r Sin[th],          0}};
 tpD = Teleparallel[diag, ys];
 tpR = Teleparallel[rot, ys];

 Dataset @ <|
   "the two tetrads build the same metric" ->
     (Simplify[tpD["Metric"] - tpR["Metric"]] === ConstantArray[0, {4, 4}]),
   "that metric is flat, R =" -> RicciScalarOf[tpD["Metric"], ys],
   "T from the diagonal tetrad" -> tpD["TorsionScalar"],
   "T from the rotated tetrad"  -> tpR["TorsionScalar"],
   "so T is not a function of the metric alone" ->
     (Simplify[tpD["TorsionScalar"] - tpR["TorsionScalar"]] =!= 0)|>]

(*::md::
## 6. Field equations by direct variation

Rather than quote the $f(\mathcal{T})$ field equations, vary the action. The Lagrangian density
$e\,f(\mathcal{T})$ contains the tetrad and its first derivatives only, so the Euler-Lagrange
operator is the first-order one, and for a homogeneous ansatz there is a single independent
variable.

For FLRW this is legitimate minisuperspace variation, exactly as in GR-02: the lapse and the
scale factor are varied, and the lapse is set to one only afterwards.
::*)

(*::code::*)
ClearAll[eulerT];

eulerT[lag_, q_Symbol] := D[lag, q[t]] - D[D[lag, q'[t]], t];

bg = Teleparallel[DiagonalMatrix[{Nl[t], a[t], a[t], a[t]}], xs];

Column[{
  Row[{"T with a lapse : ", bg["TorsionScalar"]}],
  Row[{"e = det(tetrad): ", bg["Determinant"]}]}]

(*::md::
## 7. The $f(\mathcal{T})$ Friedmann equations

Matter enters as $-e\,\rho(a)$ and $\rho'(a)$ is eliminated afterwards with the continuity
equation, the same treatment GR-02 uses. What comes out is

$$
\kappa\rho = \tfrac{1}{2}\left(f - 2\mathcal{T}f_{\mathcal{T}}\right),
$$

and setting $f = -\mathcal{T}$, which is TEGR in this sign convention, returns
$3H^2 = \kappa\rho$ and the standard acceleration equation. That is the cross-check worth
having: GR-02 reached the same equations from the metric side by varying a completely
different Lagrangian, and the two agree.
::*)

(*::code::*)
Module[{lag, toN1, contin, names, fr1, fr2, grSub, tv},
 lag = bg["Determinant"] (F[bg["TorsionScalar"]]/(2 \[Kappa]) - dens[a[t]]);

 toN1   = {Derivative[_][Nl][t] -> 0, Nl[t] -> 1};
 contin = {Derivative[1][dens][a[t]] -> -3 (dens[a[t]] + pres[a[t]])/a[t]};
 names  = {dens[a[t]] -> \[Rho], pres[a[t]] -> p};

 fr1 = Simplify[(eulerT[lag, Nl] /. toN1 /. contin /. names)/a[t]^3];
 fr2 = Simplify[(eulerT[lag, a]  /. toN1 /. contin /. names)/a[t]^2];
 grSub = {F -> Function[u, -u]};
 tv = 6 a'[t]^2/a[t]^2;

 Column[{
   Row[{Style["Friedmann I  : ", Bold], fr1, " = 0"}],
   Row[{Style["Friedmann II : ", Bold], Short[fr2, 3], " = 0"}],
   "",
   Dataset @ {
    <|"check" -> "Friedmann I is kappa rho = (f - 2 T f_T)/2",
      "ok" -> (Simplify[fr1 - ((F[tv] - 2 tv F'[tv])/(2 \[Kappa]) - \[Rho])] === 0)|>,
    <|"check" -> "TEGR, f = -T, gives 3H^2 = kappa rho",
      "ok" -> (Simplify[(fr1 /. grSub) - (3 a'[t]^2/(\[Kappa] a[t]^2) - \[Rho])] === 0)|>,
    <|"check" -> "TEGR gives the standard acceleration equation",
      "ok" -> (Simplify[((fr2 /. grSub) /. First[Solve[(fr1 /. grSub) == 0, \[Rho]]])
          - (3 (a'[t]^2 + 2 a[t] a''[t])/(\[Kappa] a[t]^2) + 3 p)] === 0)|>,
    <|"check" -> "agrees with GR-02, whose Ts = -T and whose f(Ts) = Ts is TEGR",
      "ok" -> (Simplify[((F[tv] - 2 tv F'[tv])/2 /. grSub) - tv/2] === 0)|>}}]]

(*::md::
## 8. Where this stops, and why

The natural next step is $G_{\text{eff}}$ for $f(\mathcal{T})$, by the route GR-05 used for
$f(R)$: perturb, expand, truncate, solve. Running that route here fails, and the failure is
worth more than a result would have been.

Take the perturbed diagonal tetrad in Newtonian gauge,
$e^A{}_\mu = \mathrm{diag}\!\left(1+\Psi,\; a(1-\Phi),\; a(1-\Phi),\; a(1-\Phi)\right)$,
which builds exactly the perturbed metric GR-05 used. Two things then happen.

First, $\delta\mathcal{T}$ carries **no spatial gradient at all**. The linear perturbation of
the torsion scalar comes out proportional to $H(\Psi H + \dot\Phi)$, with no $k^2$ anywhere.
That is not an artefact: $\mathcal{T}$ is quadratic in torsion, the background torsion is
purely of the $T^i{}_{0j}$ type, and the gradient-carrying perturbation components have index
structures that do not contract with it.

Second, and decisively, varying this tetrad family and truncating the way GR-05 does returns
**zero** &mdash; and it does so even in the TEGR limit, where the answer has to be the Poisson
equation $k^2\Phi/a^2 \propto \rho\delta$ that GR-05 derives without difficulty. A method that
loses general relativity cannot be trusted with anything else.

The diagnosis is the point. The diagonal tetrad has six fewer components than a general one,
and those six are precisely the local Lorentz directions that GR-05's scope note flagged. In
TEGR they are pure gauge and nothing is lost by fixing them; in a perturbative calculation
they carry constraint equations, and dropping them drops the gradient content. A trustworthy
$f(\mathcal{T})$ perturbation calculation has to vary all sixteen tetrad components, or work
from the second-order action with the full scalar sector present, and neither is a small
extension of what is here.

So the honest state of the series is:

- $f(R)$: $G_{\text{eff}}$ **derived**, in GR-05.
- $f(\mathcal{T})$: geometry, TEGR equivalence and background cosmology **derived**, here.
  $G_{\text{eff}}$ still **quoted** in GR-04's library.
- $f(Q)$: nothing derived yet.

The cell below reproduces the diagnosis rather than hiding it.
::*)

(*::code::*)
Module[{tet, tp, ans, dT, pw},
 tet = DiagonalMatrix[{aa[t, x], bb[t, x], cc[t, x], dd[t, x]}];
 tp  = Teleparallel[tet, xs];
 ans = {aa -> Function[{t, x}, 1 + eps psi[t] Exp[I kk x]],
        bb -> Function[{t, x}, a[t] (1 - eps phi[t] Exp[I kk x])],
        cc -> Function[{t, x}, a[t] (1 - eps phi[t] Exp[I kk x])],
        dd -> Function[{t, x}, a[t] (1 - eps phi[t] Exp[I kk x])]};
 dT = Simplify[
    Coefficient[Normal[Series[tp["TorsionScalar"] /. ans, {eps, 0, 1}]], eps]/Exp[I kk x]];
 pw = Exponent[Expand[dT], kk];

 Column[{
   Row[{Style["delta T at first order : ", Bold], dT}],
   Row[{Style["power of k in it       : ", Bold], pw}],
   Dataset @ {
     <|"observation" -> "delta T carries no spatial gradient",
       "ok" -> (pw <= 0)|>,
     <|"observation" -> "so the diagonal tetrad cannot produce a k^2 Poisson term",
       "ok" -> (pw < 2)|>}}]]

(*::md::
## 9. Checks

Everything the notebook claims, gathered in one place. The torsion scalar is computed twice by
routes sharing nothing but the torsion tensor; the TEGR identity is tested on three unrelated
geometries; the background cosmology is matched against GR-02, which reached it from the metric
side; and the perturbative limitation of section 8 is recorded as a check too, so that if a
later version fixes it the check will start failing and say so.
::*)

(*::code::*)
Module[{ys = {t, r, th, ph}, f1, s1, d1, tp2, tv, lag, fr1, grSub, toN1, contin, names},
 f1 = Teleparallel[DiagonalMatrix[{1, a[t], a[t], a[t]}], xs];
 s1 = Teleparallel[DiagonalMatrix[{Sqrt[1 - 2 M/r], 1/Sqrt[1 - 2 M/r], r, r Sin[th]}], ys];
 d1 = Teleparallel[DiagonalMatrix[{1, Exp[hh t], Exp[hh t], Exp[hh t]}], xs];

 lag = bg["Determinant"] (F[bg["TorsionScalar"]]/(2 \[Kappa]) - dens[a[t]]);
 toN1 = {Derivative[_][Nl][t] -> 0, Nl[t] -> 1};
 contin = {Derivative[1][dens][a[t]] -> -3 (dens[a[t]] + pres[a[t]])/a[t]};
 names = {dens[a[t]] -> \[Rho], pres[a[t]] -> p};
 fr1 = Simplify[(eulerT[lag, Nl] /. toN1 /. contin /. names)/a[t]^3];
 grSub = {F -> Function[u, -u]};
 tv = 6 a'[t]^2/a[t]^2;

 Dataset @ {
  <|"check" -> "tetrad reproduces the FLRW metric",
    "ok" -> (Simplify[f1["Metric"] - DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}]]
       === ConstantArray[0, {4, 4}])|>,
  <|"check" -> "T = 6 H^2 in flat FLRW",
    "ok" -> (Simplify[f1["TorsionScalar"] - 6 a'[t]^2/a[t]^2] === 0)|>,
  <|"check" -> "superpotential and three-term forms agree, all three tetrads",
    "ok" -> AllTrue[{f1, s1, d1},
       Simplify[#["TorsionScalar"] - #["TorsionScalarAlt"]] === 0 &]|>,
  <|"check" -> "R = -T + B for flat FLRW",
    "ok" -> (Simplify[RicciScalarOf[f1["Metric"], xs]
        - (-f1["TorsionScalar"] + BoundaryTerm[f1])] === 0)|>,
  <|"check" -> "R = -T + B for the Schwarzschild tetrad",
    "ok" -> (Simplify[RicciScalarOf[s1["Metric"], ys]
        - (-s1["TorsionScalar"] + BoundaryTerm[s1])] === 0)|>,
  <|"check" -> "R = -T + B for de Sitter",
    "ok" -> (Simplify[RicciScalarOf[d1["Metric"], xs]
        - (-d1["TorsionScalar"] + BoundaryTerm[d1])] === 0)|>,
  <|"check" -> "Friedmann I is kappa rho = (f - 2 T f_T)/2",
    "ok" -> (Simplify[fr1 - ((F[tv] - 2 tv F'[tv])/(2 \[Kappa]) - \[Rho])] === 0)|>,
  <|"check" -> "TEGR limit returns 3H^2 = kappa rho, matching GR-02",
    "ok" -> (Simplify[(fr1 /. grSub) - (3 a'[t]^2/(\[Kappa] a[t]^2) - \[Rho])] === 0)|>,
  <|"check" -> "the good-tetrad problem is real: same metric, different T",
    "ok" -> Module[{p1, p2},
      p1 = Teleparallel[DiagonalMatrix[{1, 1, r, r Sin[th]}], ys];
      p2 = Teleparallel[{{1, 0, 0, 0},
         {0, Sin[th] Cos[ph], r Cos[th] Cos[ph], -r Sin[th] Sin[ph]},
         {0, Sin[th] Sin[ph], r Cos[th] Sin[ph], r Sin[th] Cos[ph]},
         {0, Cos[th], -r Sin[th], 0}}, ys];
      Simplify[p1["Metric"] - p2["Metric"]] === ConstantArray[0, {4, 4}] &&
       Simplify[p1["TorsionScalar"] - p2["TorsionScalar"]] =!= 0]|>,
  <|"check" -> "section 8 still stands: no gradient in delta T for the diagonal tetrad",
    "ok" -> Module[{tp, ans, dT},
      tp = Teleparallel[DiagonalMatrix[{aa[t, x], bb[t, x], cc[t, x], dd[t, x]}], xs];
      ans = {aa -> Function[{t, x}, 1 + eps psi[t] Exp[I kk x]],
             bb -> Function[{t, x}, a[t] (1 - eps phi[t] Exp[I kk x])],
             cc -> Function[{t, x}, a[t] (1 - eps phi[t] Exp[I kk x])],
             dd -> Function[{t, x}, a[t] (1 - eps phi[t] Exp[I kk x])]};
      dT = Simplify[Coefficient[
         Normal[Series[tp["TorsionScalar"] /. ans, {eps, 0, 1}]], eps]/Exp[I kk x]];
      Exponent[Expand[dT], kk] <= 0]|>}]
