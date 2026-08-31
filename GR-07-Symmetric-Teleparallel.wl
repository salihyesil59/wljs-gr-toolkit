(* ::Package:: *)

(* ============================================================================
   GR-07 : Symmetric Teleparallel Geometry and f(Q)
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-07-Symmetric-Teleparallel.wl

   ============================================================================ *)

(*::md::
# Symmetric Teleparallel Geometry and $f(Q)$

GR-06 crossed from curvature into torsion. This one crosses into **non-metricity**, the third
corner of the geometric trinity, and closes the last gap in the series: the $f(Q)$ entry in
GR-04's library, the only expression left that was still quoted rather than derived.

The answer, derived in section 12, is

$$
\frac{G_{\text{eff}}}{G} = \frac{1}{f_Q}, \qquad \eta = \frac{\Phi}{\Psi} = 1,
$$

scale independent and with no gravitational slip &mdash; structurally like $f(\mathcal{T})$ and
unlike $f(R)$.

## The third geometry

General relativity uses a connection that is metric-compatible and torsion-free, and puts
gravity in its curvature. Teleparallel gravity keeps metric compatibility, drops curvature, and
uses torsion. Symmetric teleparallel gravity drops **both** curvature and torsion, and what is
left is the failure of the connection to preserve the metric:

$$
Q_{\lambda\mu\nu} = \nabla_\lambda g_{\mu\nu} = \partial_\lambda g_{\mu\nu} - \Gamma^\rho{}_{\lambda\mu}g_{\rho\nu} - \Gamma^\rho{}_{\lambda\nu}g_{\mu\rho}
$$

The non-metricity scalar is the particular quadratic combination

$$
Q = -\tfrac{1}{4}Q_{\alpha\beta\gamma}Q^{\alpha\beta\gamma} + \tfrac{1}{2}Q_{\alpha\beta\gamma}Q^{\gamma\beta\alpha} + \tfrac{1}{4}Q_\alpha Q^\alpha - \tfrac{1}{2}Q_\alpha\tilde{Q}^\alpha,
$$

with $Q_\alpha = Q_\alpha{}^\mu{}_\mu$ and $\tilde Q_\alpha = Q^\mu{}_{\alpha\mu}$, and it obeys
the same kind of identity as the torsion scalar does: $R$ and $Q$ differ by a total derivative,
so the linear theory &mdash; STEGR &mdash; is general relativity again, and the equivalence
dies as soon as $f$ is non-linear.

## Two fields, not one

The awkward part of $f(Q)$, and the reason GR-05 and GR-06 kept deferring it, is that the
connection is a **second field**. It is required to be flat and torsion-free, which is exactly
the statement that it can be written from four functions $\xi^\rho$,

$$
\Gamma^\lambda{}_{\mu\nu} = \frac{\partial x^\lambda}{\partial \xi^\rho}\,\partial_\mu\partial_\nu \xi^\rho ,
$$

so a general flat connection is four Stückelberg fields, no more and no less. Setting
$\xi^\rho = x^\rho$ makes $\Gamma$ vanish; that is the **coincident gauge**, and it is what makes
$Q$ come out as the tidy $-6H^2$ used in GR-02.

Here is the tension this notebook has to resolve. A diffeomorphism shifts $\delta\xi^\mu$ by the
diffeomorphism parameter, so the coincident gauge *is* a gauge choice, and it uses up the same
freedom that Newtonian gauge on the metric wants. You cannot have both. This notebook spends
the freedom on Newtonian gauge, keeps the connection perturbations $\delta\xi^\mu$ as
**physical fields**, and varies with respect to them along with the metric. Their equations then
close the system.

## Conventions

The quadratic form above gives $Q = -6H^2$ in flat FLRW, which is already the sign GR-02 and
GR-04 use for `Qs`: a linear $f$ is STEGR is general relativity, with no sign flip needed
anywhere. That is a happier accident than the torsion case, where GR-06's
$\mathcal{T} = +6H^2$ needs $f = -\mathcal{T}$ for the same statement.
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
| **GR-07** | a metric and a flat connection in, non-metricity and $f(Q)$ cosmology out |
| **GR-08** | the same quadratic actions in, kinetic matrices and wave speeds out |

With this one, every $G_{\text{eff}}$ the series uses is derived inside it. GR-08 then asks the
harder question of the same quadratic actions: not what the equations say, but whether the modes
appearing in them propagate at all.
::*)

(*::md::
## 1. Setup

Only the curvature backend is borrowed, and only to check the non-metricity results against
Levi-Civita ones computed independently.
::*)

(*::code::*)
ClearAll["Global`*"];

(* The volume element brings in Sqrt[-Det g] = Sqrt[a^6], and without knowing the
   scale factor is positive Wolfram will not reduce that to a^3, leaving factors of
   (a^3 - Sqrt[a^6]) sitting in expressions that are otherwise zero. *)
$Assumptions = a[t] > 0 && Nl[t] > 0;

$grMetric = ResourceFunction["MetricTensor"];
$grRicci  = ResourceFunction["RicciTensor"];

RicciScalarOf[gg_?MatrixQ, cs_List] := Module[{mt = $grMetric[gg, cs]},
  Simplify[Tr[mt["InverseMetricTensor"]["MatrixRepresentation"] .
     $grRicci[mt]["MatrixRepresentation"]]]];

xs = {t, x, y, z};

(*::md::
## 2. The non-metricity toolkit

Everything from its definition. `NonMetricity` takes a metric and a connection &mdash; the
coincident gauge is simply a connection of zeros &mdash; and returns the non-metricity tensor,
its two traces and the scalar.

`StueckelbergConnection` builds a general flat, torsion-free connection from four functions.
That it really is flat and torsion-free is checked rather than asserted: torsion is the
antisymmetric part, which vanishes by construction, and flatness is the vanishing of the
Riemann tensor built from $\Gamma$, which the next section tests on a non-trivial example.
::*)

(*::code::*)
ClearAll[NonMetricity, StueckelbergConnection, BoundaryQ];

NonMetricity[gg_?MatrixQ, gam_, cs_List] :=
 Module[{n = Length[cs], ig, qd, qu, qtr, qtil, qup, qtilup, l, m, i, p, q, al},
  ig = Inverse[gg];
  qd = Table[D[gg[[m, i]], cs[[l]]]
     - Sum[gam[[p, l, m]] gg[[p, i]] + gam[[p, l, i]] gg[[m, p]], {p, n}],
    {l, n}, {m, n}, {i, n}];
  qu = Table[Sum[ig[[l, p]] ig[[m, q]] ig[[i, al]] qd[[p, q, al]],
      {p, n}, {q, n}, {al, n}], {l, n}, {m, n}, {i, n}];
  qtr    = Table[Sum[ig[[m, i]] qd[[l, m, i]], {m, n}, {i, n}], {l, n}];
  qtil   = Table[Sum[ig[[m, i]] qd[[m, l, i]], {m, n}, {i, n}], {l, n}];
  qup    = Table[Sum[ig[[l, p]] qtr[[p]], {p, n}], {l, n}];
  qtilup = Table[Sum[ig[[l, p]] qtil[[p]], {p, n}], {l, n}];
  <|"Q" -> qd, "Trace" -> qtr, "TraceTilde" -> qtil, "Metric" -> gg, "Coordinates" -> cs,
    "Scalar" -> Simplify[
      -(1/4) Sum[qd[[l, m, i]] qu[[l, m, i]], {l, n}, {m, n}, {i, n}]
       + (1/2) Sum[qd[[l, m, i]] qu[[i, m, l]], {l, n}, {m, n}, {i, n}]
       + (1/4) Sum[qtr[[l]] qup[[l]], {l, n}]
       - (1/2) Sum[qtr[[l]] qtilup[[l]], {l, n}]]|>];

(* the general flat torsion-free connection, from four Stueckelberg fields *)
StueckelbergConnection[xi_List, cs_List] :=
 Module[{n = Length[cs], jinv, l, m, i, p},
  jinv = Inverse[Table[D[xi[[p]], cs[[m]]], {p, n}, {m, n}]];
  Table[Sum[jinv[[l, p]] D[xi[[p]], cs[[m]], cs[[i]]], {p, n}], {l, n}, {m, n}, {i, n}]];

(* B = nabla_mu (Q^mu - Qtilde^mu) *)
BoundaryQ[nn_Association] :=
 Module[{gg = nn["Metric"], cs = nn["Coordinates"], n, ig, sq, vec, m, p},
  n = Length[cs]; ig = Inverse[gg]; sq = Sqrt[-Det[gg]];
  vec = Table[Sum[ig[[m, p]] (nn["Trace"][[p]] - nn["TraceTilde"][[p]]), {p, n}], {m, n}];
  Simplify[(1/sq) Sum[D[sq vec[[m]], cs[[m]]], {m, n}]]];

$zeroConnection = ConstantArray[0, {4, 4, 4}];

(*::md::
## 3. $Q$ in flat FLRW, and the STEGR identity

In the coincident gauge the non-metricity is just $\partial_\lambda g_{\mu\nu}$, and for flat
FLRW the scalar comes out as $-6H^2$ &mdash; built from first derivatives of the metric, where
the Ricci scalar needs second ones. That is the whole of symmetric teleparallel gravity in one
observation, and it is the same observation the torsion scalar makes in GR-06.

The identity that makes STEGR equivalent to general relativity is

$$
R = Q - \nabla_\mu\!\left(Q^\mu - \tilde{Q}^\mu\right),
$$

and like every identity it has to hold for cases it was not designed for. The three below are a
homogeneous cosmology, the same with a lapse, and an inhomogeneous perturbed metric.
::*)

(*::code::*)
Module[{gF, gL, gP, cases},
 gF = DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}];
 gL = DiagonalMatrix[{-Nl[t]^2, a[t]^2, a[t]^2, a[t]^2}];
 gP = DiagonalMatrix[{-(1 + 2 uu[t, x]), a[t]^2 (1 - 2 vw[t, x]),
    a[t]^2 (1 - 2 vw[t, x]), a[t]^2 (1 - 2 vw[t, x])}];
 cases = {{"flat FLRW", gF}, {"FLRW with a lapse", gL}, {"perturbed FLRW", gP}};

 Dataset[
  Function[c,
   Module[{nn = NonMetricity[c[[2]], $zeroConnection, xs], rr},
    rr = RicciScalarOf[c[[2]], xs];
    <|"metric" -> c[[1]],
      "Q" -> nn["Scalar"],
      "R = Q - B" -> (Simplify[rr - (nn["Scalar"] - BoundaryQ[nn])] === 0)|>]] /@ cases]]

(*::md::
## 4. The connection is a field, and the coincident gauge is a choice

A flat torsion-free connection carries four functions, not more. The cell below builds one from
a non-trivial $\xi$, confirms it is torsion-free, confirms the identity still holds with a
connection switched on, and shows the thing that matters: the coincident gauge is reachable by
a coordinate change, so it is gauge, so it competes with Newtonian gauge for the same
freedom.

Under a diffeomorphism with parameter $\varepsilon^\mu$ the metric perturbation shifts by a Lie
derivative and $\delta\xi^\mu$ shifts by $-\varepsilon^\mu$. Spend the freedom setting
$\delta\xi^\mu = 0$ and you are in coincident gauge with a metric that is not Newtonian; spend
it the other way, as this notebook does, and $\delta\xi^\mu$ becomes four physical fields to be
varied.
::*)

(*::code::*)
Module[{xiT, gamT, nn, gF, tors},
 gF   = DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}];
 xiT  = {t + ff[t], x + gg1[x], y, z};
 gamT = Simplify[StueckelbergConnection[xiT, xs]];
 tors = Simplify[Table[gamT[[l, m, i]] - gamT[[l, i, m]], {l, 4}, {m, 4}, {i, 4}]];
 nn   = NonMetricity[gF, gamT, xs];
 Column[{
   Row[{Style["Gamma^0_00 = ", Bold], gamT[[1, 1, 1]]}],
   Row[{Style["Gamma^1_11 = ", Bold], gamT[[2, 2, 2]]}],
   Dataset @ {
     <|"check" -> "the Stueckelberg connection is torsion-free by construction",
       "ok" -> (tors === ConstantArray[0, {4, 4, 4}])|>,
     <|"check" -> "the identity R = Q - B survives a non-zero connection",
       "ok" -> (Simplify[RicciScalarOf[gF, xs] - (nn["Scalar"] - BoundaryQ[nn])] === 0)|>}}]]

(*::md::
## 5. Background cosmology

With the coincident gauge on the background and a lapse in the metric, minisuperspace variation
gives the $f(Q)$ Friedmann equations. They should be the ones GR-02 obtained from the metric
side, which is the same cross-check GR-06 ran for torsion.
::*)

(*::code::*)
ClearAll[eulerT];
eulerT[lag_, q_Symbol] := D[lag, q[t]] - D[D[lag, q'[t]], t];

Module[{gL, nn, lag, toN1, contin, names, fr1, fr2, qv},
 gL = DiagonalMatrix[{-Nl[t]^2, a[t]^2, a[t]^2, a[t]^2}];
 nn = NonMetricity[gL, $zeroConnection, xs];
 lag = Sqrt[-Det[gL]] (F[nn["Scalar"]]/(2 \[Kappa]) - dens[a[t]]);

 toN1   = {Derivative[_][Nl][t] -> 0, Nl[t] -> 1};
 contin = {Derivative[1][dens][a[t]] -> -3 (dens[a[t]] + pres[a[t]])/a[t]};
 names  = {dens[a[t]] -> \[Rho], pres[a[t]] -> p};
 qv = -6 a'[t]^2/a[t]^2;

 fr1 = Simplify[(eulerT[lag, Nl] /. toN1 /. contin /. names)/a[t]^3];
 fr2 = Simplify[(eulerT[lag, a]  /. toN1 /. contin /. names)/a[t]^2];

 Column[{
   Row[{Style["Q with a lapse : ", Bold], nn["Scalar"]}],
   Row[{Style["Friedmann I    : ", Bold], fr1, " = 0"}],
   Dataset @ {
    <|"check" -> "Friedmann I is kappa rho = (f - 2 Q f_Q)/2",
      "ok" -> (Simplify[fr1 - ((F[qv] - 2 qv F'[qv])/(2 \[Kappa]) - \[Rho])] === 0)|>,
    <|"check" -> "STEGR, f = Q, gives 3H^2 = kappa rho",
      "ok" -> (Simplify[(fr1 /. F -> Identity) - (3 a'[t]^2/(\[Kappa] a[t]^2) - \[Rho])] === 0)|>,
    <|"check" -> "STEGR gives the standard acceleration equation",
      "ok" -> (Simplify[((fr2 /. F -> Identity)
          /. First[Solve[(fr1 /. F -> Identity) == 0, \[Rho]]])
          - (3 (a'[t]^2 + 2 a[t] a''[t])/(\[Kappa] a[t]^2) + 3 p)] === 0)|>,
    <|"check" -> "no sign flip needed: this is already the GR-02 convention for Qs",
      "ok" -> (Simplify[qv + 6 a'[t]^2/a[t]^2] === 0)|>}}]]

(*::md::
## 6. The Euler-Lagrange operator, twice over

Two things have to be right here, and the first one is the trap GR-06 documents: the obvious
Wolfram spelling of the field-theory Euler-Lagrange operator can silently discard the spatial
term, which is exactly the one carrying $k^2$. It does so in WLJS and not under `wolframscript`,
on the same engine, for reasons neither notebook can account for &mdash; so the operator is
shown below but never asserted. The single-variable operator used for the background above is
unaffected, and both operators built here work by pattern replacement instead, which behaves
the same everywhere.

The second is specific to $f(Q)$ and easy to walk past. The non-metricity tensor contains the
connection, and the connection contains $\partial_\mu\partial_\nu\xi^\rho$, so
$\mathcal{L}^{(2)}$ carries **second** derivatives of the Stückelberg fields &mdash; the metric
perturbation appears only with first derivatives, but $\xi$ does not. An Euler-Lagrange operator
that stops at first derivatives therefore gets the connection equations wrong while leaving the
metric equations untouched, which is a nasty failure mode: the STEGR acceptance test lives in
the metric sector and passes either way.

`EulerLagrange` below runs to second order,

$$
\frac{\partial L}{\partial q} - \partial_\mu\frac{\partial L}{\partial(\partial_\mu q)} + \partial_\mu\partial_\nu\frac{\partial L}{\partial(\partial_\mu\partial_\nu q)},
$$

and the difference is visible immediately: with it, the STEGR connection equations come out
identically zero, which is what diffeomorphism invariance demands, where the first-order
operator leaves a spurious constraint behind.

The table below is a demonstration, not a checklist. Rows that come back zero are the point:
each shows an operator falling short on the test that motivated the next one. The assertions
live in the checks section, and they cover only `EulerLagrange2D` and `EulerLagrange`, whose
results do not vary between environments.
::*)

(*::code::*)
ClearAll[EulerLagrangeNaive, EulerLagrange2D];

EulerLagrangeNaive[lag_, q_] :=
  D[lag, q[t, x]] - D[D[lag, Derivative[1, 0][q][t, x]], t]
    - D[D[lag, Derivative[0, 1][q][t, x]], x];

EulerLagrange2D[lag_, q_] :=
 Module[{u0, u1, u2, fwd, back, ls},
  fwd  = {Derivative[1, 0][q][t, x] -> u1, Derivative[0, 1][q][t, x] -> u2, q[t, x] -> u0};
  back = {u0 -> q[t, x], u1 -> Derivative[1, 0][q][t, x], u2 -> Derivative[0, 1][q][t, x]};
  ls   = lag /. fwd;
  (D[ls, u0] /. back) - D[D[ls, u1] /. back, t] - D[D[ls, u2] /. back, x]];

(* the one actually used: second order, because Q carries d_mu d_nu xi *)
ClearAll[EulerLagrange];
EulerLagrange[lag_, q_] :=
 Module[{u, ord, fwd, back, ls},
  ord = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}};
  fwd = Table[If[o === {0, 0}, q[t, x] -> u[0, 0],
     Derivative[o[[1]], o[[2]]][q][t, x] -> u[o[[1]], o[[2]]]], {o, ord}];
  back = Table[If[o === {0, 0}, u[0, 0] -> q[t, x],
     u[o[[1]], o[[2]]] -> Derivative[o[[1]], o[[2]]][q][t, x]], {o, ord}];
  ls = lag /. fwd;
  (D[ls, u[0, 0]] /. back)
   - D[D[ls, u[1, 0]] /. back, t] - D[D[ls, u[0, 1]] /. back, x]
   + D[D[ls, u[2, 0]] /. back, {t, 2}] + D[D[ls, u[1, 1]] /. back, t, x]
   + D[D[ls, u[0, 2]] /. back, {x, 2}]];

(* A demonstration, not a checklist. The answers by hand are -d_x^2 B for the
   first test and d_t^2 B for the second. The naive row carries no verdict on
   purpose: what it returns depends on the environment. *)
Module[{toy, toy2},
 toy  = D[uA[t, x], x] D[uB[t, x], x];
 toy2 = D[uA[t, x], {t, 2}] uB[t, x];
 Dataset @ {
   <|"operator" -> "naive, first order", "test" -> "(d_x A)(d_x B)",
     "result" -> EulerLagrangeNaive[toy, uA],
     "note" -> "environment dependent: zero in WLJS, correct under wolframscript"|>,
   <|"operator" -> "EulerLagrange2D", "test" -> "(d_x A)(d_x B)",
     "result" -> EulerLagrange2D[toy, uA],
     "note" -> "correct"|>,
   <|"operator" -> "EulerLagrange2D", "test" -> "(d_t^2 A) B",
     "result" -> EulerLagrange2D[toy2, uA],
     "note" -> "blind to second derivatives, which is why the next row exists"|>,
   <|"operator" -> "EulerLagrange", "test" -> "(d_t^2 A) B",
     "result" -> EulerLagrange[toy2, uA],
     "note" -> "correct"|>}]

(*::md::
## 7. The second order action, in fourteen free functions

Ten metric perturbations and four Stückelberg perturbations, all free functions of $t$ and $x$,
truncated at $O(\epsilon^2)$ at every step so that no general symbolic metric is ever inverted.
The linear field equations are the Euler-Lagrange equations of $\mathcal{L}^{(2)}$, and the
perturbations stay general until after the variation &mdash; a plane wave substituted early
would put $e^{2ikx}$ into a quadratic Lagrangian and the phases would stop cancelling.

To second order the connection is
$\Gamma^\lambda{}_{\mu\nu} = \epsilon\,\partial_\mu\partial_\nu\zeta^\lambda - \epsilon^2\,\partial_\rho\zeta^\lambda\,\partial_\mu\partial_\nu\zeta^\rho$.
::*)

(*::code::*)
ClearAll[trunc2, lagQ1, lagQ2, $Q0];

trunc2[e_] := Normal[Series[e, {eps, 0, 2}]];

Module[{bgG, ibgG, hm, gfull, igfull, mmg, sqrtg, zt, gam, qd, qu, qtr, qtil,
        qup, qtilup, qScalar, dpar, lg},
 bgG  = DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}];
 ibgG = DiagonalMatrix[{-1, 1/a[t]^2, 1/a[t]^2, 1/a[t]^2}];
 hm   = Table[hq[Min[i, j], Max[i, j]][t, x], {i, 4}, {j, 4}];   (* symmetric: ten functions *)

 gfull  = bgG + eps hm;
 igfull = ibgG - eps ibgG . hm . ibgG + eps^2 ibgG . hm . ibgG . hm . ibgG;
 mmg    = ibgG . hm;
 sqrtg  = Sqrt[-Det[bgG]] (1 + eps Tr[mmg]/2 + eps^2 (Tr[mmg]^2 - 2 Tr[mmg . mmg])/8);

 dpar[f_, mu_] := If[mu <= 2, D[f, {t, x}[[mu]]], 0];

 zt  = Table[zz[l][t, x], {l, 4}];
 gam = Table[eps dpar[dpar[zt[[l]], m], i]
     - eps^2 Sum[dpar[zt[[l]], p] dpar[dpar[zt[[p]], m], i], {p, 4}],
    {l, 4}, {m, 4}, {i, 4}];

 qd = Table[trunc2[dpar[gfull[[m, i]], l]
     - Sum[gam[[p, l, m]] gfull[[p, i]] + gam[[p, l, i]] gfull[[m, p]], {p, 4}]],
    {l, 4}, {m, 4}, {i, 4}];
 qu = Table[trunc2[Sum[igfull[[l, p]] igfull[[m, q]] igfull[[i, al]] qd[[p, q, al]],
      {p, 4}, {q, 4}, {al, 4}]], {l, 4}, {m, 4}, {i, 4}];
 qtr    = Table[trunc2[Sum[igfull[[m, i]] qd[[l, m, i]], {m, 4}, {i, 4}]], {l, 4}];
 qtil   = Table[trunc2[Sum[igfull[[m, i]] qd[[m, l, i]], {m, 4}, {i, 4}]], {l, 4}];
 qup    = Table[trunc2[Sum[igfull[[l, p]] qtr[[p]], {p, 4}]], {l, 4}];
 qtilup = Table[trunc2[Sum[igfull[[l, p]] qtil[[p]], {p, 4}]], {l, 4}];

 qScalar = trunc2[
    -(1/4) Sum[qd[[l, m, i]] qu[[l, m, i]], {l, 4}, {m, 4}, {i, 4}]
     + (1/2) Sum[qd[[l, m, i]] qu[[i, m, l]], {l, 4}, {m, 4}, {i, 4}]
     + (1/4) Sum[qtr[[l]] qup[[l]], {l, 4}]
     - (1/2) Sum[qtr[[l]] qtilup[[l]], {l, 4}]];

 $Q0 = Simplify[qScalar /. eps -> 0];
 lg  = trunc2[sqrtg F[qScalar]];
 lagQ1 = Coefficient[lg, eps, 1];
 lagQ2 = Coefficient[lg, eps, 2];

 Dataset @ {
  <|"check" -> "the series inverse metric is exact to O(eps^2)",
    "ok" -> (Simplify[trunc2[gfull . igfull] - IdentityMatrix[4]] === ConstantArray[0, {4, 4}])|>,
  <|"check" -> "the series volume element is exact to O(eps^2)",
    "ok" -> (Simplify[trunc2[Sqrt[-Det[gfull]]] - sqrtg] === 0)|>,
  <|"check" -> "background Q = -6 H^2",
    "ok" -> (Simplify[$Q0 + 6 a'[t]^2/a[t]^2] === 0)|>,
  <|"check" -> "the second order Lagrangian carries metric gradients",
    "ok" -> (Length[Cases[lagQ2, Derivative[0, 1][hq[__]][t, x], Infinity]] > 0)|>,
  <|"check" -> "and carries the Stueckelberg fields",
    "ok" -> (Length[Cases[lagQ2, zz[_][t, x] | Derivative[__][zz[_]][t, x], Infinity]] > 0)|>}]

(*::md::
## 8. The scalar ansatz and the matter source

Newtonian gauge on the metric, and the two scalar connection modes kept:

| | |
|---|---|
| $\psi$ | $g_{00} = -(1+2\psi)$ |
| $\phi$ | $g_{ij} = a^2(1-2\phi)\delta_{ij}$ |
| $z_0$ | $\delta\xi^0$, a scalar connection mode |
| $z_1$ | $\delta\xi^1$, the other one |

Matter couples to the metric and **not** to the connection, so the four connection equations
carry no source at all. That is not an assumption, it is what $S_m[g]$ gives, and it is what
makes the system close: the connection sector is pure constraint.

The constant in $\delta S_m/\delta g_{\mu\nu} = \lambda\sqrt{-g}\,T^{\mu\nu}$ is fixed by
demanding the background reproduce section 5, and comes out $\lambda = -\kappa$.
::*)

(*::code::*)
ClearAll[metVars, conVars, qAnsatz, matterSource, backgroundOnShell];

metVars = Flatten[Table[If[i <= j, hq[i, j], Nothing], {i, 4}, {j, 4}]];
conVars = Table[zz[l], {l, 4}];

qAnsatz = Join[
  {hq[1, 1] -> Function[{t, x}, -2 psi[t] Exp[I k x]],
   hq[2, 2] -> Function[{t, x}, -2 a[t]^2 phi[t] Exp[I k x]],
   hq[3, 3] -> Function[{t, x}, -2 a[t]^2 phi[t] Exp[I k x]],
   hq[4, 4] -> Function[{t, x}, -2 a[t]^2 phi[t] Exp[I k x]]},
  (# -> Function[{t, x}, 0] &) /@ {hq[1, 2], hq[1, 3], hq[1, 4], hq[2, 3], hq[2, 4], hq[3, 4]},
  {zz[1] -> Function[{t, x}, z0[t] Exp[I k x]],
   zz[2] -> Function[{t, x}, z1[t] Exp[I k x]],
   zz[3] -> Function[{t, x}, 0], zz[4] -> Function[{t, x}, 0]}];

Module[{ph, trunc1, gpert, uUp, rhoT, tUp, sqg, bgMet, b00, b11},
 ph = Exp[I k x];
 trunc1[e_] := Normal[Series[e, {eps, 0, 1}]];
 gpert = DiagonalMatrix[{-(1 + 2 eps psi[t] ph), a[t]^2 (1 - 2 eps phi[t] ph),
     a[t]^2 (1 - 2 eps phi[t] ph), a[t]^2 (1 - 2 eps phi[t] ph)}];
 uUp   = {1 - eps psi[t] ph, eps vv[t] ph, 0, 0};
 rhoT  = rho[t] (1 + eps del[t] ph);
 tUp   = Table[trunc1[rhoT uUp[[i]] uUp[[j]]], {i, 4}, {j, 4}];
 sqg   = trunc1[Sqrt[-Det[gpert]]];

 matterSource = Association[
    Flatten@Table[hq[i, j] -> trunc1[(-\[Kappa]) sqg tUp[[i, j]] If[i === j, 1, 2]],
      {i, 4}, {j, i, 4}]];

 bgMet = Association[# -> EulerLagrange[lagQ1, #] & /@ metVars];
 b00 = Simplify[bgMet[hq[1, 1]] - (matterSource[hq[1, 1]] /. eps -> 0)];
 b11 = Simplify[bgMet[hq[2, 2]] - (matterSource[hq[2, 2]] /. eps -> 0)];
 backgroundOnShell = Simplify[First@Solve[{b00 == 0, b11 == 0}, {rho[t], a''[t]}]];

 Column[{
   Row[{Style["background 0-0 : ", Bold], Short[b00, 2], " = 0"}],
   Dataset @ {
    <|"check" -> "the four-velocity is normalised to O(eps)",
      "ok" -> (Simplify[trunc1[Sum[gpert[[i, j]] uUp[[i]] uUp[[j]], {i, 4}, {j, 4}]] + 1] === 0)|>,
    <|"check" -> "lambda = -kappa reproduces the section 5 Friedmann equation",
      "ok" -> (Simplify[(b00 /. F -> Identity)
         /. First@Solve[3 a'[t]^2/a[t]^2 == \[Kappa] rho[t], rho[t]]] === 0)|>,
    <|"check" -> "matter does not source the connection equations",
      "ok" -> (! MemberQ[Keys[matterSource], zz[_]])|>}}]]

(*::md::
## 9. The STEGR acceptance test

Vary fourteen ways, substitute, subtract the source, drop time derivatives of the amplitudes,
and put the background on shell. Then set $f = Q$ and demand general relativity back:

- the $00$ equation must be the Poisson equation, carrying $k^2$;
- the transverse equation must force $\Phi = \Psi$;
- the connection modes must drop out of the metric equations entirely, because in STEGR the
  $\xi$ dependence of $\sqrt{-g}\,Q$ is a total derivative;
- the pressure equation must be consistent with dust.

If any of that failed, nothing after it would be worth reading.
::*)

(*::code::*)
ClearAll[linMet, linCon, quasiStatic];

quasiStatic = {Derivative[_][psi][t] -> 0, Derivative[_][phi][t] -> 0,
               Derivative[_][z0][t] -> 0, Derivative[_][z1][t] -> 0,
               Derivative[_][del][t] -> 0, Derivative[_][vv][t] -> 0};

Module[{ph, subst, onShell, gMet, gCon},
 ph = Exp[I k x];
 subst[e_] := Simplify[(e /. qAnsatz)/ph];
 onShell[e_] := Simplify[Simplify[e /. quasiStatic] /. backgroundOnShell];

 gMet = Association[# -> subst[EulerLagrange[lagQ2, #]] & /@ metVars];
 gCon = Association[# -> subst[EulerLagrange[lagQ2, #]] & /@ conVars];

 linMet = Association[# -> onShell[gMet[#] - Coefficient[matterSource[#], eps, 1]/ph]
     & /@ {hq[1, 1], hq[1, 2], hq[2, 2], hq[3, 3], hq[4, 4]}];
 linCon = Association[# -> onShell[gCon[#]] & /@ {zz[1], zz[2]}];
];

Module[{st = {F -> Identity}, m00, m0x, mxx, myy, c1, c2},
 m00 = Simplify[linMet[hq[1, 1]] /. st]; m0x = Simplify[linMet[hq[1, 2]] /. st];
 mxx = Simplify[linMet[hq[2, 2]] /. st]; myy = Simplify[linMet[hq[3, 3]] /. st];
 c1 = Simplify[linCon[zz[1]] /. st];     c2 = Simplify[linCon[zz[2]] /. st];
 Column[{
  Row[{Style["STEGR  0-0       : ", Bold], m00, " = 0"}],
  Row[{Style["STEGR  transverse: ", Bold], myy, " = 0"}],
  Row[{Style["STEGR  xx        : ", Bold], mxx, " = 0"}],
  "",
  Dataset @ {
   <|"acceptance test" -> "0-0 is the Poisson equation, carrying k^2",
     "ok" -> (Exponent[Expand[m00], k] === 2)|>,
   <|"acceptance test" -> "the transverse equation forces phi = psi",
     "ok" -> (myy =!= 0 && Simplify[myy /. phi[t] -> psi[t]] === 0)|>,
   <|"acceptance test" -> "the connection modes leave the metric equations",
     "ok" -> AllTrue[{m00, m0x, mxx, myy}, FreeQ[#, z0 | z1] &]|>,
   <|"acceptance test" -> "the dust pressure equation is consistent",
     "ok" -> (mxx === 0 || Simplify[mxx /. phi[t] -> psi[t]] === 0)|>,
   <|"acceptance test" -> "the connection equations vanish identically, as diffeomorphism invariance demands",
     "ok" -> (Simplify[c1] === 0 && Simplify[c2] === 0)|>}}]]

(*::md::
## 10. Closing the system

Five scalar unknowns &mdash; $\psi$, $\phi$, $z_0$, $z_1$ and the peculiar velocity $v$ &mdash;
and five equations. The $0x$ equation gives $v$; the $xx$, transverse and connection equations
are all source free and fix $\phi$, $z_0$ and $z_1$ in terms of $\psi$; the $00$ equation is
then read off against

$$\frac{k^2}{a^2}\Psi = -4\pi G_{\text{eff}}\,\rho_m\delta_m .$$

The $xx$ equation is the one it is tempting to leave out, and leaving it out is exactly what
keeps $z_0$ undetermined and the answer contaminated by it. A result that still depends on a
connection mode is a result that has not closed.
::*)

(*::code::*)
ClearAll[geffExact, geffSubHorizon, slipSubHorizon];

Module[{abbrev, nm, e00, e0x, exx, eyy, c1, vsol, sol, poi, dsol},
 abbrev = {Derivative[3][F][$Q0] -> fQQQ, Derivative[2][F][$Q0] -> fQQ,
           Derivative[1][F][$Q0] -> fQ, F[$Q0] -> f0};
 nm[e_] := Simplify[Simplify[e] /. abbrev];

 e00 = nm[linMet[hq[1, 1]]]; e0x = nm[linMet[hq[1, 2]]];
 exx = nm[linMet[hq[2, 2]]]; eyy = nm[linMet[hq[3, 3]]];
 c1  = nm[linCon[zz[1]]];

 vsol = First@Solve[e0x == 0, vv[t]];
 sol  = Simplify[First@Solve[
    {(exx /. vsol) == 0, (eyy /. vsol) == 0, (c1 /. vsol) == 0}, {phi[t], z0[t], z1[t]}]];
 poi  = Simplify[(e00 /. vsol) /. sol];
 dsol = First@Solve[poi == 0, del[t]];

 geffExact = Simplify[
   ((k^2 psi[t]/a[t]^2)/(-(\[Kappa]/2) rho[t] (del[t] /. dsol)))
    /. backgroundOnShell /. abbrev];
 geffSubHorizon = Simplify[Limit[geffExact, k -> Infinity]];
 slipSubHorizon = Simplify[
   Limit[Simplify[((phi[t] /. sol)/psi[t]) /. backgroundOnShell /. abbrev], k -> Infinity]];

 Column[{
   Row[{Style["G_eff / G, sub-horizon     : ", Bold], geffSubHorizon}],
   Row[{Style["slip  Phi/Psi, sub-horizon : ", Bold], slipSubHorizon}]}]]

(*::md::
## 11. What it says

$$
\frac{G_{\text{eff}}}{G} = \frac{1}{f_Q}, \qquad \eta = \frac{\Phi}{\Psi} = 1
$$

with no $k$ anywhere. Gravity is rescaled by a single time-dependent factor and nothing else
happens: growth and lensing feel the same modification, and every wavenumber feels it equally.

Set against the other two families the pattern is clean.

| | $G_{\text{eff}}/G$ | slip | scale dependent |
|---|---|---|---|
| $f(R)$, GR-05 | $\dfrac{1}{f_R}\dfrac{1+4m}{1+3m}$ | $\dfrac{1+2m}{1+4m}$ | **yes**, running to $4/3$ |
| $f(\mathcal{T})$, GR-06 | $-\dfrac{1}{f_{\mathcal{T}}}$ | $1$ | no |
| $f(Q)$, here | $\dfrac{1}{f_Q}$ | $1$ | no |

$f(R)$ carries a scalaron with a Compton wavelength, so it has a scale to compare $k$ against;
the two teleparallel families do not, and their modification is a pure rescaling. That is the
sharpest structural difference between them, and it is the one a survey covering a range of
scales could actually see.

The expression is exactly what GR-04's library has been quoting for $f(Q)$. Nothing in the
series is quoted now.
::*)

(*::code::*)
Dataset @ {
  <|"theory" -> "f(R)",  "G_eff/G" -> "(1/fR)(1+4m)/(1+3m)", "slip" -> "(1+2m)/(1+4m)",
    "scale dependent" -> True,  "derived in" -> "GR-05"|>,
  <|"theory" -> "f(T)",  "G_eff/G" -> "-1/fT", "slip" -> "1",
    "scale dependent" -> False, "derived in" -> "GR-06"|>,
  <|"theory" -> "f(Q)",  "G_eff/G" -> "1/fQ",  "slip" -> "1",
    "scale dependent" -> False, "derived in" -> "GR-07"|>}

(*::md::
## 12. Checks

Geometry, background, and then the perturbative result. The last group is the one that matters.
::*)

(*::code::*)
Module[{gF, gL, gP, nF, nL, nP, qv, lag, fr1, toN1, contin, names, stg},
 gF = DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}];
 gL = DiagonalMatrix[{-Nl[t]^2, a[t]^2, a[t]^2, a[t]^2}];
 gP = DiagonalMatrix[{-(1 + 2 uu[t, x]), a[t]^2 (1 - 2 vw[t, x]),
    a[t]^2 (1 - 2 vw[t, x]), a[t]^2 (1 - 2 vw[t, x])}];
 nF = NonMetricity[gF, $zeroConnection, xs];
 nL = NonMetricity[gL, $zeroConnection, xs];
 nP = NonMetricity[gP, $zeroConnection, xs];
 qv = -6 a'[t]^2/a[t]^2;
 lag = Sqrt[-Det[gL]] (F[nL["Scalar"]]/(2 \[Kappa]) - dens[a[t]]);
 toN1   = {Derivative[_][Nl][t] -> 0, Nl[t] -> 1};
 contin = {Derivative[1][dens][a[t]] -> -3 (dens[a[t]] + pres[a[t]])/a[t]};
 names  = {dens[a[t]] -> \[Rho], pres[a[t]] -> p};
 fr1 = Simplify[(eulerT[lag, Nl] /. toN1 /. contin /. names)/a[t]^3];
 stg = {f0 -> $Q0, fQ -> 1, fQQ -> 0, fQQQ -> 0};

 Dataset @ {
  <|"check" -> "Q = -6 H^2 in flat FLRW",
    "ok" -> (Simplify[nF["Scalar"] - qv] === 0)|>,
  <|"check" -> "R = Q - B on three metrics",
    "ok" -> AllTrue[{{nF, gF}, {nL, gL}, {nP, gP}},
       Simplify[RicciScalarOf[#[[2]], xs] - (#[[1]]["Scalar"] - BoundaryQ[#[[1]]])] === 0 &]|>,
  <|"check" -> "Friedmann I is kappa rho = (f - 2 Q f_Q)/2",
    "ok" -> (Simplify[fr1 - ((F[qv] - 2 qv F'[qv])/(2 \[Kappa]) - \[Rho])] === 0)|>,
  <|"check" -> "STEGR background gives 3H^2 = kappa rho",
    "ok" -> (Simplify[(fr1 /. F -> Identity) - (3 a'[t]^2/(\[Kappa] a[t]^2) - \[Rho])] === 0)|>,
  <|"check" -> "EulerLagrange2D gets the gradient term right",
    "ok" -> (Simplify[EulerLagrange2D[D[uA[t, x], x] D[uB[t, x], x], uA]
       + D[uB[t, x], {x, 2}]] === 0)|>,
  <|"check" -> "L2 really does carry second time derivatives of the connection",
    "ok" -> (Length[Cases[lagQ2, Derivative[2, _][zz[_]][t, x], Infinity]] > 0)|>,
  <|"check" -> "so a first-order operator is not enough, and the second-order one is right",
    "ok" -> (Simplify[EulerLagrange[D[uA[t, x], {t, 2}] uB[t, x], uA]
          - D[uB[t, x], {t, 2}]] === 0 &&
       Simplify[EulerLagrange2D[D[uA[t, x], {t, 2}] uB[t, x], uA]
          - D[uB[t, x], {t, 2}]] =!= 0)|>,
  <|"check" -> "STEGR: the connection equations vanish identically",
    "ok" -> (Simplify[linCon[zz[1]] /. F -> Identity] === 0 &&
             Simplify[linCon[zz[2]] /. F -> Identity] === 0)|>,
  <|"check" -> "STEGR: the 0-0 equation carries k^2",
    "ok" -> (Exponent[Expand[Simplify[linMet[hq[1, 1]] /. F -> Identity]], k] === 2)|>,
  <|"check" -> "STEGR: no slip",
    "ok" -> Module[{e = Simplify[linMet[hq[3, 3]] /. F -> Identity]},
      e =!= 0 && Simplify[e /. phi[t] -> psi[t]] === 0]|>,
  <|"check" -> "STEGR: the connection modes leave the metric equations",
    "ok" -> AllTrue[{hq[1, 1], hq[1, 2], hq[2, 2], hq[3, 3]},
       FreeQ[Simplify[linMet[#] /. F -> Identity], z0 | z1] &]|>,
  <|"check" -> "the closed system no longer depends on the connection modes",
    "ok" -> FreeQ[geffSubHorizon, z0 | z1]|>,
  <|"check" -> "G_eff = 1/fQ, scale independent",
    "ok" -> (Simplify[geffSubHorizon - 1/fQ] === 0 && FreeQ[geffSubHorizon, k])|>,
  <|"check" -> "STEGR gives G_eff = G",
    "ok" -> (Simplify[geffSubHorizon /. stg] === 1)|>,
  <|"check" -> "no gravitational slip",
    "ok" -> (Simplify[slipSubHorizon - 1] === 0)|>}]

(*::md::
## 13. What is left

Nothing in the series is quoted any more. Every $G_{\text{eff}}$ in GR-04's library is derived:
$f(R)$ in GR-05, $f(\mathcal{T})$ in GR-06, $f(Q)$ here.

What remains are the approximations, and they are shared by all three. Quasi-static and
sub-horizon both fail on scales approaching the horizon, so none of these expressions should be
trusted for the largest modes a survey measures. And both teleparallel families carry known
perturbative pathologies: extra modes propagate, and their kinetic terms can degenerate on
cosmological backgrounds, which is the strong-coupling problem. A correct linear calculation is
necessary but not sufficient, and this notebook has only done the linear calculation.
::*)

(*::code::*)
Dataset @ <|
  "derived, all of it" -> "f(R) in GR-05, f(T) in GR-06, f(Q) here",
  "still approximate"  -> "quasi-static and sub-horizon; both fail near the horizon",
  "not addressed"      -> "strong coupling, where linear theory may not describe the extra modes",
  "cross-checks"       -> "background matches GR-02; the STEGR limit reproduces general relativity"|>
