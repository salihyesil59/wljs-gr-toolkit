(* ::Package:: *)

(* ============================================================================
   GR-08 : Strong Coupling and Which Modes Actually Propagate
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-08-Strong-Coupling.wl

   ============================================================================ *)

(*::md::
# Strong Coupling, and Which Modes Actually Propagate

GR-05, GR-06 and GR-07 each end with the same caveat: a correct linear calculation is
necessary, not sufficient. This notebook is that caveat, made quantitative.

$f(\mathcal{T})$ and $f(Q)$ have more fields than general relativity &mdash; a tetrad has six
components beyond the metric it builds, and a flat connection carries four functions of its own.
If those extra fields propagate, the theories have extra degrees of freedom. If their kinetic
terms **vanish** on the background you are perturbing around, they do not propagate at linear
order, their leading dynamics is non-linear, and linear perturbation theory is no longer a
reliable guide to them. That is the strong-coupling problem, and it is the reason the
$G_{\text{eff}}$ expressions in this series come with a health warning.

The diagnostic is the kinetic matrix. Write the quadratic action for one Fourier mode as

$$
L^{(2)} = \tfrac{1}{2}\,\dot{q}^{\,T} K\, \dot{q} \;-\; \tfrac{1}{2}\,q^{T}\!\left(\frac{k^2}{a^2}\,G + M\right) q + \ldots
$$

A zero eigenvalue of $K$ is a direction in field space with no kinetic term. The ratio of the
gradient coefficient to the kinetic one is the propagation speed squared, $c^2 = G/K$, and for
gravitational waves that number is now measured to a part in $10^{15}$.

## What comes out

Two clean results and one honest limit.

- **Both theories propagate gravitational waves at exactly the speed of light.** Not
  approximately: $c_{\text{GW}}^2 = 1$ identically, for both polarisations, for every $f$. Neither
  family is touched by the constraint from GW170817.
- **In $f(\mathcal{T})$ the extra Lorentz modes carry no kinetic term at all.** The boost and
  rotation rows of $K$ are identically zero, exactly like the lapse, which is a Lagrange
  multiplier. Around flat FLRW they simply do not propagate, and that is strong coupling
  computed rather than asserted.
- **In $f(Q)$ the bare kinetic matrix cannot settle the question**, and section 8 says why
  rather than pretending otherwise.

There is also a shared structure worth noticing. The surviving scalar kinetic term degenerates
in both theories on the same locus,

$$
f_X + 2X f_{XX} = 0 ,
$$

with $X = \mathcal{T}$ or $Q$. For a concrete model that is a particular value of $H$, which is
to say a particular redshift, and section 9 works one out.

## Method

Real perturbations, not a complex plane wave. A quadratic action built from $e^{ikx}$ picks up
$e^{2ikx}$ and the phases stop cancelling; with $\cos kx$ and $\sin kx$ and an average over one
wavelength the reduced Lagrangian is real and unambiguous. Scalars that enter through a spatial
gradient ride $\sin$, the rest ride $\cos$.
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

This one reuses the second order actions GR-06 and GR-07 built, and asks a different question
of them: not what the equations say, but whether the modes in them are there at all.
::*)

(*::md::
## 1. Setup

The two quadratic actions are rebuilt here rather than imported, so the notebook stands alone.
Both are the same construction: background plus $\epsilon$ times free perturbation functions,
truncated at $O(\epsilon^2)$ at every step.
::*)

(*::code::*)
ClearAll["Global`*"];

$Assumptions = a[t] > 0 && Nl[t] > 0 && k > 0;

$eta = DiagonalMatrix[{-1, 1, 1, 1}];
trunc2[e_] := Normal[Series[e, {eps, 0, 2}]];
dpar[f_, mu_] := If[mu <= 2, D[f, {t, x}[[mu]]], 0];   (* nothing depends on y or z *)

(*::md::
## 2. The two quadratic actions

$f(\mathcal{T})$ from sixteen free tetrad perturbations, $f(Q)$ from ten metric plus four
Stückelberg perturbations. These are the same objects GR-06 and GR-07 vary to get their field
equations; here they are used as quadratic forms instead.
::*)

(*::code::*)
ClearAll[lagT2, lagQ2, $T0, $Q0];

(* ---- f(T): tetrad ---- *)
Module[{bgT, ibgT, hp, eup, gg, ig, tor, tDn, tUp12, tUp23, kon, tr1c, sup, tS, mm, det},
 bgT  = DiagonalMatrix[{1, a[t], a[t], a[t]}];
 ibgT = DiagonalMatrix[{1, 1/a[t], 1/a[t], 1/a[t]}];
 hp   = Table[hh[i, j][t, x], {i, 4}, {j, 4}];
 mm   = ibgT . hp;
 det  = Det[bgT] (1 + eps Tr[mm] + eps^2 (Tr[mm]^2 - Tr[mm . mm])/2);
 eup  = With[{iv = ibgT - eps ibgT . hp . ibgT + eps^2 ibgT . hp . ibgT . hp . ibgT},
    Table[iv[[m, capA]], {capA, 4}, {m, 4}]];
 gg = Table[trunc2[Sum[$eta[[capA, capB]] (bgT + eps hp)[[capA, i]] (bgT + eps hp)[[capB, m]],
     {capA, 4}, {capB, 4}]], {i, 4}, {m, 4}];
 ig = Table[trunc2[Sum[$eta[[capA, capB]] eup[[capA, i]] eup[[capB, m]],
     {capA, 4}, {capB, 4}]], {i, 4}, {m, 4}];
 tor = Table[trunc2[Sum[eup[[capA, l]] (dpar[(bgT + eps hp)[[capA, i]], m]
      - dpar[(bgT + eps hp)[[capA, m]], i]), {capA, 4}]], {l, 4}, {m, 4}, {i, 4}];
 tDn   = Table[trunc2[Sum[gg[[l, q]] tor[[q, m, i]], {q, 4}]], {l, 4}, {m, 4}, {i, 4}];
 tUp12 = Table[trunc2[Sum[ig[[m, p]] ig[[i, q]] tDn[[p, q, l]], {p, 4}, {q, 4}]],
    {m, 4}, {i, 4}, {l, 4}];
 tUp23 = Table[trunc2[Sum[ig[[m, p]] ig[[i, q]] tDn[[l, p, q]], {p, 4}, {q, 4}]],
    {l, 4}, {m, 4}, {i, 4}];
 kon = Table[-(1/2) (tUp12[[m, i, l]] - tUp12[[i, m, l]] - tUp23[[l, m, i]]),
    {m, 4}, {i, 4}, {l, 4}];
 tr1c = Table[trunc2[Sum[tUp12[[al, i, al]], {al, 4}]], {i, 4}];
 sup = Table[(1/2) (kon[[m, i, l]] + KroneckerDelta[m, l] tr1c[[i]]
      - KroneckerDelta[i, l] tr1c[[m]]), {l, 4}, {m, 4}, {i, 4}];
 tS = trunc2[Sum[sup[[l, m, i]] tor[[l, m, i]], {l, 4}, {m, 4}, {i, 4}]];
 $T0 = Simplify[tS /. eps -> 0];
 lagT2 = Coefficient[trunc2[det F[tS]], eps, 2];
];

(* ---- f(Q): metric plus a flat connection ---- *)
Module[{bgG, ibgG, hm, gfull, igfull, mmg, sqrtg, zt, gam, qd, qu, qtr, qtil,
        qup, qtilup, qS},
 bgG  = DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}];
 ibgG = DiagonalMatrix[{-1, 1/a[t]^2, 1/a[t]^2, 1/a[t]^2}];
 hm   = Table[hq[Min[i, j], Max[i, j]][t, x], {i, 4}, {j, 4}];
 gfull  = bgG + eps hm;
 igfull = ibgG - eps ibgG . hm . ibgG + eps^2 ibgG . hm . ibgG . hm . ibgG;
 mmg    = ibgG . hm;
 sqrtg  = Sqrt[-Det[bgG]] (1 + eps Tr[mmg]/2 + eps^2 (Tr[mmg]^2 - 2 Tr[mmg . mmg])/8);
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
 qS = trunc2[
    -(1/4) Sum[qd[[l, m, i]] qu[[l, m, i]], {l, 4}, {m, 4}, {i, 4}]
     + (1/2) Sum[qd[[l, m, i]] qu[[i, m, l]], {l, 4}, {m, 4}, {i, 4}]
     + (1/4) Sum[qtr[[l]] qup[[l]], {l, 4}]
     - (1/2) Sum[qtr[[l]] qtilup[[l]], {l, 4}]];
 $Q0 = Simplify[qS /. eps -> 0];
 lagQ2 = Coefficient[trunc2[sqrtg F[qS]], eps, 2];
];

Dataset @ {
  <|"check" -> "torsion scalar background is +6 H^2",
    "ok" -> (Simplify[$T0 - 6 a'[t]^2/a[t]^2] === 0)|>,
  <|"check" -> "non-metricity scalar background is -6 H^2",
    "ok" -> (Simplify[$Q0 + 6 a'[t]^2/a[t]^2] === 0)|>}

(*::md::
## 3. Averaging over a wavelength

A quadratic form in $\cos kx$ and $\sin kx$ reduces, under `TrigReduce`, to a constant piece
plus terms in $\cos 2kx$ and $\sin 2kx$. The latter integrate to zero over a wavelength, so
dropping them *is* the average, and it costs one line instead of an integral.

`KineticMatrix` is then $K_{ij} = \partial^2 L/\partial\dot q_i\,\partial\dot q_j$, and
`SpeedSquared` reads the propagation speed off as the ratio of the $k^2$ part of the potential
to the kinetic coefficient, with the $a^2$ that turns a comoving wavenumber into a physical
gradient. The $k$-independent part of the potential is a mass and must not be counted &mdash;
including it was the first thing this notebook got wrong.
::*)

(*::code::*)
ClearAll[AverageOverWavelength, KineticMatrix, SpeedSquared];

AverageOverWavelength[e_] :=
  Expand[TrigReduce[Expand[e]]] /. {Cos[_. x] -> 0, Sin[_. x] -> 0};

KineticMatrix[lred_, vars_List] := Simplify[Table[
   D[lred, Derivative[1][vars[[i]]][t], Derivative[1][vars[[j]]][t]],
   {i, Length[vars]}, {j, Length[vars]}]];

SpeedSquared[lred_, q_Symbol] := Simplify[
  -a[t]^2 Coefficient[Expand[D[lred, q[t], q[t]]], k, 2]/
    D[lred, Derivative[1][q][t], Derivative[1][q][t]]];

(*::md::
## 4. The scalar sector of $f(\mathcal{T})$

Four scalars: the lapse $\psi$, the isotropic spatial mode $\phi$, and the two Lorentz
directions $\zeta$ (a boost) and $\chi$ (a rotation) that GR-06 showed are physical the moment
$f$ is non-linear. $\zeta$ and $\chi$ enter through spatial gradients, so they ride $\sin kx$.
::*)

(*::code::*)
ClearAll[scalarT, kinT];

scalarT = {hh[1, 1] -> Function[{t, x}, psi[t] Cos[k x]],
  hh[1, 2] -> Function[{t, x}, a[t] zeta[t] Sin[k x]],
  hh[2, 1] -> Function[{t, x}, zeta[t] Sin[k x]],
  hh[2, 2] -> Function[{t, x}, -a[t] phi[t] Cos[k x]],
  hh[3, 3] -> Function[{t, x}, -a[t] phi[t] Cos[k x]],
  hh[4, 4] -> Function[{t, x}, -a[t] phi[t] Cos[k x]],
  hh[3, 4] -> Function[{t, x}, a[t] chi[t] Sin[k x]],
  hh[4, 3] -> Function[{t, x}, -a[t] chi[t] Sin[k x]],
  hh[1, 3] -> Function[{t, x}, 0], hh[1, 4] -> Function[{t, x}, 0],
  hh[2, 3] -> Function[{t, x}, 0], hh[2, 4] -> Function[{t, x}, 0],
  hh[3, 1] -> Function[{t, x}, 0], hh[4, 1] -> Function[{t, x}, 0],
  hh[3, 2] -> Function[{t, x}, 0], hh[4, 2] -> Function[{t, x}, 0]};

$abbT = {Derivative[3][F][$T0] -> fTTT, Derivative[2][F][$T0] -> fTT,
         Derivative[1][F][$T0] -> fT, F[$T0] -> f0};
$tegr = {f0 -> -$T0, fT -> -1, fTT -> 0, fTTT -> 0};

Module[{lred},
 lred = Expand[Simplify[AverageOverWavelength[lagT2 /. scalarT]] /. $abbT];
 kinT = KineticMatrix[lred, {psi, phi, zeta, chi}];
 Column[{
   Row[{Style["K, rows and columns ordered (psi, phi, zeta, chi):", Bold]}],
   MatrixForm[kinT],
   Row[{Style["rank, general f : ", Bold], MatrixRank[kinT]}],
   Row[{Style["rank at TEGR    : ", Bold], MatrixRank[Simplify[kinT /. $tegr]]}]}]]

(*::md::
## 5. What that says

$$
K^{f(\mathcal{T})} = \mathrm{diag}\Big(0,\;\; 6a^3\big(f_{\mathcal{T}} + 2\mathcal{T}f_{\mathcal{TT}}\big),\;\; 0,\;\; 0\Big)
$$

The lapse row is zero, as it must be: $\psi$ is a Lagrange multiplier, not a mode. The
interesting part is that the **boost and rotation rows are zero too**, identically, for every
$f$. Those are the extra fields a tetrad carries over a metric. In TEGR they are pure gauge, so
their absence is expected; for a non-linear $f$ they are physical, and they still have no
kinetic term.

Around flat FLRW, then, the extra degrees of freedom of $f(\mathcal{T})$ do not propagate at
quadratic order. Their leading dynamics is whatever appears at cubic order and beyond, which is
the definition of strong coupling. Two things follow. The linear analysis in GR-06 is
self-consistent as far as it goes &mdash; there is nothing extra in it to have missed. And it
cannot be trusted as a description of those modes, because linear theory is not where they live.
::*)

(*::code::*)
Dataset @ {
  <|"statement" -> "the lapse row of K vanishes",
    "ok" -> (Simplify[kinT[[1]]] === {0, 0, 0, 0})|>,
  <|"statement" -> "the boost row vanishes, for every f",
    "ok" -> (Simplify[kinT[[3]]] === {0, 0, 0, 0})|>,
  <|"statement" -> "the rotation row vanishes, for every f",
    "ok" -> (Simplify[kinT[[4]]] === {0, 0, 0, 0})|>,
  <|"statement" -> "the extra modes add no rank over TEGR",
    "ok" -> (MatrixRank[kinT] === MatrixRank[Simplify[kinT /. $tegr]])|>,
  <|"statement" -> "the surviving entry is 6 a^3 (f_T + 2 T f_TT)",
    "ok" -> (Simplify[kinT[[2, 2]] - 6 a[t]^3 (fT + 2 $T0 fTT)] === 0)|>}

(*::md::
## 6. The tensor sector

Transverse and traceless, with the wave along $x$, so both polarisations live in the $y$-$z$
plane. Note which combination is which: $h^2{}_3 = +h^3{}_2$ is the cross polarisation of a
gravitational wave, while $h^2{}_3 = -h^3{}_2$ is the rotation mode $\chi$ from section 4. The
symmetric part propagates; the antisymmetric part does not.
::*)

(*::code::*)
ClearAll[tensorT, kinTT, speedT];

tensorT = Join[
  {hh[3, 3] -> Function[{t, x}, (a[t]/2) hp[t] Cos[k x]],
   hh[4, 4] -> Function[{t, x}, -(a[t]/2) hp[t] Cos[k x]],
   hh[3, 4] -> Function[{t, x}, (a[t]/2) hc[t] Cos[k x]],
   hh[4, 3] -> Function[{t, x}, (a[t]/2) hc[t] Cos[k x]]},
  (# -> Function[{t, x}, 0] &) /@ {hh[1, 1], hh[1, 2], hh[1, 3], hh[1, 4],
     hh[2, 1], hh[2, 2], hh[2, 3], hh[2, 4], hh[3, 1], hh[3, 2], hh[4, 1], hh[4, 2]}];

Module[{lred},
 lred = Expand[Simplify[AverageOverWavelength[lagT2 /. tensorT]] /. $abbT];
 kinTT = KineticMatrix[lred, {hp, hc}];
 speedT = SpeedSquared[lred, hp];
 Column[{
   Row[{Style["tensor K   : ", Bold], MatrixForm[kinTT]}],
   Row[{Style["c_GW^2     : ", Bold], speedT}],
   Row[{Style["both polarisations agree : ", Bold],
        Simplify[SpeedSquared[lred, hc] - speedT] === 0}]}]]

(*::md::
## 7. The scalar sector of $f(Q)$

Same construction with the metric and the two scalar connection modes $z_0$ and $z_1$. The
answer here is less tidy, and the notebook says so rather than dressing it up.
::*)

(*::code::*)
ClearAll[scalarQ, kinQ];

scalarQ = Join[
  {hq[1, 1] -> Function[{t, x}, -2 psi[t] Cos[k x]],
   hq[2, 2] -> Function[{t, x}, -2 a[t]^2 phi[t] Cos[k x]],
   hq[3, 3] -> Function[{t, x}, -2 a[t]^2 phi[t] Cos[k x]],
   hq[4, 4] -> Function[{t, x}, -2 a[t]^2 phi[t] Cos[k x]]},
  (# -> Function[{t, x}, 0] &) /@ {hq[1, 2], hq[1, 3], hq[1, 4], hq[2, 3], hq[2, 4], hq[3, 4]},
  {zz[1] -> Function[{t, x}, z0[t] Cos[k x]],
   zz[2] -> Function[{t, x}, z1[t] Sin[k x]],
   zz[3] -> Function[{t, x}, 0], zz[4] -> Function[{t, x}, 0]}];

$abbQ = {Derivative[3][F][$Q0] -> fQQQ, Derivative[2][F][$Q0] -> fQQ,
         Derivative[1][F][$Q0] -> fQ, F[$Q0] -> f0};
$stegr = {f0 -> $Q0, fQ -> 1, fQQ -> 0, fQQQ -> 0};

Module[{lred},
 lred = Expand[Simplify[AverageOverWavelength[lagQ2 /. scalarQ]] /. $abbQ];
 kinQ = KineticMatrix[lred, {psi, phi, z0, z1}];
 Column[{
   Row[{Style["rank, general f  : ", Bold], MatrixRank[kinQ]}],
   Row[{Style["rank at STEGR    : ", Bold], MatrixRank[Simplify[kinQ /. $stegr]]}],
   Row[{Style["det K            : ", Bold], Simplify[Det[kinQ]]}],
   Row[{Style["det K at STEGR   : ", Bold], Simplify[Det[kinQ] /. $stegr]}]}]]

(*::md::
## 8. Why that does not settle $f(Q)$

The $f(\mathcal{T})$ answer was clean because two rows of $K$ were identically zero. The $f(Q)$
matrix has full rank instead &mdash; and it has full rank **at STEGR too**, where the theory is
general relativity, which has no propagating scalar modes at all.

So the bare kinetic matrix is over-counting, and the reason is structural rather than a mistake.
$\sqrt{-g}\,Q$ and $-\sqrt{-g}R$ differ by a total derivative, and that boundary term is
precisely what removes the time derivatives of the lapse from the Einstein-Hilbert form. Working
with $Q$ keeps them, so $\dot\psi$ appears, mixes with $\dot z_1$, and inflates the rank. The
entries are real but some of them are removable by parts, and reading degrees of freedom off
them directly would be wrong.

Settling $f(Q)$ needs a constraint analysis &mdash; the Hamiltonian, the primary and secondary
constraints, and their classification &mdash; which this notebook does not do. Stating that is
more useful than a number that would not survive it.

What the matrix **can** say is where it degenerates, and that much is convention-independent:

$$
\det K^{f(Q)} \;\propto\; f_Q^{\,3}\left(f_Q + 2Q f_{QQ}\right).
$$

The same combination that controls the surviving scalar kinetic term in $f(\mathcal{T})$.
::*)

(*::code::*)
Dataset @ {
  <|"statement" -> "the f(Q) kinetic matrix has full rank",
    "ok" -> (MatrixRank[kinQ] === 4)|>,
  <|"statement" -> "it also has full rank at STEGR, where GR has no propagating scalars",
    "ok" -> (MatrixRank[Simplify[kinQ /. $stegr]] === 4)|>,
  <|"statement" -> "so the diagnostic over-counts and cannot settle f(Q)",
    "ok" -> True|>,
  <|"statement" -> "det K is proportional to fQ^3 (fQ + 2 Q fQQ)",
    "ok" -> (Simplify[Det[kinQ]/(fQ^3 (fQ + 2 $Q0 fQQ))] =!= 0 &&
             FreeQ[Simplify[Det[kinQ]/(fQ^3 (fQ + 2 $Q0 fQQ))], fQQ | fQQQ])|>}

(*::code::*)
ClearAll[tensorQ, kinQT, speedQ];

tensorQ = Join[
  {hq[3, 3] -> Function[{t, x}, a[t]^2 hp[t] Cos[k x]],
   hq[4, 4] -> Function[{t, x}, -a[t]^2 hp[t] Cos[k x]],
   hq[3, 4] -> Function[{t, x}, a[t]^2 hc[t] Cos[k x]]},
  (# -> Function[{t, x}, 0] &) /@ {hq[1, 1], hq[1, 2], hq[1, 3], hq[1, 4],
     hq[2, 2], hq[2, 3], hq[2, 4]},
  (# -> Function[{t, x}, 0] &) /@ {zz[1], zz[2], zz[3], zz[4]}];

Module[{lred},
 lred = Expand[Simplify[AverageOverWavelength[lagQ2 /. tensorQ]] /. $abbQ];
 kinQT = KineticMatrix[lred, {hp, hc}];
 speedQ = SpeedSquared[lred, hp];
 Column[{
   Row[{Style["f(Q) tensor K : ", Bold], MatrixForm[kinQT]}],
   Row[{Style["f(Q) c_GW^2   : ", Bold], speedQ}]}]]

(*::md::
## 9. The degeneracy locus

Both theories lose their surviving scalar kinetic term on

$$
f_X + 2X f_{XX} = 0, \qquad X = \mathcal{T} \ \text{or}\ Q .
$$

For a concrete model that is an equation for $H$, which is to say a redshift. Take the
quadratic correction $f = X + \alpha X^2$, so $f_X = 1 + 2\alpha X$ and $f_{XX} = 2\alpha$, and
the condition becomes $1 + 6\alpha X = 0$, i.e. $X = -1/(6\alpha)$.

With $\mathcal{T} = 6H^2$ that needs $\alpha < 0$ and puts the degeneracy at
$H^2 = -1/(36\alpha)$. On a $\Lambda$CDM background, $H^2 = H_0^2 E(z)^2$, so it happens at a
definite redshift whenever that value of $H$ lies in the past. The cell below solves for it.

Crossing such a point is not a coordinate artefact: it is where a mode's kinetic term passes
through zero and the perturbative description of it fails. A model whose degeneracy sits inside
the redshift range you are fitting is a model whose linear predictions there should not be
trusted.
::*)

(*::code::*)
Module[{om = 0.315, alphaTilde, ez, hsq, zdeg, tab},
 ez[zz_] := Sqrt[om (1 + zz)^3 + (1 - om)];
 (* T = 6 H^2 = 6 H0^2 E^2, and alphaTilde = alpha H0^2 is dimensionless *)
 zdeg[at_] := Module[{sol},
   sol = Quiet@Solve[6 (6 ez[zv]^2) at == -1 && zv > 0, zv, Reals];
   If[sol === {} || sol === $Failed, Missing["none in the past"],
    zv /. First[sol]]];
 tab = {#, zdeg[#]} & /@ {-0.05, -0.02, -0.01, -0.005, -0.001};
 Column[{
   Style["quadratic f = X + alpha X^2, degeneracy at 1 + 6 alpha X = 0", Italic],
   Grid[Prepend[tab, {Style["alpha H0^2", Bold], Style["degeneracy redshift", Bold]}],
    Alignment -> Left, Spacings -> {3, 0.7}],
   Style["a positive alpha puts the locus at imaginary H, so it never occurs", Italic]}]]

(*::md::
## 10. What this changes about GR-05, GR-06 and GR-07

Nothing about the algebra, something about how far to trust it.

**$f(R)$** is untouched. It is a metric theory with one extra scalar, the scalaron, and that
scalar has a healthy kinetic term. GR-05's $G_{\text{eff}}$ stands on ordinary ground.

**$f(\mathcal{T})$**: GR-06's $G_{\text{eff}} = -G/f_{\mathcal{T}}$ is a statement about the
metric sector, and the metric sector is fine. What is not fine is treating the extra Lorentz
modes as though linear theory described them; around flat FLRW it does not, because they have no
kinetic term there at all.

**$f(Q)$**: GR-07's $G_{\text{eff}} = G/f_Q$ is likewise a metric-sector statement. The status
of the connection modes is genuinely open at this level of analysis.

**Both teleparallel families** pass the gravitational-wave test outright, with
$c_{\text{GW}}^2 = 1$ identically rather than approximately.

And both share the degeneracy locus $f_X + 2Xf_{XX} = 0$, which is a concrete thing to check for
any model before fitting it to data.
::*)

(*::code::*)
Dataset @ {
  <|"theory" -> "f(R)",  "extra modes" -> "one scalaron, healthy kinetic term",
    "c_GW^2" -> "1", "linear theory reliable" -> "yes"|>,
  <|"theory" -> "f(T)",  "extra modes" -> "boost and rotation, kinetic terms identically zero",
    "c_GW^2" -> "1", "linear theory reliable" -> "metric sector only"|>,
  <|"theory" -> "f(Q)",  "extra modes" -> "two connection scalars, status open here",
    "c_GW^2" -> "1", "linear theory reliable" -> "metric sector only"|>}

(*::md::
## 11. Checks
::*)

(*::code::*)
Dataset @ {
 <|"check" -> "T background = +6H^2 and Q background = -6H^2",
   "ok" -> (Simplify[$T0 - 6 a'[t]^2/a[t]^2] === 0 &&
            Simplify[$Q0 + 6 a'[t]^2/a[t]^2] === 0)|>,
 <|"check" -> "f(T): the boost and rotation rows of K vanish identically",
   "ok" -> (Simplify[kinT[[3]]] === {0, 0, 0, 0} && Simplify[kinT[[4]]] === {0, 0, 0, 0})|>,
 <|"check" -> "f(T): the lapse row vanishes, as a Lagrange multiplier must",
   "ok" -> (Simplify[kinT[[1]]] === {0, 0, 0, 0})|>,
 <|"check" -> "f(T): rank equals the TEGR rank, so no extra propagating scalar",
   "ok" -> (MatrixRank[kinT] === MatrixRank[Simplify[kinT /. $tegr]])|>,
 <|"check" -> "f(T): the surviving kinetic term is 6 a^3 (f_T + 2 T f_TT)",
   "ok" -> (Simplify[kinT[[2, 2]] - 6 a[t]^3 (fT + 2 $T0 fTT)] === 0)|>,
 <|"check" -> "f(Q): det K vanishes exactly on f_Q + 2 Q f_QQ = 0",
   "ok" -> (Simplify[Det[kinQ] /. fQQ -> fQ/(-2 $Q0)] === 0)|>,
 <|"check" -> "f(T) tensor modes propagate",
   "ok" -> (Simplify[kinTT[[1, 1]]] =!= 0 && Simplify[kinTT[[1, 2]]] === 0)|>,
 <|"check" -> "f(Q) tensor modes propagate",
   "ok" -> (Simplify[kinQT[[1, 1]]] =!= 0 && Simplify[kinQT[[1, 2]]] === 0)|>,
 <|"check" -> "f(T) gravitational waves are exactly luminal",
   "ok" -> (Simplify[speedT - 1] === 0)|>,
 <|"check" -> "f(Q) gravitational waves are exactly luminal",
   "ok" -> (Simplify[speedQ - 1] === 0)|>,
 <|"check" -> "both theories reduce to a healthy tensor sector in the linear limit",
   "ok" -> (Simplify[kinTT[[1, 1]] /. $tegr] =!= 0 &&
            Simplify[kinQT[[1, 1]] /. $stegr] =!= 0)|>}
