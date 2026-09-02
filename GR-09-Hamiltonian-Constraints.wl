(* ::Package:: *)

(* ============================================================================
   GR-09 : Hamiltonian Constraint Analysis and the Degrees of Freedom of f(Q)
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-09-Hamiltonian-Constraints.wl

   ============================================================================ *)

(*::md::
# Hamiltonian Constraints, and How Many Modes $f(Q)$ Really Has

GR-08 ended with an admission. The kinetic matrix said the $f(Q)$ scalar sector had rank four
where STEGR had rank three, which would mean an extra propagating mode &mdash; but that matrix is
built from $\partial^2 L/\partial\dot q^2$, and a Lagrangian carrying **second** time derivatives
does not have its dynamics in that object at all. The honest conclusion there was that the
question needed an Ostrogradsky reduction followed by a Dirac constraint analysis, and that
until someone did it the number was not known.

This notebook does it. The answer is that general $f(Q)$ carries **one propagating scalar degree
of freedom**, where general relativity carries none, and the mechanism is visible in the
constraint algebra: the number of constraints does not change, but two of them stop commuting.

$$
\text{GR: } 4 \text{ first class} \qquad\longrightarrow\qquad
f(Q): \; 6 \text{ first class},\; 6 \text{ second class}
$$

Two constraints that generate gauge transformations in general relativity become second class
once $f_{QQ}\neq 0$. A broken gauge symmetry leaves half of its pair behind as a physical mode.

## The counting rule

For a system with $n$ coordinates, $F$ first-class constraints and $S$ second-class ones,

$$
\text{DOF} \;=\; \frac{2n - 2F - S}{2}.
$$

First-class constraints cost twice: once for the constraint, once for the gauge direction it
generates. Second-class constraints cost once each.

## Two traps that make the obvious implementation wrong

Both of these produced confident, wrong numbers before the acceptance test caught them, and both
are worth stating plainly because neither is obvious.

**Primary constraints come from the null space of the velocity Hessian, not from momenta that
happen to be velocity free.** It is tempting to scan the momenta $p_i = \partial L/\partial\dot
q^i$ and call $q^i$ constrained when $p_i$ contains no velocity. A single cross term
$\dot q^1\dot q^2$ defeats that: every momentum then depends on a velocity, while the Hessian is
still degenerate and the constraint is still there. The constraints are the null combinations
$v_i\,(p_i - \partial L/\partial\dot q^i)$.

**A constraint is new only when it is linearly independent of the ones already found.** The
Dirac chain regenerates old constraints multiplied by background factors like $a(t)$. Testing new
against old with an equality test never recognises them, the chain never terminates, and the
count runs off to nonsense. Since the quadratic action makes every constraint a linear form in
$(q,p)$, the right test is a rank test on their coefficient vectors.
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
| **GR-09** | a quadratic action in, primary and secondary constraints and a mode count out |
::*)

(*::md::
## 1. The constraint counter

`DiracCount` takes a Lagrangian quadratic in the velocities, the list of coordinates, and runs
the Dirac algorithm: velocity Hessian, primary constraints from its null space, canonical
Hamiltonian, consistency conditions, secondary constraints, then the classification into first
and second class from the rank of the matrix of Poisson brackets.

The consistency conditions are linear in the undetermined multipliers, $A\,u + b = 0$. Where $A$
can be inverted the condition fixes a multiplier and produces nothing new; the genuine secondary
constraints are the components of $b$ along the **left null space** of $A$, which no choice of
multiplier can absorb.
::*)

(*::code::*)
ClearAll[DiracCount, OstrogradskyReduce, ReduceByParts];

DiracCount[lagIn_, qsl_List, opts : OptionsPattern[]] :=
 Module[{n, toQ, vel, lag, hmat, cvec, rk, nulls, chi, ind, sub, vSol, hc,
   basis, covec, mat, addIf, pb, us, hT, evolve, allCon, newCon, stage,
   cmat, rk2, fc, dof, log = {}},
  n = Length[qsl];
  toQ = Table[qsl[[i]][t] -> qq[i], {i, n}];
  vel = Table[Derivative[1][qsl[[i]]][t] -> vv[i], {i, n}];
  lag = Expand[lagIn /. vel /. toQ];
  If[! FreeQ[lag, Derivative[_][_][t]],
   Return[<|"ok" -> False, "why" -> "the Lagrangian is not first order"|>]];

  (* primary constraints from the null space of the velocity Hessian *)
  hmat = Simplify[Table[D[lag, vv[i], vv[j]], {i, n}, {j, n}]];
  If[! FreeQ[hmat, vv[_]],
   Return[<|"ok" -> False, "why" -> "the Lagrangian is not quadratic in velocities"|>]];
  cvec = Simplify[Table[D[lag, vv[i]] /. Table[vv[m] -> 0, {m, n}], {i, n}]];
  rk = MatrixRank[hmat];
  nulls = Simplify[NullSpace[hmat]];
  chi = Simplify[# . (Table[pp[i], {i, n}] - cvec)] & /@ nulls;

  (* invert the Hessian on a maximal non-degenerate block; the null directions
     are undetermined velocities and drop out of the canonical Hamiltonian *)
  ind = If[rk === 0, {},
    SelectFirst[Subsets[Range[n], {rk}], Simplify[Det[hmat[[#, #]]]] =!= 0 &]];
  sub = If[rk === 0, {},
    Quiet@First@Solve[Thread[(Table[pp[i], {i, n}] - cvec)[[ind]]
        == hmat[[ind, ind]] . Table[vv[i], {i, ind}]], Table[vv[i], {i, ind}]]];
  vSol = Join[If[rk === 0, {}, sub], Table[vv[i] -> 0, {i, Complement[Range[n], ind]}]];
  hc = Simplify[(Sum[pp[i] vv[i], {i, n}] - lag) /. vSol];
  If[! FreeQ[hc, vv[_]],
   Return[<|"ok" -> False, "why" -> "the canonical Hamiltonian is not velocity free"|>]];

  (* every constraint is a linear form in (q, p): keep only independent ones *)
  basis = Join[Table[qq[i], {i, n}], Table[pp[i], {i, n}]];
  covec[c_] := Simplify[Join[Table[D[c, v], {v, basis}],
     {c /. Table[v -> 0, {v, basis}]}]];
  mat = {};
  addIf[c_] := Module[{v = covec[c], m2},
    If[Simplify[c] === 0, Return[False]];
    m2 = Append[mat, v];
    If[mat === {},
     If[MatrixRank[m2] > 0, mat = m2; True, False],
     If[MatrixRank[m2] > MatrixRank[mat], mat = m2; True, False]]];

  allCon = Select[chi, addIf];
  chi = allCon;
  pb[f_, g_] := Simplify[Sum[D[f, qq[i]] D[g, pp[i]] - D[f, pp[i]] D[g, qq[i]], {i, n}]];
  us = Table[uu[i], {i, Length[chi]}];
  hT = hc + Sum[us[[i]] chi[[i]], {i, Length[chi]}];
  evolve[f_] := Simplify[pb[f, hT] + D[f, t]];

  newCon = allCon; stage = 0;
  While[stage < 8 && newCon =!= {},
   stage++;
   Module[{cons, amat, bvec, ln, sec},
    cons = evolve /@ newCon;
    amat = Simplify[Table[D[c, uu[j]], {c, cons}, {j, Length[us]}]];
    bvec = Simplify[cons /. Table[uu[j] -> 0, {j, Length[us]}]];
    ln = NullSpace[Transpose[amat]];
    sec = If[ln === {}, {}, DeleteCases[Simplify[ln . bvec], 0]];
    sec = Select[sec, addIf];
    AppendTo[log, "stage " <> ToString[stage] <> ": " <> ToString[Length[sec]] <> " new"];
    allCon = Join[allCon, sec]; newCon = sec]];
  If[newCon =!= {},
   Return[<|"ok" -> False, "why" -> "the constraint chain did not terminate"|>]];

  cmat = Simplify[Outer[pb, allCon, allCon, 1]];
  rk2 = MatrixRank[cmat];
  fc = Length[allCon] - rk2;
  dof = (2 n - 2 fc - rk2)/2;
  <|"ok" -> True, "Coordinates" -> n, "HessianRank" -> rk, "Primary" -> Length[chi],
    "Total" -> Length[allCon], "SecondClass" -> rk2, "FirstClass" -> fc,
    "DOF" -> dof, "Chain" -> log, "Constraints" -> allCon|>];

(* trade zdot for a new coordinate and impose the relation with a multiplier *)
OstrogradskyReduce[lag_, zs_List] :=
 Module[{l = lag, extra = {}, w, lm},
  Do[w = Symbol["w$" <> ToString[z]]; lm = Symbol["lam$" <> ToString[z]];
   l = Expand[l /. {Derivative[2][z][t] -> Derivative[1][w][t],
                    Derivative[1][z][t] -> w[t]}] + lm[t] (Derivative[1][z][t] - w[t]);
   extra = Join[extra, {w, lm}], {z, zs}];
  {Expand[l], extra}];

(*::md::
## 2. Baseline: general relativity must give zero

Nothing below is trustworthy until the machinery reproduces the one case where the answer is
known. Vacuum general relativity with a cosmological constant, on de Sitter, in the scalar
sector, has **no propagating scalar mode**.

The scalar sector is kept complete &mdash; lapse $\psi$, curvature $\phi$, shift $b$ and
anisotropy $e$ &mdash; and no gauge is fixed:

$$
N = 1 + \epsilon\,\psi\cos kx,\qquad N_x = \epsilon\, a\, b \sin kx,
$$
$$
h_{ij} = a^2\big[(1-2\epsilon\phi\cos kx)\,\delta_{ij} - 2\epsilon\, e \cos kx\,\delta_{ix}\delta_{jx}\big].
$$

This matters more than it looks. **Fixing Newtonian gauge in the action before varying deletes
the momentum constraint**, because that constraint is what varying the shift produces, and
Newtonian gauge sets the shift to zero. Do the count that way and general relativity comes out
with one propagating scalar instead of none &mdash; off by exactly one, because exactly one
constraint was thrown away. That is a legitimate thing to do for $G_{\text{eff}}$, where the
missing equation is supplied from the covariant field equations, and a fatal thing to do here.

Writing the action in ADM form matters too. In $\sqrt{h}\,N\,(K_{ij}K^{ij} - K^2 + {}^{(3)}R)$
the lapse and shift carry no time derivatives at all, so their momenta are primary constraints by
construction. Reaching a first-order Lagrangian instead by integrating $\ddot\phi$ away by parts
generates a $\dot\psi\dot\phi$ cross term, hands the lapse a velocity, and destroys the very
structure the count depends on.
::*)

(*::code::*)
ClearAll[psi, phi, b, ee, z0, z1, hh, k, eps];
$Assumptions = hh > 0 && k > 0;
ser[e_] := Normal[Series[e, {eps, 0, 2}]];
avg[e_] := Simplify[(k/(2 Pi)) Integrate[e, {x, 0, 2 Pi/k}]];
aa[t_] := Exp[hh t];
cc := Cos[k x]; ss := Sin[k x];
ys = {x, y, z};

nlapse = 1 + eps psi[t] cc;
nshift = {eps aa[t] b[t] ss, 0, 0};
h3 = DiagonalMatrix[{aa[t]^2 (1 - 2 eps phi[t] cc - 2 eps ee[t] cc),
                     aa[t]^2 (1 - 2 eps phi[t] cc), aa[t]^2 (1 - 2 eps phi[t] cc)}];
h3inv = ser[Inverse[h3]];

ch3 = Module[{p, m, n, s}, Table[ser[(1/2) Sum[h3inv[[p, s]] (D[h3[[s, m]], ys[[n]]]
    + D[h3[[s, n]], ys[[m]]] - D[h3[[m, n]], ys[[s]]]), {s, 3}]], {p, 3}, {m, 3}, {n, 3}]];
rm3 = Module[{p, m, n, q, s}, Table[ser[D[ch3[[p, m, q]], ys[[n]]]
    - D[ch3[[p, m, n]], ys[[q]]] + Sum[ch3[[p, n, s]] ch3[[s, m, q]]
    - ch3[[p, q, s]] ch3[[s, m, n]], {s, 3}]], {p, 3}, {m, 3}, {n, 3}, {q, 3}]];
rc3 = Module[{m, n, p}, Table[ser[Sum[rm3[[p, m, p, n]], {p, 3}]], {m, 3}, {n, 3}]];
r3 = Module[{m, n}, ser[Sum[h3inv[[m, n]] rc3[[m, n]], {m, 3}, {n, 3}]]];

dnab = Module[{m, n, s}, Table[ser[D[nshift[[n]], ys[[m]]]
    - Sum[ch3[[s, m, n]] nshift[[s]], {s, 3}]], {m, 3}, {n, 3}]];
kext = Module[{m, n}, Table[ser[(1/(2 nlapse)) (D[h3[[m, n]], t]
    - dnab[[m, n]] - dnab[[n, m]])], {m, 3}, {n, 3}]];
kmix = Module[{m, n, s}, Table[ser[Sum[h3inv[[m, s]] kext[[s, n]], {s, 3}]], {m, 3}, {n, 3}]];
ktr = Module[{m}, ser[Sum[kmix[[m, m]], {m, 3}]]];
ksq = Module[{m, n}, ser[Sum[kmix[[m, n]] kmix[[n, m]], {m, 3}, {n, 3}]]];

lagADM = ser[Sqrt[Det[h3]] nlapse (ksq - ktr^2 + r3 - 6 hh^2)];
l2GR = Simplify[avg[Coefficient[lagADM, eps, 2]]];
grResult = DiracCount[l2GR, {psi, phi, b, ee}];
grResult

(*::md::
The lapse and the shift give two primary constraints, the chain adds the Hamiltonian and
momentum constraints, all four are first class, and

$$
\text{DOF} = \frac{2\cdot 4 - 2\cdot 4 - 0}{2} = 0.
$$

Textbook general relativity. The machinery can now be pointed at something whose answer is not
known in advance.
::*)

(*::md::
## 3. The complete scalar sector of $f(Q)$

The metric is the same one. The connection is the general flat, torsion-free connection built
from Stückelberg fields, as in GR-07, expanded about the coincident gauge:

$$
\xi^0 = t + \epsilon\, z_0 \cos kx,\qquad \xi^x = x + \epsilon\, z_1 \sin kx,
$$
$$
\Gamma^\lambda{}_{\mu\nu} = \frac{\partial x^\lambda}{\partial\xi^\rho}\,
\partial_\mu\partial_\nu \xi^\rho .
$$

Because $\Gamma$ carries $\partial\partial\xi$, the non-metricity scalar carries **second time
derivatives** of $z_0$ and $z_1$. That is exactly what made GR-08's kinetic matrix the wrong
object, and it is handled here by the Ostrogradsky reduction: introduce $w = \dot z$ as an
independent coordinate and impose $\dot z = w$ with a multiplier, enlarging the phase space
honestly instead of pretending the second derivatives are not there.
::*)

(*::code::*)
xs = {t, x, y, z};
gfull = ser@ArrayFlatten[{{{{-(nlapse^2) + nshift . h3inv . nshift}}, {nshift}},
                          {Transpose[{nshift}], h3}}];
ginv = ser[Inverse[gfull]];
sqg = ser[Sqrt[-Det[gfull]]];

xi = {t + eps z0[t] cc, x + eps z1[t] ss, y, z};
jac = Table[D[xi[[p]], xs[[m]]], {p, 4}, {m, 4}];
jinv = ser[Inverse[jac]];
gam = Module[{l, m, i, p}, Table[ser[Sum[jinv[[l, p]] D[xi[[p]], xs[[m]], xs[[i]]], {p, 4}]],
   {l, 4}, {m, 4}, {i, 4}]];

qd = Module[{l, m, i, p}, Table[ser[D[gfull[[m, i]], xs[[l]]]
    - Sum[gam[[p, l, m]] gfull[[p, i]] + gam[[p, l, i]] gfull[[m, p]], {p, 4}]],
   {l, 4}, {m, 4}, {i, 4}]];
qu = Module[{l, m, i, p, q, al}, Table[ser[Sum[
     ginv[[l, p]] ginv[[m, q]] ginv[[i, al]] qd[[p, q, al]], {p, 4}, {q, 4}, {al, 4}]],
   {l, 4}, {m, 4}, {i, 4}]];
qtr = Module[{l, m, i}, Table[ser[Sum[ginv[[m, i]] qd[[l, m, i]], {m, 4}, {i, 4}]], {l, 4}]];
qtl = Module[{l, m, i}, Table[ser[Sum[ginv[[m, i]] qd[[m, l, i]], {m, 4}, {i, 4}]], {l, 4}]];
qup = Module[{l, p}, Table[ser[Sum[ginv[[l, p]] qtr[[p]], {p, 4}]], {l, 4}]];
qtu = Module[{l, p}, Table[ser[Sum[ginv[[l, p]] qtl[[p]], {p, 4}]], {l, 4}]];
qscalar = Module[{l, m, i}, ser[
   -(1/4) Sum[qd[[l, m, i]] qu[[l, m, i]], {l, 4}, {m, 4}, {i, 4}]
   + (1/2) Sum[qd[[l, m, i]] qu[[i, m, l]], {l, 4}, {m, 4}, {i, 4}]
   + (1/4) Sum[qtr[[l]] qup[[l]], {l, 4}] - (1/2) Sum[qtr[[l]] qtu[[l]], {l, 4}]]];

$Q0 = Simplify[qscalar /. eps -> 0];
$Q0

(*::md::
$Q_0 = -6H^2$, the same convention as GR-02, GR-04 and GR-07.

## 4. The background equations, read off before averaging

The de Sitter background is only a solution of $f(Q)$ for an $f$ satisfying one condition, and
that condition has to be imposed or the quadratic action is built around a point that is not a
solution. Leave it out and the linear branch $f''=0$ &mdash; which *is* general relativity
&mdash; comes out with a propagating mode.

The condition has to be read off the linear term **before** the wavelength average. Every linear
term carries a factor $\cos kx$ or $\sin kx$, whose mean over a wavelength vanishes identically,
so the averaged linear term is zero whether the background is on shell or not.
::*)

(*::code::*)
abbrev = {Derivative[3][F][$Q0] -> fQQQ, Derivative[2][F][$Q0] -> fQQ,
          Derivative[1][F][$Q0] -> fQ, F[$Q0] -> f0};
lagQ = ser[sqg F[qscalar]] /. abbrev;
elOp[l_, q_] := Simplify[D[l, q[t]] - D[D[l, Derivative[1][q][t]], t]
    + D[D[l, Derivative[2][q][t]], {t, 2}]];
linTerm = Coefficient[lagQ, eps, 1];
bgEnergy = Simplify[elOp[linTerm, psi] /. x -> 0];
bgEvolve = Simplify[elOp[linTerm, phi] /. x -> 0];
bgRule = First@Solve[bgEnergy == 0, f0];
{bgEnergy == 0, bgRule, Simplify[bgEvolve /. bgRule] === 0}

(*::md::
$$
f_0 = -12 H^2 f_Q \;=\; 2\,Q_0\, f_Q ,
$$

one condition, and the second background equation follows from it. STEGR with a cosmological
constant, $f = Q - 2\Lambda$ with $\Lambda = 3H^2$, satisfies it: $f_0 = -6H^2 - 6H^2 = -12H^2$
and $f_Q = 1$.

## 5. The acceptance test, and then the answer

Three cases go through the same pipeline. Two of them are general relativity in disguise and must
return zero.
::*)

(*::code::*)
runCase[rules_] := Module[{l2, need, lred, extra},
  l2 = Simplify[avg[Coefficient[lagQ, eps, 2]] /. bgRule /. rules];
  need = Select[{z0, z1}, ! FreeQ[l2, Derivative[2][#][t]] &];
  {lred, extra} = If[need === {}, {l2, {}}, OstrogradskyReduce[l2, need]];
  DiracCount[lred, Join[{psi, phi, b, ee, z0, z1}, extra]]];

stegrResult = runCase[{fQ -> 1, fQQ -> 0, fQQQ -> 0}];
linearResult = runCase[{fQQ -> 0, fQQQ -> 0}];
fqResult = runCase[{}];

Dataset @ {
 <|"case" -> "general relativity, metric ADM sector", "n" -> grResult["Coordinates"],
   "constraints" -> grResult["Total"], "first class" -> grResult["FirstClass"],
   "second class" -> grResult["SecondClass"], "DOF" -> grResult["DOF"]|>,
 <|"case" -> "STEGR, f = Q - 2 Lambda", "n" -> stegrResult["Coordinates"],
   "constraints" -> stegrResult["Total"], "first class" -> stegrResult["FirstClass"],
   "second class" -> stegrResult["SecondClass"], "DOF" -> stegrResult["DOF"]|>,
 <|"case" -> "linear branch, f'' = 0", "n" -> linearResult["Coordinates"],
   "constraints" -> linearResult["Total"], "first class" -> linearResult["FirstClass"],
   "second class" -> linearResult["SecondClass"], "DOF" -> linearResult["DOF"]|>,
 <|"case" -> "general f(Q)", "n" -> fqResult["Coordinates"],
   "constraints" -> fqResult["Total"], "first class" -> fqResult["FirstClass"],
   "second class" -> fqResult["SecondClass"], "DOF" -> fqResult["DOF"]|>}

(*::md::
## 6. What the table says

The constraint **count** is identical in the linear branch and in general $f(Q)$: twelve, in both
cases, eight primary and four secondary. Nothing is created or destroyed. What changes is the
algebra:

| | first class | second class | DOF |
|---|---|---|---|
| $f'' = 0$ (general relativity) | 8 | 4 | **0** |
| general $f(Q)$ | 6 | 6 | **1** |

Two constraints move from first class to second. In the counting rule a first-class constraint
removes two phase-space dimensions and a second-class one removes one, so moving two constraints
across the line frees up $2\times 2 - 2\times 1 = 2$ dimensions, which is one degree of freedom.

Physically: general relativity's scalar sector is pure gauge, and $f_{QQ}\neq 0$ breaks one of
those gauge symmetries. The mode that was gauge in general relativity becomes physical. This is
the same story as $f(R)$, where the extra scalar is the scalaron, but the mechanism here runs
through the connection rather than through a conformal factor.

Note also that this is a **discontinuous** limit in the same sense as the strong-coupling
statement of GR-08: the extra mode does not smoothly go away as $f_{QQ}\to 0$, it disappears
abruptly when the constraint becomes first class again. A theory that is arbitrarily close to
general relativity in its background and linear response still has an extra mode in its
constraint structure.
::*)

(*::md::
## 7. What this does and does not settle

**Settled.** On a de Sitter background, in the scalar sector, for one Fourier mode, general
$f(Q)$ propagates one scalar mode and general relativity propagates none. The count is
reproduced with concrete numbers as well as symbolically, so it does not rest on a symbolic rank
being evaluated correctly.

**Not settled.** The background here is de Sitter with constant $H$, chosen because it is an
exact solution and keeps the algebra finite. Constraint structure can be background dependent
&mdash; that is precisely the strong-coupling worry of GR-08 &mdash; and a background with
$\dot H \neq 0$ could in principle change the rank of the bracket matrix. The vector and tensor
sectors are not analysed here; GR-08 covers the tensor sector and finds it healthy, with exactly
luminal waves. And the full non-perturbative count, with no expansion in $\epsilon$ at all, is a
larger problem than this notebook attempts: what is counted here is the degrees of freedom of the
linearised theory about this background.

**Corrected.** GR-08 §7 and §8 said the question could not be settled with the tools in that
notebook and that the rank-four kinetic matrix over-counted for reasons that were only half
understood. The first part was right. The second is now complete: the kinetic matrix is not
merely the wrong normalisation, it is the wrong object, because with $\ddot z$ in the Lagrangian
the phase space is larger than the coordinates suggest and $\dot z$ is an independent variable.
Rank four was a symptom of that, not a mode count.
::*)

(*::md::
## 8. $f(R)$, and a second route to the whole calculation

Everything above reaches a first-order Lagrangian by writing the action in ADM form, where the
lapse and shift carry no velocities. $f(R)$ cannot be reached that way, and the reason is worth
being precise about, because getting it wrong produces a plausible number.

In ADM variables the four-dimensional Ricci scalar is

$$
{}^{(4)}R \;=\; K_{ij}K^{ij} - K^2 + {}^{(3)}R \;+\; \text{total derivatives}.
$$

In general relativity those total derivatives multiply a **constant** and integrate away, which
is why $\sqrt{h}\,N(K_{ij}K^{ij} - K^2 + {}^{(3)}R)$ is a legitimate starting point. In $f(R)$
they multiply $f'(R)$, which is not constant, and dropping them changes the theory. It is the
same trap the ordinary minisuperspace reduction hits — see the companion library's
`theory/curvature.py`, which avoids it with a Lagrange multiplier for exactly this reason.

So this section does not integrate anything by parts by hand. It builds the full four-metric
back out of the ADM pieces already defined above,

$$
g_{00} = -(N^2 - N_iN^i), \qquad g_{0i} = N_i, \qquad g_{ij} = h_{ij},
$$

computes ${}^{(4)}R$ from it directly, and hands the result — second time derivatives and all —
to `OstrogradskyReduce`. Nothing is discarded, so nothing can be discarded wrongly.

**The scalar-tensor form.** With a Lagrange multiplier, $f(R)$ is

$$
L = \sqrt{-g}\,\big[\,F\,{}^{(4)}R - V(F)\,\big], \qquad F = f'(\chi), \qquad V = F\chi - f(\chi),
$$

so that $V'(F) = \chi$ and $V''(F) = 1/f''(\chi)$. The background must solve its own equations or
the constraint structure is not the theory's: on de Sitter with $F = f_0$ constant, varying $F$
gives $V'(f_0) = {}^{(4)}R = 12H^2$ and the Friedmann constraint gives $V(f_0) = 6f_0H^2$. Both
are imposed. $V''$ is left free, and it is the whole story:

- $V''$ finite and non-zero is a genuine $f(R)$, and the scalaron propagates;
- $f'' \to 0$ — general relativity — is $V'' \to \infty$, an infinitely heavy scalar, which does
  not propagate at all.

The control below runs the same pipeline with the scalaron switched off, so the covariant route
is checked against the answer §2 already established by the ADM route.
::*)

(*::code::*)
ClearAll[ff, ffld, vv2, f0fr, g4, g4inv, ch4, ric4, rr4, sq4, nsUp, nsSq];

(* nshift is N_i with a lower index -- see how dnab contracts it above *)
nsUp = Table[Sum[h3inv[[i, j]] nshift[[j]], {j, 3}], {i, 3}];
nsSq = ser[Sum[nshift[[i]] nsUp[[i]], {i, 3}]];

g4 = Table[0, {4}, {4}];
g4[[1, 1]] = ser[-(nlapse^2 - nsSq)];
Do[g4[[1, i + 1]] = nshift[[i]]; g4[[i + 1, 1]] = nshift[[i]], {i, 3}];
Do[g4[[i + 1, j + 1]] = h3[[i, j]], {i, 3}, {j, 3}];
g4 = ser[g4];

co4 = {t, x, y, z};
g4inv = ser[Inverse[g4]];

ch4 = Table[ser[(1/2) Sum[g4inv[[p, s]] (D[g4[[s, m]], co4[[n]]]
     + D[g4[[s, n]], co4[[m]]] - D[g4[[m, n]], co4[[s]]]), {s, 4}]],
   {p, 4}, {m, 4}, {n, 4}];

ric4 = Table[ser[
    Sum[D[ch4[[s, m, n]], co4[[s]]], {s, 4}]
  - Sum[D[ch4[[s, m, s]], co4[[n]]], {s, 4}]
  + Sum[ch4[[s, s, p]] ch4[[p, m, n]], {s, 4}, {p, 4}]
  - Sum[ch4[[s, n, p]] ch4[[p, m, s]], {s, 4}, {p, 4}]],
  {m, 4}, {n, 4}];

rr4 = ser[Sum[g4inv[[m, n]] ric4[[m, n]], {m, 4}, {n, 4}]];
sq4 = ser[nlapse Sqrt[Det[h3]]];

(* F = f0fr + eps ff[t] cos kx, with the background conditions imposed *)
ffld = f0fr + eps ff[t] cc;
vpot = 6 f0fr hh^2 + 12 hh^2 (ffld - f0fr) + (1/2) vv2 (ffld - f0fr)^2;

lagFR = ser[sq4 (ffld rr4 - vpot)];

countCov[lag_, vars_] := Module[{l2, need, lred, extra},
  l2 = Simplify[avg[Coefficient[lag, eps, 2]]];
  need = Select[vars, ! FreeQ[l2, Derivative[2][#][t]] &];
  {lred, extra} = If[need === {}, {l2, {}}, OstrogradskyReduce[l2, need]];
  DiracCount[lred, Join[vars, extra]]];

grCovResult = countCov[ser[sq4 (f0fr rr4 - 6 f0fr hh^2)], {psi, phi, b, ee}];
frResult     = countCov[lagFR, {psi, phi, b, ee, ff}];

Dataset @ {
 <|"case" -> "general relativity, covariant route", "coords" -> grCovResult["Coordinates"],
   "constraints" -> grCovResult["Total"], "1st" -> grCovResult["FirstClass"],
   "2nd" -> grCovResult["SecondClass"], "DOF" -> grCovResult["DOF"]|>,
 <|"case" -> "general f(R)", "coords" -> frResult["Coordinates"],
   "constraints" -> frResult["Total"], "1st" -> frResult["FirstClass"],
   "2nd" -> frResult["SecondClass"], "DOF" -> frResult["DOF"]|>}

(*::md::
$f(R)$ propagates **one** scalar mode, the scalaron, and general relativity by the same route
propagates none. The difference is exactly one, which is the answer every textbook gives, and it
is here derived rather than quoted.

The bookkeeping is worth reading. Both cases carry four first-class constraints — the four
diffeomorphisms — and eight second class. What differs is only the number of coordinates: nine
against eight, the extra one being $\delta F$. So

$$
\text{GR} : \frac{2\cdot 8 - 2\cdot 4 - 8}{2} = 0, \qquad
f(R) : \frac{2\cdot 9 - 2\cdot 4 - 8}{2} = 1 .
$$

The eight second-class constraints are themselves informative. They are what removes the
Ostrogradsky variables that ${}^{(4)}R$'s second time derivatives forced into the phase space.
General relativity is not a fourth-order theory, and the constraint algorithm discovers that on
its own — the second derivatives were an artefact of the covariant form, and the algorithm takes
them straight back out. That is a stronger check on this route than the mode count alone: a
mistake in the four-metric or in the Ostrogradsky reduction would show up as a control that no
longer returns zero.
::*)

(*::md::
## 9. Checks
::*)

(*::code::*)
numRules = {hh -> 1, k -> 1, fQ -> 1, fQQ -> 1/3, fQQQ -> 0};
numFq = runCase[numRules];
numLin = runCase[{hh -> 1, k -> 1, fQ -> 1, fQQ -> 0, fQQQ -> 0}];

Dataset @ {
 <|"check" -> "background non-metricity scalar is -6 H^2",
   "ok" -> (Simplify[$Q0 + 6 hh^2] === 0)|>,
 <|"check" -> "the de Sitter background condition is f0 = 2 Q0 fQ",
   "ok" -> (Simplify[(f0 /. bgRule) - 2 $Q0 fQ] === 0)|>,
 <|"check" -> "STEGR satisfies that background condition",
   "ok" -> (Simplify[((f0 /. bgRule) - f0) /. {f0 -> $Q0 - 6 hh^2, fQ -> 1}] === 0)|>,
 <|"check" -> "the second background equation follows from the first",
   "ok" -> (Simplify[bgEvolve /. bgRule] === 0)|>,
 <|"check" -> "every constraint chain terminated",
   "ok" -> (AllTrue[{grResult, stegrResult, linearResult, fqResult}, TrueQ[#["ok"]] &])|>,
 <|"check" -> "ACCEPTANCE: general relativity has no propagating scalar mode",
   "ok" -> (grResult["DOF"] === 0)|>,
 <|"check" -> "ACCEPTANCE: STEGR has no propagating scalar mode",
   "ok" -> (stegrResult["DOF"] === 0)|>,
 <|"check" -> "ACCEPTANCE: the linear branch f'' = 0 has no propagating scalar mode",
   "ok" -> (linearResult["DOF"] === 0)|>,
 <|"check" -> "general f(Q) propagates exactly one scalar mode",
   "ok" -> (fqResult["DOF"] === 1)|>,
 <|"check" -> "the constraint count is unchanged; only the classification moves",
   "ok" -> (fqResult["Total"] === linearResult["Total"] &&
            fqResult["FirstClass"] =!= linearResult["FirstClass"])|>,
 <|"check" -> "exactly two constraints move from first class to second",
   "ok" -> (linearResult["FirstClass"] - fqResult["FirstClass"] === 2 &&
            fqResult["SecondClass"] - linearResult["SecondClass"] === 2)|>,
 <|"check" -> "the connection modes really do carry second time derivatives",
   "ok" -> (! FreeQ[Simplify[avg[Coefficient[lagQ, eps, 2]] /. bgRule],
            Derivative[2][z0][t] | Derivative[2][z1][t]])|>,
 <|"check" -> "numerical spot check, general f(Q) gives one mode",
   "ok" -> (numFq["DOF"] === 1)|>,
 <|"check" -> "numerical spot check, the linear branch gives none",
   "ok" -> (numLin["DOF"] === 0)|>,
 <|"check" -> "the 4D Ricci scalar is 12 H^2 on the de Sitter background",
   "ok" -> (Simplify[(rr4 /. eps -> 0) - 12 hh^2] === 0)|>,
 <|"check" -> "ACCEPTANCE: the covariant route also gives general relativity no scalar mode",
   "ok" -> (grCovResult["DOF"] === 0)|>,
 <|"check" -> "ACCEPTANCE: f(R) propagates exactly one scalar mode",
   "ok" -> (frResult["DOF"] === 1)|>,
 <|"check" -> "f(R) differs from general relativity by one coordinate and nothing else",
   "ok" -> (frResult["Coordinates"] - grCovResult["Coordinates"] === 1 &&
            frResult["FirstClass"] === grCovResult["FirstClass"] &&
            frResult["SecondClass"] === grCovResult["SecondClass"])|>,
 <|"check" -> "the covariant form really does carry second time derivatives",
   "ok" -> (! FreeQ[Simplify[avg[Coefficient[lagFR, eps, 2]]],
            Derivative[2][phi][t] | Derivative[2][ee][t]])|>}
