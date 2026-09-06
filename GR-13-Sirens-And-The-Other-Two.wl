(*::md::
# GR-13 — Sirens, and the other two families

GR-11 asked what a gravitational wave does on the way here, answered it for $f(R)$, and ended
with a flat statement that the answer does **not** carry over: the teleparallel families have
tensor coefficients $f_\mathcal{T}$ and $f_Q$ which are not tied to unity the way $f_R$ is, so
whether sirens can see them is open. This closes it, and the answer is the opposite of $f(R)$'s.

## Where the machinery comes from

GR-08 built the second-order Lagrangians for both families — sixteen tetrad perturbations for
$f(\mathcal{T})$, ten metric plus four Stückelberg for $f(Q)$ — and checked them against the
TEGR and STEGR limits. It then asked those objects one question, the propagation *speed*, and
got $c_{\rm GW}^2 = 1$ exactly.

This notebook asks a different question of the same objects: not the ratio of the gradient term
to the kinetic one, but the **size of the kinetic one**, whose running is the friction. So it
loads GR-08 rather than rebuilding it. A second copy of that index algebra is a second thing to
keep correct, and the first thing it would do on drifting is agree with nothing.

## What comes out

$$\text{f(\mathcal{T}):}\quad \frac{\partial^2 L}{\partial\dot h^2} = -\tfrac{1}{2}f_\mathcal{T}a^3,
\qquad
\text{f(Q):}\quad \frac{\partial^2 L}{\partial\dot h^2} = +\tfrac{1}{2}f_Q a^3,$$

with $c_{\rm GW}^2 = 1$ falling out again as a by-product — GR-08's result, reached from the
coefficient rather than the ratio. Normalising so that general relativity gives one, the
effective Planck mass is $M^2 = -f_\mathcal{T}$ and $M^2 = f_Q$, and GR-11's amplitude argument
applies unchanged:

$$\frac{d_L^{\rm GW}(z)}{d_L^{\rm EM}(z)} = \sqrt{\frac{M^2(0)}{M^2(z)}}.$$

So all three families say the same sentence with a different noun. What differs is how large the
noun is allowed to be, and there the three part company completely.

## The contrast, and why it exists

For $f(R)$ the deviation is bounded by $|f_{R0}|/2 < 5\times10^{-7}$, because GR-12 shows the
Solar System pins $|f_{R0}|$ below $10^{-6}$. Nothing of the sort pins $f_\mathcal{T}$ or $f_Q$,
and the reason is in this series already: the Cassini bound bites $f(R)$ because GR-05's slip
runs to $\eta = 1/2$ on small scales, and it is silent about the other two because GR-06 and
GR-07 find $\eta = 1$ **exactly**, with no scale dependence and nothing to screen.

Fitted to the same data, the $f(\mathcal{T})$ power law moves the siren distance by around
**ten percent**. That is five orders of magnitude above $f(R)$ and comfortably inside what a
third-generation detector is built to measure.

## Scope

Flat FLRW, one Fourier mode, linear order, tensor sector. GR-08's caveat travels with its
Lagrangians: in $f(\mathcal{T})$ the extra Lorentz modes carry no kinetic term around this
background, so linear theory does not describe them. That bounds the *scalar* sector; the
tensor modes here have a healthy kinetic term wherever $f_\mathcal{T} \neq 0$.
::*)

(*::md::
## 1. Loading GR-08

GR-08 begins with `ClearAll["Global`*"]`, so it has to be loaded before anything else is
defined here — otherwise it deletes it. The cell returns a string rather than GR-08's own
closing `Dataset`, so that its checks are not silently counted a second time as this notebook's.
::*)

(*::code::*)
$gr08 = SelectFirst[
   {FileNameJoin[{Directory[], "GR-08-Strong-Coupling.wl"}],
    If[StringQ[$InputFileName] && $InputFileName =!= "",
      FileNameJoin[{DirectoryName[$InputFileName], "GR-08-Strong-Coupling.wl"}], ""],
    Quiet @ Check[FileNameJoin[{NotebookDirectory[], "GR-08-Strong-Coupling.wl"}], ""]},
   StringQ[#] && FileExistsQ[#] &];

Get[$gr08];
"GR-08 loaded from " <> ToString[$gr08]

(*::md::
## 2. The tensor coefficient

`lagT2` and `lagQ2` are GR-08's second-order Lagrangians in general perturbations; `tensorT` and
`tensorQ` are its transverse-traceless substitutions, with the wave along $x$ and both
polarisations in the $y$–$z$ plane. Averaging over a wavelength is the same one-line `TrigReduce`
trick.

What is new is reading the kinetic coefficient itself instead of dividing it out. Note the sign
convention that travels with each: GR-06 puts $\mathcal{T} = +6H^2$ so TEGR is
$f = -\mathcal{T}$ and $f_\mathcal{T} = -1$, while GR-07 puts $Q = -6H^2$ so STEGR is $f = +Q$
and $f_Q = +1$. Both therefore have to give $M^2 = 1$ after normalisation, and that is the check
that the signs were carried correctly rather than cancelled by accident.
::*)

(*::code::*)
ClearAll[reduceTensor, kineticOf, gradientOf];

reduceTensor[lag_, sub_, abb_] :=
  Expand[Simplify[AverageOverWavelength[lag /. sub]] /. abb];

kineticOf[lred_]  := Simplify[D[lred, hp'[t], hp'[t]]];
gradientOf[lred_] := Simplify[Coefficient[Expand[D[lred, hp[t], hp[t]]], k, 2]];

lredT = reduceTensor[lagT2, tensorT, $abbT];
lredQ = reduceTensor[lagQ2, tensorQ, $abbQ];

kinT = kineticOf[lredT];  gradT = gradientOf[lredT];
kinQ = kineticOf[lredQ];  gradQ = gradientOf[lredQ];

speedTfromCoeff = Simplify[-a[t]^2 gradT/kinT];
speedQfromCoeff = Simplify[-a[t]^2 gradQ/kinQ];

(* normalised so that general relativity gives 1 *)
m2T = Simplify[2 kinT/a[t]^3];
m2Q = Simplify[2 kinQ/a[t]^3];

Dataset @ {
  <|"family" -> "f(T)", "kinetic" -> kinT, "gradient k^2" -> gradT,
    "c_GW^2" -> speedTfromCoeff, "M^2" -> m2T, "M^2 in the GR limit" -> Simplify[m2T /. $tegr]|>,
  <|"family" -> "f(Q)", "kinetic" -> kinQ, "gradient k^2" -> gradQ,
    "c_GW^2" -> speedQfromCoeff, "M^2" -> m2Q, "M^2 in the GR limit" -> Simplify[m2Q /. $stegr]|>}

(*::md::
## 3. The same siren formula

GR-11's argument was about the shape of the equation, not about $f(R)$: a tensor mode with
kinetic coefficient $M^2a^3$ and $c_{\rm GW} = 1$ obeys, in conformal time,

$$h'' + \left(2\frac{a'}{a} + \frac{(M^2)'}{M^2}\right)h' + k^2h = 0,$$

whose friction is exactly $2z'/z$ with $z = a M$. So $h \propto 1/(aM)$, the wave carries
$aMh$ unchanged, and the distance an observer assuming general relativity infers is off by
$\sqrt{M^2(0)/M^2(z)}$. Nothing in that used the metric sector, so it transfers.

The reduction is redone here rather than cited, because "it transfers" is the kind of claim
this series checks.
::*)

(*::code::*)
ClearAll[confEq, zed];

zed = A[s] Sqrt[G[s]];
confEq[f_] := D[f[s], {s, 2}] + (2 A'[s]/A[s] + G'[s]/G[s]) D[f[s], s] + k^2 f[s];

reducedEq = Simplify[confEq[Function[ss, u[ss]/(A[ss] Sqrt[G[ss]])]] zed];

frictionGone = Simplify[Coefficient[Expand[reducedEq], u'[s]]] === 0;
leftIsWKB = Simplify[reducedEq - (u''[s] + (k^2 - D[zed, {s, 2}]/zed) u[s])] === 0;

Column[{
  Row[{Style["h = u/(a M) removes the friction: ", Bold], frictionGone}],
  Row[{Style["leaving u'' + (k^2 - z''/z) u = 0: ", Bold], leftIsWKB}]}]

(*::md::
## 4. A model that is actually fitted

The power law is the one the companion library fits to data, and the one
`cosmofit-export.wls` exports. In GR-06's convention, where TEGR is $f = -\mathcal{T}$,

$$f(\mathcal{T}) = -\mathcal{T} + A_0\,\mathcal{T}^{\,b},$$

which is the same theory as GR-02's $\mathcal{T}_s + A_0(-\mathcal{T}_s)^b$ with
$\mathcal{T}_s = -\mathcal{T}$. Its background follows from GR-06's own Friedmann equation
$\kappa\rho = \tfrac{1}{2}(f - 2\mathcal{T}f_\mathcal{T})$, and $A_0$ is fixed by $E(0) = 1$.

That construction is independent of GR-02's `HubbleFunction`, so it can be checked against it:
the two numbers below are what `cosmofit-export.wls` printed from GR-02, and CosmoFit reproduces
them to $0.0\mathrm{e}0$.

$f(Q)$ needs no separate calculation. Writing $u = -Q = 6H^2$, GR-07's power law
$Q + A_0(-Q)^b$ is $-u + A_0u^b$, the same function of the same positive quantity, with the same
$f - 2Xf_X$ and therefore the same $M^2$. That is the background degeneracy the README already
notes, now extended to the tensor sector: **sirens cannot separate $f(\mathcal{T})$ from
$f(Q)$ either**, though both separate from general relativity.
::*)

(*::code::*)
ClearAll[fpl, fplT, constraint, closureA0, Ez, M2norm];

Om = 3/10;

fpl[tt_, a0_, b_] := -tt + a0 tt^b;
fplT[tt_, a0_, b_] := D[fpl[uu, a0, b], uu] /. uu -> tt;

(* kap rho = (f - 2 T fT)/2, with T = 6 E^2 in H0 = 1 units *)
constraint[ez_, zz_, a0_, b_] := With[{tt = 6 ez^2},
   (fpl[tt, a0, b] - 2 tt fplT[tt, a0, b])/2 - 3 Om (1 + zz)^3];

closureA0[b_] := Module[{sols = Solve[constraint[1, 0, a0, b] == 0, a0]},
  If[Length[sols] === 1, a0 /. First[sols], $Failed]];

(* The bracket has to reach past the answer. E ~ Sqrt[Om (1+z)^3] is already 199
   at z = 50, so a fixed ceiling of 60 leaves FindRoot stopping at the edge of
   its own search region -- which it reports, and which still returns a number
   that a check can pass on. Scale the ceiling with the matter-only value. *)
Ez[zz_, b_] := Module[{a0 = closureA0[b], hi},
  hi = 10 Sqrt[Om (1 + zz)^3] + 10;
  ee /. FindRoot[constraint[ee, zz, a0, b] == 0, {ee, 1, 1/50, hi},
     WorkingPrecision -> 30]];

M2norm[zz_, b_] := Module[{a0 = closureA0[b]}, -fplT[6 Ez[zz, b]^2, a0, b]];

fromGR02 = {1/2 -> 1.341568202099828, 1 -> 1.813326892482546};

Dataset @ Table[
  <|"z" -> N[zz], "E(z) here" -> N[Ez[zz, 1/5], 16],
    "E(z) from GR-02" -> (zz /. fromGR02),
    "agree to" -> ScientificForm[N[Abs[Ez[zz, 1/5]/(zz /. fromGR02) - 1]], 2]|>,
  {zz, {1/2, 1}}]

(*::md::
## 5. How large is it?

$M^2 \to 1$ at high redshift for $b < 1$, since the correction carries $\mathcal{T}^{\,b-1}$ and
$\mathcal{T}$ grows. So the ratio saturates, exactly as $f(R)$'s did — and unlike $f(R)$'s, it
saturates somewhere visible.
::*)

(*::code::*)
sirenRatio[zz_, b_] := Sqrt[M2norm[0, b]/M2norm[zz, b]];

sizeTable = Table[
  <|"b" -> N[b], "A0 (closure)" -> N[closureA0[b], 6],
    "M^2(0)" -> N[M2norm[0, b], 6],
    "ratio at z = 1" -> N[sirenRatio[1, b], 6],
    "ratio at z = 3" -> N[sirenRatio[3, b], 6],
    "saturation" -> N[Sqrt[M2norm[0, b]], 6]|>,
  {b, {0, 1/20, 1/10, 1/5, 3/10}}];

Dataset[sizeTable]

(*::md::
Against $f(R)$, where GR-12's bound caps the same quantity at $5\times10^{-7}$:

| | $1 - d_L^{\rm GW}/d_L^{\rm EM}$, largest |
|---|---|
| $f(R)$, $\|f_{R0}\| < 10^{-6}$ | $5\times10^{-7}$ |
| $f(\mathcal{T})$ or $f(Q)$, $b = 0.1$ | $4\times10^{-2}$ |
| $f(\mathcal{T})$ or $f(Q)$, $b = 0.2$ | $1\times10^{-1}$ |

Five orders of magnitude, and the reason is not that the teleparallel families are wilder. It is
that the laboratory constrains them differently. $f(R)$ is caught because its slip runs to
$\eta = 1/2$ on small scales (GR-05), which Cassini excludes, which forces the screening of
GR-12, which forces $|f_{R0}| < 10^{-6}$, which is the very same number that bounds the siren
effect. $f(\mathcal{T})$ and $f(Q)$ have $\eta = 1$ identically (GR-06, GR-07): the chain never
starts, and nothing local pins $f_\mathcal{T}$ or $f_Q$ near one.
::*)

(*::md::
## 6. The constraint this does **not** settle

Section 5 says nothing local pins $f_\mathcal{T}$. That is true of its *value* and false of its
*rate*, and the difference matters enough to state before the ten percent is believed.

GR-06 derives $G_{\rm eff}/G = -1/f_\mathcal{T}$, which in this notebook's normalisation is
$G_{\rm eff} = 1/M^2$. So a running $M^2$ is a running gravitational constant:

$$\alpha_M \equiv \frac{d\ln M^2}{d\ln a} = -\frac{1}{H}\frac{\dot G_{\rm eff}}{G_{\rm eff}}.$$

Lunar laser ranging bounds $|\dot G/G| < 1.6\times10^{-13}\,\mathrm{yr}^{-1}$, and
$H_0 \approx 7.2\times10^{-11}\,\mathrm{yr}^{-1}$, so $|\alpha_M| \lesssim 2\times10^{-3}$ today.
The model of section 5 is nowhere near that, as the cell below measures.

**Whether that excludes it is a question this series cannot answer**, and the reason is exactly
GR-12. Lunar ranging measures $G$ inside the Solar System, so what enters is the *local*
$G_{\rm eff}$, not the cosmological one. For $f(R)$ those two differ by many orders of magnitude,
because the chameleon pins the field locally — that is the whole content of GR-12. Whether
$f(\mathcal{T})$ or $f(Q)$ do anything analogous needs the static spherically symmetric solution
in those theories, the teleparallel counterpart of GR-12, which GR-06 and GR-07 do not build.

So the honest version of section 5's table is narrower than it looks. A siren would see ten
percent **if** the cosmological running survives into the local system — and if it does, the same
running already exceeds what lunar ranging allows by a factor of about seventy at $b = 0.2$, and
by fifteen even at $b = 0.05$. Either way the two experiments are looking at the same function,
which is the useful part: whatever bounds one bounds the other, and $f(R)$'s $10^{-6}$ has no
counterpart here yet in either direction.
::*)

(*::code::*)
ClearAll[alphaM];

(* d ln a = -dz/(1+z), so at z = 0 the running is minus the z-derivative. A
   difference quotient, with the step halved underneath to show it is resolved. *)
alphaM[b_, h_] := -(Log[M2norm[h, b]] - Log[M2norm[0, b]])/h;

llrBound = 22/10000;          (* |Gdot/G| < 1.6e-13 / yr against H0 = 7.2e-11 / yr *)

alphaTable = Table[
  <|"b" -> N[b],
    "alpha_M(0), h = 1/100" -> N[alphaM[b, 1/100], 6],
    "same, h = 1/200" -> N[alphaM[b, 1/200], 6],
    "times the LLR bound" -> N[Abs[alphaM[b, 1/200]]/llrBound, 4]|>,
  {b, {1/20, 1/10, 1/5, 3/10}}];

Column[{
  Dataset[alphaTable],
  Row[{Style["lunar laser ranging allows |alpha_M| up to about ", Bold], N[llrBound]}]}]

(*::md::
## 7. One value of $b$ the background cannot see, and neither can the closure

$f - 2\mathcal{T}f_\mathcal{T} = \mathcal{T} + A_0(1-2b)\mathcal{T}^{\,b}$, so at
$b = \tfrac{1}{2}$ the correction drops out of the Friedmann constraint entirely, for **any**
$A_0$. The background is then matter-only whatever $A_0$ does, $E(0) = 1$ has no solution, and
the closure fails rather than returning something.

It is worth saying out loud what that would look like if the closure were not checked: an
exporter or a fit that took $b$ near $\tfrac12$ would be working with a model whose free
parameter does nothing to the expansion history and everything to $M^2$ — the correction is
still there in $f_\mathcal{T} = -1 + A_0b\mathcal{T}^{\,b-1}$. Background data would be flat in
$A_0$ and the fit would wander.
::*)

(*::code::*)
halfIsDegenerate = closureA0[1/2] === $Failed;

Column[{
  Row[{Style["f - 2 T fT = ", Bold],
       Simplify[fpl[tt, a0, b] - 2 tt fplT[tt, a0, b]]}],
  Row[{Style["the A0 term carries (1 - 2b), so it vanishes at b = 1/2: ", Bold],
       Simplify[(fpl[tt, a0, 1/2] - 2 tt fplT[tt, a0, 1/2]) - tt] === 0}],
  Row[{Style["and the closure has no solution there: ", Bold], halfIsDegenerate}],
  Row[{Style["while M^2 still depends on A0: ", Bold],
       ! FreeQ[Simplify[-fplT[tt, a0, 1/2]], a0]}]}]

(*::md::
## 7. Checks
::*)

(*::code::*)
ratios = Association[#["b"] -> #["saturation"] & /@ sizeTable];
alphas = Association[#["b"] -> #["same, h = 1/200"] & /@ alphaTable];

Dataset @ {
 <|"check" -> "GR-08 was loaded and supplied its second-order Lagrangians",
   "ok" -> (Head[lagT2] =!= Symbol && Head[lagQ2] =!= Symbol)|>,
 <|"acceptance test" -> "the f(T) tensor kinetic coefficient is -fT a^3 / 2",
   "ok" -> (Simplify[kinT + fT a[t]^3/2] === 0)|>,
 <|"acceptance test" -> "the f(Q) tensor kinetic coefficient is +fQ a^3 / 2",
   "ok" -> (Simplify[kinQ - fQ a[t]^3/2] === 0)|>,
 <|"acceptance test" -> "c_GW^2 = 1 for both, read off the coefficients",
   "ok" -> (Simplify[speedTfromCoeff - 1] === 0 && Simplify[speedQfromCoeff - 1] === 0)|>,
 <|"check" -> "which is what GR-08 got by its own route",
   "ok" -> (Simplify[speedTfromCoeff - speedT] === 0)|>,
 <|"acceptance test" -> "M^2 is 1 in both GR limits, so the two sign conventions cancel",
   "ok" -> (Simplify[m2T /. $tegr] === 1 && Simplify[m2Q /. $stegr] === 1)|>,
 <|"check" -> "M^2 is -fT for torsion and +fQ for non-metricity",
   "ok" -> (Simplify[m2T + fT] === 0 && Simplify[m2Q - fQ] === 0)|>,
 <|"check" -> "the tensor mode has a healthy kinetic term away from fT = 0",
   "ok" -> (Simplify[kinT /. fT -> -1] === a[t]^3/2)|>,
 <|"acceptance test" -> "h = u/(a M) removes the friction exactly",
   "ok" -> frictionGone|>,
 <|"check" -> "and leaves u'' + (k^2 - z''/z) u = 0",
   "ok" -> leftIsWKB|>,
 <|"acceptance test" -> "the background here agrees with GR-02's, to twelve digits",
   "ok" -> AllTrue[{1/2, 1}, Abs[Ez[#, 1/5]/(# /. fromGR02) - 1] < 10^-12 &]|>,
 <|"check" -> "b = 0 is exactly general relativity, so M^2 = 1 and the ratio is 1",
   "ok" -> (Simplify[M2norm[0, 0] - 1] === 0 && Simplify[sirenRatio[3, 0] - 1] === 0)|>,
 <|"check" -> "M^2 tends back to 1 at high redshift, so the effect saturates",
   "ok" -> (Abs[M2norm[50, 1/5] - 1] < 1/100)|>,
 <|"check" -> "the effect grows with b",
   "ok" -> (ratios[0.05] < ratios[0.1] < ratios[0.2] < ratios[0.3])|>,
 <|"acceptance test" -> "at b = 0.2 the siren distance moves by more than five percent",
   "ok" -> (ratios[0.2] - 1 > 1/20)|>,
 <|"acceptance test" -> "which beats the f(R) bound of |fR0|/2 by five orders of magnitude",
   "ok" -> ((ratios[0.2] - 1)/(10^-6/2) > 10^4)|>,
 <|"check" -> "the running alpha_M is resolved: halving the step barely moves it",
   "ok" -> AllTrue[alphaTable,
      Abs[#["same, h = 1/200"]/#["alpha_M(0), h = 1/100"] - 1] < 1/20 &]|>,
 <|"check" -> "alpha_M grows with b and vanishes as b does",
   "ok" -> (alphas[0.05] < alphas[0.1] < alphas[0.2] < alphas[0.3]
            && alphas[0.05] < 1/10)|>,
 <|"acceptance test" -> "the same running is far above what lunar ranging allows",
   "ok" -> (Abs[alphas[0.2]]/llrBound > 20)|>,
 <|"check" -> "even the mildest exponent in the table exceeds that bound",
   "ok" -> (Abs[alphas[0.05]]/llrBound > 5)|>,
 <|"check" -> "at b = 1/2 the correction leaves the Friedmann constraint",
   "ok" -> (Simplify[(fpl[tt, a0, 1/2] - 2 tt fplT[tt, a0, 1/2]) - tt] === 0)|>,
 <|"check" -> "so the closure fails there rather than returning a value",
   "ok" -> halfIsDegenerate|>,
 <|"check" -> "though M^2 still depends on A0 at that b, which is why it would mislead",
   "ok" -> (! FreeQ[Simplify[-fplT[tt, a0, 1/2]], a0])|>}
