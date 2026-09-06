(*::md::
# GR-12 — Screening, and the one number the series still quotes

Every other result in this series is derived. One is not. GR-05, GR-11, the README and the
companion fitting library all lean on $|f_{R0}| \lesssim 10^{-6}$, the Solar System bound on
$f(R)$, and every one of them takes it on trust. This notebook derives it.

It matters more after GR-11, which used that bound twice: once to conclude that standard sirens
cannot see $f(R)$, and once to say lensing beats them by a factor of two. Both conclusions are
only as good as the number underneath.

## Why $f(R)$ needs rescuing at all

GR-05 already contains the problem. Its gravitational slip is
$\eta = \Phi/\Psi = (1+2m)/(1+4m)$, and on small scales — $m \to \infty$, a scalaron whose
Compton wavelength dwarfs the system — it goes to $\tfrac{1}{2}$. In a static weak field that
slip *is* the parametrised post-Newtonian $\gamma$. Cassini measures
$\gamma - 1 = (2.1 \pm 2.3)\times10^{-5}$.

So $f(R)$ with a light scalaron is not marginally disfavoured; it is off by four orders of
magnitude. Either the theory is dead or something switches the scalaron off near a dense body.

## What does the switching

The scalaron's mass depends on the local density, and inside a dense body it becomes enormous.
The field is then pinned at its high-density minimum through the interior, and only a **thin
shell** near the surface is free to move and radiate a fifth force. That is the chameleon
mechanism, and this notebook works it out for a uniform sphere, in closed form:

$$\frac{x_s}{r_b} = \sqrt{1 - \varepsilon}, \qquad
\frac{A}{A_{\rm linear}} = 1 - (1-\varepsilon)^{3/2}, \qquad
\varepsilon \equiv \frac{|f_{R,\rm bg}|}{\Phi_N}.$$

A screened core exists only for $\varepsilon < 1$: the field can move by $\tfrac{2}{3}|\Phi_N|$
and no more, so if the background value it has to climb to exceeds that, there is nothing to
pin and the body is bare.

## Scope

Static, spherically symmetric, weak field, uniform density, one body in a uniform background.
The pinned-core idealisation is not exact and section 6 measures how far off it is rather than
waving at it. Chameleon screening only; Vainshtein and symmetron mechanisms are different
theories, not different approximations.
::*)

(*::md::
## 1. Setup

The same backend, and the same convention as the rest of the series: $\kappa = 8\pi G$, so the
Newtonian limit of the $00$ equation is $\nabla^2\Phi_N = \tfrac{\kappa}{2}\rho$.
::*)

(*::code::*)
ClearAll["Global`*"];

$grMetric = ResourceFunction["MetricTensor"];
$grRicci  = ResourceFunction["RicciTensor"];
$grChr    = ResourceFunction["ChristoffelSymbols"];

(*::md::
## 2. The scalaron equation, from the trace

Contracting the $f(R)$ field equations with $g^{\mu\nu}$ gives

$$f_R R - 2f + 3\Box f_R = \kappa T,$$

which is the whole reason $f(R)$ has a propagating scalar: the trace is a wave equation for
$f_R$ rather than the algebraic constraint it is in general relativity, where $f_R = 1$ kills
the $\Box$ term and leaves $-R = \kappa T$.

Given the field equations, the identity is two contractions:
$g^{\mu\nu}g_{\mu\nu} = 4$ turns $-\tfrac{1}{2}fg_{\mu\nu}$ into $-2f$ and $g_{\mu\nu}\Box f_R$
into $4\Box f_R$, and $g^{\mu\nu}\nabla_\mu\nabla_\nu f_R = \Box f_R$ removes one of them again,
leaving $3\Box f_R$.

Those two contractions are what the cell below checks, on the static spherically symmetric weak
field the rest of the notebook uses. They are also where an implementation error would actually
live — the algebra above them is one line. Contracting the assembled field equations instead and
calling `Simplify` on the result is the same statement and does not finish in ten minutes.
::*)

(*::code::*)
coords = {t, r, th, ph};

gStatic = DiagonalMatrix[{
   -(1 + 2 eps pot[r]), 1 - 2 eps psi[r], r^2 (1 - 2 eps psi[r]),
   r^2 Sin[th]^2 (1 - 2 eps psi[r])}];

mt  = $grMetric[gStatic, coords];
gdn = mt["MatrixRepresentation"];
gup = mt["InverseMetricTensor"]["MatrixRepresentation"];
chr = $grChr[mt]["TensorRepresentation"];
ric = $grRicci[mt]["MatrixRepresentation"];
rsc = Sum[gup[[i, j]] ric[[i, j]], {i, 4}, {j, 4}];

ClearAll[hess, boxOf];
hess[f_] := hess[f] = Table[
   D[f, coords[[m]], coords[[q]]] - Sum[chr[[l, m, q]] D[f, coords[[l]]], {l, 4}],
   {m, 4}, {q, 4}];
boxOf[f_] := boxOf[f] = Sum[gup[[m, q]] hess[f][[m, q]], {m, 4}, {q, 4}];

metricTrace = Simplify[Sum[gup[[i, j]] gdn[[i, j]], {i, 4}, {j, 4}]];

(* an ordinary scalar, so the contraction is tested rather than a special case *)
boxContraction = Simplify[
   Sum[gup[[i, j]] hess[q[r]][[i, j]], {i, 4}, {j, 4}] - boxOf[q[r]]];

traceIdentity = (metricTrace === 4) && (boxContraction === 0);

Column[{
  Row[{Style["g^ab g_ab = ", Bold], metricTrace}],
  Row[{Style["g^ab (Hess f)_ab - Box f = ", Bold], boxContraction}],
  Row[{Style["so the trace is fR R - 2f + 3 Box fR - kap T : ", Bold], traceIdentity}]}]

(*::md::
## 3. The Newtonian limit, and how far the field can move

Linearise around a background where the field sits at its minimum. Writing
$f_R = f_{R,\rm bg} + \delta f_R$ and keeping the trace equation to first order,

$$\nabla^2 \delta f_R = m^2\,\delta f_R - \frac{\kappa}{3}\delta\rho ,$$

with $m^2 = \tfrac{1}{3}(f_R/f_{RR} - R)$. Well inside the scalaron's Compton wavelength the
mass term drops and this is a Poisson equation with $-\kappa/3$ where gravity has $+\kappa/2$.
So, sourced by the same matter,

$$\delta f_R = -\frac{2}{3}\Phi_N .$$

That single relation is the whole of screening. The field cannot move further than
$\tfrac{2}{3}|\Phi_N|$, no matter how dense the body, because that is all the source can push
it. If the background value it would have to climb from is larger than that, it never reaches
its interior minimum and the body is unscreened.

Both statements are checked below by solving the two Poisson equations for the same uniform
sphere and taking the ratio, rather than by comparing coefficients by eye.
::*)

(*::code::*)
ClearAll[poisson];

(* interior of a uniform sphere of radius rb, regular at the origin and matched
   to the exterior vacuum solution *)
poisson[source_, rb_] := Module[{cin, cout, win, wout, sol},
  win  = source rr^2/6 + cin;
  wout = -cout/rr;
  sol  = First @ Solve[{(win /. rr -> rb) == (wout /. rr -> rb),
                        (D[win, rr] /. rr -> rb) == (D[wout, rr] /. rr -> rb)}, {cin, cout}];
  {Simplify[win /. sol], Simplify[wout /. sol]}];

{phiIn, phiOut} = poisson[kap rhoIn/2, rb];        (* Laplacian Phi = kap rho / 2 *)
{dfIn, dfOut}   = poisson[-kap rhoIn/3, rb];       (* Laplacian df = -kap rho / 3 *)

ratioIn  = Simplify[dfIn/phiIn];
ratioOut = Simplify[dfOut/phiOut];

Column[{
  Row[{Style["Phi_N inside  : ", Bold], phiIn}],
  Row[{Style["delta fR inside: ", Bold], dfIn}],
  Row[{Style["ratio, inside : ", Bold], ratioIn}],
  Row[{Style["ratio, outside: ", Bold], ratioOut}],
  Row[{Style["both equal -2/3: ", Bold],
       Simplify[{ratioIn + 2/3, ratioOut + 2/3}] === {0, 0}}]}]

(*::md::
## 4. Unscreened is already excluded

With the scalaron light, the extra scalar force adds a third to Newtonian gravity, and the two
potentials split: $\Psi = \Psi_N(1+\delta)$, $\Phi = \Psi_N(1-\delta)$ with $\delta = 1/3$. So

$$\gamma = \frac{\Phi}{\Psi} = \frac{1-\delta}{1+\delta} = \frac{1}{2},$$

which is GR-05's $\eta \to 1/2$ on small scales, arrived at from the static side. Cassini's
$\gamma - 1 = (2.1 \pm 2.3)\times10^{-5}$ is four orders of magnitude away.

Turning that into a requirement: an unscreened fraction $S$ of the fifth force gives
$\delta = S/3$ and $|\gamma - 1| = 2\delta/(1+\delta)$, so $S$ has to be below a few parts in
$10^5$.
::*)

(*::code::*)
gammaOf[d_] := (1 - d)/(1 + d);
deltaOf[S_] := S/3;

cassini = 23/1000000;                                   (* 2.3e-5, one sigma *)
deltaMax = delta /. First @ Solve[2 delta/(1 + delta) == cassini, delta];
sMax = 3 deltaMax;

Column[{
  Row[{Style["gamma, fully unscreened : ", Bold], gammaOf[deltaOf[1]]}],
  Row[{Style["matches GR-05's eta at large m: ", Bold],
       Simplify[gammaOf[deltaOf[1]] - Limit[(1 + 2 m)/(1 + 4 m), m -> Infinity]] === 0}],
  Row[{Style["largest fifth-force fraction Cassini allows: ", Bold],
       ScientificForm[N[sMax], 3]}]}]

(*::md::
## 5. The thin shell

Inside a dense body the scalaron is heavy and the field is pinned at its interior minimum. Take
that literally for now: $\delta f_R$ sits at the minimum for $r < r_s$ and only the shell
$r_s < r < r_b$ responds to matter.

In units $x = r/r_b$ and $w = (1 - f_R)/(1 - f_{R,\rm bg})$ — so $w \to 1$ far away and $w = 0$
at the interior minimum — the shell obeys $w'' + \tfrac{2}{x}w' = S$ with $S = 2/\varepsilon$
and $\varepsilon = |f_{R,\rm bg}|/\Phi_N$. Matching to a pinned core at $x_s$ ($w = w' = 0$)
and to the massless exterior $w = 1 - A/x$ at the surface fixes everything.
::*)

(*::code::*)
ClearAll[shell, xs, A];

Ssrc = 2/epsv;
shellGen = Ssrc x^2/6 + c1 + c2/x;

innerMatch = First @ Solve[
   {(shellGen /. x -> xs) == 0, (D[shellGen, x] /. x -> xs) == 0}, {c1, c2}];
shellSol = Simplify[shellGen /. innerMatch];

surfaceMatch = Simplify @ Solve[
   {(shellSol /. x -> 1) == 1 - A, (D[shellSol, x] /. x -> 1) == A}, {A, xs}];

xsClaim = Sqrt[1 - epsv];
aClaim  = (2/(3 epsv)) (1 - (1 - epsv)^(3/2));

Column[{
  Row[{Style["shell solution w = ", Bold], shellSol}],
  Row[{Style["obeys w'' + 2w'/x = S : ", Bold],
       Simplify[D[shellSol, {x, 2}] + (2/x) D[shellSol, x] - Ssrc] === 0}],
  Row[{Style["pinned at xs, with zero slope: ", Bold],
       Simplify[{shellSol /. x -> xs, D[shellSol, x] /. x -> xs}] === {0, 0}}],
  Row[{Style["matching gives: ", Bold], surfaceMatch}],
  Row[{Style["so xs = Sqrt[1-eps] and A = (2/3eps)(1-(1-eps)^(3/2)): ", Bold],
       Simplify[{(shellSol /. x -> 1) - (1 - aClaim),
                 (D[shellSol, x] /. x -> 1) - aClaim} /. xs -> xsClaim] === {0, 0}}],
  Row[{Style["screening factor A / A_linear = ", Bold],
       Simplify[aClaim/(2/(3 epsv))]}]}]

(*::md::
The two limits are the ones to check. As $\varepsilon \to 1$ the core shrinks to nothing,
$A/A_{\rm linear} \to 1$, and the body is bare. As $\varepsilon \to 0$ the core fills the body,
the shell thins, and the factor goes to $\tfrac{3}{2}\varepsilon$ — linear in how small the
background field is, which is what makes the mechanism work.
::*)

(*::code::*)
screenFactor[e_] := 1 - (1 - e)^(3/2);

Dataset @ Table[
  <|"eps" -> N[e], "xs / rb" -> N[Sqrt[1 - e], 6],
    "shell thickness" -> N[1 - Sqrt[1 - e], 6],
    "A / A_linear" -> N[screenFactor[e], 6],
    "3 eps / 2" -> N[3 e/2, 6]|>,
  {e, {9/10, 1/2, 3/10, 1/10, 3/100, 1/100, 1/1000}}]

(*::md::
## 6. How good is a pinned core?

Not exact, and the honest way to say so is to measure what was dropped. The full interior
equation is $w'' + \tfrac{2}{x}w' = \tfrac{2}{\varepsilon}\left(1 - \tfrac{1}{D\sqrt{w}}\right)$
with $D$ the density contrast between body and background; section 5 kept the $1$ and threw away
the rest. So the neglected piece is a fraction $1/(D\sqrt{w})$ of the source, and it only
matters where $w \lesssim 1/D^2$ — that is, in a boundary layer hugging the pinned core.

Near $x_s$ the closed form goes as $w \approx (x - x_s)^2/\varepsilon$, so the layer reaches out
to $x - x_s = \sqrt{\varepsilon}/D$, and as a fraction of the shell that is

$$\frac{\sqrt{\varepsilon}}{D\left(1 - \sqrt{1-\varepsilon}\right)}.$$

Worth doing this algebraically rather than by integrating the full equation. The interior
minimum is an *equilibrium*, so a numerical solution started exactly on it never leaves; started
just off it, where it leaves depends on rounding. That is a genuinely awkward two-point boundary
value problem, and it answers a question this section can answer exactly.
::*)

(*::code::*)
(* the coefficient is pulled with SeriesCoefficient, not with Coefficient on the
   expanded series: Normal re-expands in powers of x, after which asking for the
   coefficient of (x - xs)^2 truthfully returns zero *)
nearCore = Simplify[Normal[Series[shellSol /. xs -> xsClaim, {x, xsClaim, 2}]]];
nearCoreCoeff = Simplify[
   SeriesCoefficient[shellSol /. xs -> xsClaim, {x, xsClaim, 2}],
   Assumptions -> 0 < epsv < 1];
layerWidth = Sqrt[epsv]/dCon;
layerFraction[e_, d_] := Sqrt[e]/(d (1 - Sqrt[1 - e]));

Column[{
  Row[{Style["w near the core, to leading order: ", Bold], nearCore}],
  Row[{Style["leading coefficient is 1/eps: ", Bold], nearCoreCoeff}],
  Dataset @ Table[
    <|"body" -> b[[1]], "density contrast D" -> ScientificForm[N[b[[2]]], 2],
      "eps" -> N[b[[3]]],
      "boundary layer / shell" -> ScientificForm[N[layerFraction[b[[3]], b[[2]]]], 3]|>,
    {b, {{"laboratory sphere", 10^30, 1/100},
         {"the Sun", 10^30, 1/1000},
         {"the Galaxy", 10^6, 1/10}}}]}]

(*::md::
## 7. The bound

Now assemble it. The Solar System is not screened on its own account — it is screened because
the **Galaxy** is, and a body sitting inside a screened halo sees a local $f_R$ far below the
cosmological one. So the requirement is that the Galaxy has a thin shell at all, which by
section 5 means $\varepsilon_{\rm gal} < 1$:

$$|f_{R0}| \;<\; \Phi_{\rm gal} \;\approx\; 10^{-6}.$$

That is the number GR-05, GR-11, the README and CosmoFit have all been quoting, and it is not a
coincidence that it looks like a potential: it *is* a potential, the depth of the well the
scalaron has to climb out of.

For contrast, the requirement if the Galaxy did **not** screen — the Sun bare in a cosmological
background — is far tighter, and is the reason the mechanism is needed rather than merely
convenient.
::*)

(*::code::*)
phiGalaxy = 10^-6;
phiSun = 21/10^7;

epsFromS[S_] := e /. First @ Solve[screenFactor[e] == S && 0 < e < 1, e, Reals];

galaxyBound = phiGalaxy;                                   (* eps_gal < 1 *)
sunBoundBare = Quiet @ Check[epsFromS[N[sMax]] phiSun, (2/3) N[sMax] phiSun];

Column[{
  Row[{Style["Galaxy potential : ", Bold], ScientificForm[N[phiGalaxy], 2]}],
  Row[{Style["|fR0| < Phi_gal  : ", Bold], ScientificForm[N[galaxyBound], 2]}],
  Row[{Style["if the Sun had to screen itself, unaided: ", Bold],
       ScientificForm[N[(2/3) sMax phiSun], 3]}],
  Row[{Style["how much tighter that would be: ", Bold],
       ScientificForm[N[galaxyBound/((2/3) sMax phiSun)], 3], " times"}]}]

(*::md::
## 8. What it costs the rest of the series

GR-11 concluded that standard sirens cannot see $f(R)$ because the deviation is bounded by
$|f_{R0}|/2$. That conclusion now rests on something derived rather than quoted. Filling in
$|f_{R0}| < 10^{-6}$:

| | size at $z = 0$ |
|---|---|
| lensing, $\Sigma - 1$ | $< 10^{-6}$ |
| sirens, $1 - d_L^{\rm GW}/d_L^{\rm EM}$ | $< 5\times10^{-7}$ |
| growth, $G_{\rm eff}/G - 1$ on small scales | $< \tfrac{4}{3}\times10^{-6}$ |

all of them the same $f_R$, all capped by the same potential well.

And there is a catch worth keeping in view, because it cuts the other way. Screening is a
**local** statement: it suppresses the fifth force near dense bodies and does nothing at all in
the voids and filaments where GR-04 and GR-05 do their work. A model can be invisible in the
Solar System and still move $f\sigma_8$ by percent, which is exactly why the cosmological fits
in the companion library are worth doing rather than being foreclosed by the laboratory. The
bound derived here constrains $|f_{R0}|$, not the shape of $f$.
::*)

(*::md::
## 9. Checks
::*)

(*::code::*)
Dataset @ {
 <|"acceptance test" -> "the two contractions behind the trace identity hold",
   "ok" -> traceIdentity|>,
 <|"check" -> "the metric trace is 4 in four dimensions",
   "ok" -> (metricTrace === 4)|>,
 <|"acceptance test" -> "the massless scalaron responds as -2/3 of the Newtonian potential",
   "ok" -> (Simplify[{ratioIn + 2/3, ratioOut + 2/3}] === {0, 0})|>,
 <|"check" -> "the same ratio holds inside and outside the body",
   "ok" -> (Simplify[ratioIn - ratioOut] === 0)|>,
 <|"acceptance test" -> "an unscreened scalaron gives PPN gamma = 1/2",
   "ok" -> (gammaOf[deltaOf[1]] === 1/2)|>,
 <|"check" -> "which is GR-05's slip in the small-scale limit",
   "ok" -> (Simplify[gammaOf[deltaOf[1]] - Limit[(1 + 2 m)/(1 + 4 m), m -> Infinity]] === 0)|>,
 <|"check" -> "Cassini allows a fifth-force fraction of a few parts in 10^5",
   "ok" -> (10^-5 < N[sMax] < 10^-4)|>,
 <|"check" -> "the shell solution obeys the shell equation",
   "ok" -> (Simplify[D[shellSol, {x, 2}] + (2/x) D[shellSol, x] - Ssrc] === 0)|>,
 <|"check" -> "and is pinned with zero slope at the core radius",
   "ok" -> (Simplify[{shellSol /. x -> xs, D[shellSol, x] /. x -> xs}] === {0, 0})|>,
 <|"acceptance test" -> "matching gives xs = Sqrt[1 - eps]",
   "ok" -> (Simplify[{(shellSol /. x -> 1) - (1 - aClaim),
                      (D[shellSol, x] /. x -> 1) - aClaim} /. xs -> xsClaim] === {0, 0})|>,
 <|"acceptance test" -> "the screening factor is 1 - (1-eps)^(3/2)",
   "ok" -> (Simplify[aClaim/(2/(3 epsv)) - (1 - (1 - epsv)^(3/2))] === 0)|>,
 <|"check" -> "no core survives once eps reaches 1, and the body is bare",
   "ok" -> (Simplify[xsClaim /. epsv -> 1] === 0 && screenFactor[1] === 1)|>,
 <|"check" -> "for a deeply screened body the factor is linear in eps",
   "ok" -> (Simplify[Limit[screenFactor[e]/(3 e/2), e -> 0]] === 1)|>,
 <|"check" -> "screening is monotone: a smaller background field screens harder",
   "ok" -> (screenFactor[9/10] > screenFactor[3/10] > screenFactor[1/100])|>,
 <|"check" -> "the field near the pinned core rises as (x - xs)^2 / eps",
   "ok" -> (Simplify[nearCoreCoeff - 1/epsv] === 0)|>,
 <|"acceptance test" -> "the pinned-core idealisation is negligible for a real body",
   "ok" -> (layerFraction[1/1000, 10^30] < 10^-20)|>,
 <|"check" -> "and is still small for the Galaxy, where the contrast is only 10^6",
   "ok" -> (layerFraction[1/10, 10^6] < 10^-5)|>,
 <|"acceptance test" -> "requiring the Galaxy to screen gives |fR0| < 10^-6",
   "ok" -> (galaxyBound <= 10^-6)|>,
 <|"check" -> "an unaided Sun would be bound far more tightly than the Galaxy route",
   "ok" -> (N[(2/3) sMax phiSun] < N[galaxyBound]/100)|>}
