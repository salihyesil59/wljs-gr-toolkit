(*::md::
# GR-11 — How a gravitational wave travels

GR-08 asked whether the tensor modes propagate and how fast, and answered
$c_{\rm GW}^2 = 1$ identically for $f(\mathcal{T})$ and $f(Q)$. Speed is half the story. A wave
also gets **damped** on the way here, and if it is damped by more than the expansion alone, the
amplitude that arrives is smaller than the distance would suggest — so a standard siren reports
the wrong distance.

That is a different kind of observable from anything else in the series. GR-03 measures
integrals of $E(z)$; GR-04 and GR-05 measure how structure grows. This one measures what
happens to a wave in transit, and it is sensitive to something neither of the others is.

## What this notebook adds

$f(R)$ is the gap. GR-08 covers the tensor sector of the two teleparallel families and does not
touch $f(R)$ at all. So the tensor mode is taken through the same engine GR-05 and GR-10 use —
perturb the field equations, expand in $\epsilon$ — and what comes out is

$$\ddot h + \left(3H + \frac{\dot f_R}{f_R}\right)\dot h + \frac{k^2}{a^2}h = 0.$$

Two things to notice. The $k^2/a^2$ coefficient is exactly one relative to $\ddot h$, so
$c_{\rm GW}^2 = 1$ **identically** — $f(R)$ joins the other two, and GW170817 says nothing about
any of them. And the friction is not $3H$: there is an extra $\dot f_R/f_R$, the running of the
same $f_R$ that GR-05 found rescaling $G_{\rm eff}$ and $\Sigma$.

## Where it ends up

Reducing the equation gives a conserved amplitude, $a\sqrt{f_R}\,h$ rather than $a h$, and
therefore

$$\frac{d_L^{\rm GW}(z)}{d_L^{\rm EM}(z)} = \sqrt{\frac{f_R(0)}{f_R(z)}}.$$

For Hu–Sawicki this is monotonic and saturates, and the whole effect is bounded by
$|f_{R0}|/2$ — **the same parameter the Solar System already bounds at $10^{-6}$.** So the
answer for $f(R)$ is a clean negative: sirens cannot see it, and cannot see it for a reason
rather than by accident. Lensing does slightly better, by a factor of two, on the same
parameter.

## Scope

Metric perturbation theory, so $f(R)$ only; $f(\mathcal{T})$ and $f(Q)$ need the tetrad and
connection perturbations GR-06 and GR-07 build, and section 9 says exactly what would have to
be done. Flat FLRW, one Fourier mode, linear order. The propagation effect only; a real siren
calculation also has to ask whether the source's own emission is modified, which is a question
about screening near the binary rather than about the wave.
::*)

(*::md::
## 1. Setup

The same curvature backend as the rest of the series, and the same trick: `F` is left undefined,
so the derivation holds for every $f(R)$ at once and $f_R$, $f_{RR}$ appear on their own.
::*)

(*::code::*)
ClearAll["Global`*"];

$grMetric = ResourceFunction["MetricTensor"];
$grRicci  = ResourceFunction["RicciTensor"];
$grChr    = ResourceFunction["ChristoffelSymbols"];

$coords = {t, x, y, z};
$wave   = Exp[I k x];

(*::md::
## 2. The transverse traceless mode

The wave runs along $x$, so both polarisations live in the $y$–$z$ block: $h_+$ on the diagonal
with opposite signs, $h_\times$ off it. Transverse because the perturbation has no $x$ index,
traceless because the two diagonal entries cancel.

That is what makes the tensor sector easy where the scalar sector was not. A traceless
perturbation of a conformally flat background does not change the Ricci scalar, so
$\delta R = 0$ — checked below, not assumed — and therefore $\delta f_R = 0$ as well. The
scalaron that made GR-10 fourth order simply is not excited.
::*)

(*::code::*)
gPert = {
  {-1, 0, 0, 0},
  {0, a[t]^2, 0, 0},
  {0, 0, a[t]^2 (1 + eps hp[t] $wave), a[t]^2 eps hc[t] $wave},
  {0, 0, a[t]^2 eps hc[t] $wave, a[t]^2 (1 - eps hp[t] $wave)}};

mt  = $grMetric[gPert, $coords];
gdn = mt["MatrixRepresentation"];
gup = mt["InverseMetricTensor"]["MatrixRepresentation"];
chr = $grChr[mt]["TensorRepresentation"];
ric = $grRicci[mt]["MatrixRepresentation"];
rsc = Sum[gup[[i, j]] ric[[i, j]], {i, 4}, {j, 4}];

MatrixForm[gPert]

(*::md::
## 3. Field equations

Matter is a perfect fluid with both density and pressure. It carries no anisotropic stress, so
it should not source the tensor modes at all — but it is put in rather than left out, so that
"no source" is something the calculation says rather than something the setup assumed.
::*)

(*::code::*)
ClearAll[hess, boxOf];
hess[f_] := hess[f] = Table[
   D[f, $coords[[m]], $coords[[q]]]
     - Sum[chr[[l, m, q]] D[f, $coords[[l]]], {l, 4}], {m, 4}, {q, 4}];
boxOf[f_] := boxOf[f] = Sum[gup[[m, q]] hess[f][[m, q]], {m, 4}, {q, 4}];

uUp = {1, 0, 0, 0};
uDn = Table[Sum[gdn[[i, j]] uUp[[j]], {j, 4}], {i, 4}];
tmn = Table[(rho[t] + pres[t]) uDn[[i]] uDn[[j]] + pres[t] gdn[[i, j]], {i, 4}, {j, 4}];

fieldEq = With[{fv = F[rsc], fr = F'[rsc]},
  Table[fr ric[[i, j]] - (1/2) fv gdn[[i, j]] + gdn[[i, j]] boxOf[fr]
      - hess[fr][[i, j]] - kap tmn[[i, j]], {i, 4}, {j, 4}]];

AbsoluteTiming[
 bgR    = Simplify[Normal[Series[rsc, {eps, 0, 0}]]];
 deltaR = Simplify[Normal[Series[rsc, {eps, 0, 1}]] - bgR];
 order1 = Table[Coefficient[Normal[Series[fieldEq[[i, j]], {eps, 0, 1}]], eps], {i, 3, 4}, {j, 3, 4}];
 bgxx   = Simplify[Normal[Series[fieldEq[[2, 2]], {eps, 0, 0}]]];
]

(*::code::*)
Column[{
  Row[{Style["background R = ", Bold], bgR}],
  Row[{Style["delta R at first order = ", Bold], deltaR}],
  Row[{Style["so the scalaron is not excited: ", Bold], deltaR === 0}]}]

(*::md::
## 4. Putting the background on shell

The $34$ component at first order is not yet the wave equation. It still contains
(background $\mathcal{E}_{yy}$) $\times\, h$, because $g_{34} = \epsilon a^2 h_\times$ multiplies
the background equation. Those pieces vanish only when the background equations hold, so they
have to be imposed rather than hoped for: solve the spatial background equation for the pressure
and substitute.

Skip this and the answer is a mess of $f_0$, $\rho$ and $p$ terms that look like matter sourcing
the tensor mode. It is not; it is the background equation in disguise.
::*)

(*::code::*)
onShell = First @ Solve[bgxx == 0, pres[t]];

abbrev = {Derivative[3][F][bgR] -> fRRR, Derivative[2][F][bgR] -> fRR,
          Derivative[1][F][bgR] -> fR, F[bgR] -> f0};

eqCross = Simplify[Expand[(order1[[1, 2]] /. onShell)/$wave]];
eqPlus  = Simplify[Expand[((order1[[1, 1]] - order1[[2, 2]]) /. onShell)/(2 $wave)]];

Column[{
  Row[{Style["cross equation: ", Bold], Simplify[eqCross /. abbrev], " = 0"}],
  Row[{Style["free of rho and p, so matter does not source it: ", Bold],
       FreeQ[eqCross, rho | pres]}],
  Row[{Style["both polarisations obey the same equation: ", Bold],
       Simplify[(eqCross /. hc -> hh) - (eqPlus /. hp -> hh)] === 0}]}]

(*::md::
## 5. The wave equation

Written out, the $34$ equation is exactly $\tfrac{a^2}{2}$ times

$$f_R\left(\ddot h + 3H\dot h + \frac{k^2}{a^2}h\right) + \dot f_R\,\dot h = 0,$$

that is

$$\ddot h + \left(3H + \frac{\dot f_R}{f_R}\right)\dot h + \frac{k^2}{a^2}h = 0.$$

The gradient term and the acceleration term carry the *same* factor $f_R$, which is what
$c_{\rm GW}^2 = 1$ means. It holds for every $f$, exactly, not to some order — the same
conclusion GR-08 reached for $f(\mathcal{T})$ and $f(Q)$ by a completely different route.

What is not the same as general relativity is the friction. Writing
$\alpha_M \equiv \dfrac{d\ln f_R}{d\ln a}$, the damping is $(3 + \alpha_M)H$ instead of $3H$.
::*)

(*::code::*)
target = (a[t]^2/2) (
   F'[bgR] (hc''[t] + 3 (a'[t]/a[t]) hc'[t] + k^2 hc[t]/a[t]^2)
   + D[F'[bgR], t] hc'[t]);

waveEqMatches = Simplify[eqCross - target] === 0;

cGW2 = Simplify[
  Coefficient[Expand[eqCross /. abbrev], hc[t]] a[t]^2/
  (k^2 Coefficient[Expand[eqCross /. abbrev], hc''[t]])];

grLimit = Simplify[(eqCross /. {Derivative[_][F][_] -> 0, F -> (# &)}) ];

Column[{
  Row[{Style["equation is fR (h'' + 3H h' + k^2 h/a^2) + fR' h', times a^2/2 : ", Bold],
       waveEqMatches}],
  Row[{Style["c_GW^2 = ", Bold], cGW2}],
  Row[{Style["in general relativity (f = R) the friction is just 3H: ", Bold],
       Simplify[(eqCross /. abbrev /. {fR -> 1, fRR -> 0, fRRR -> 0})
         - (a[t]^2/2)(hc''[t] + 3 (a'[t]/a[t]) hc'[t] + k^2 hc[t]/a[t]^2)] === 0}]}]

(*::md::
## 6. What is conserved

Friction changes the amplitude, so the question is what quantity the wave carries unchanged.
Conformal time is where that is exact. With $dt = a\,ds$ the equation becomes

$$h'' + \left(2\frac{a'}{a} + \frac{f_R'}{f_R}\right)h' + k^2 h = 0,$$

and the friction coefficient is precisely $2z'/z$ with $z = a\sqrt{f_R}$ — which is the
signature of an equation that becomes friction-free under $h = u/z$.

The conversion is done in the direction that needs no chain rule inside a replacement: start
from the conformal form, go to cosmic time, and check it is the same equation. Doing it the
other way round produces $\tau$-derivatives after the substitution rules have already fired, and
they survive into the answer looking like extra physics.
::*)

(*::code::*)
conformal = hh2 + (2 A1/A0 + G1/G0) hh1 + k^2 hh0;

toCosmic = {hh2 -> a[t] (a'[t] hd[t] + a[t] hdd[t]), hh1 -> a[t] hd[t], hh0 -> h[t],
            A1/A0 -> a'[t], G1/G0 -> a[t] FR'[t]/FR[t]};

cosmicFromConformal = Simplify[conformal //. toCosmic];
cosmicDirect = a[t]^2 (hdd[t] + 3 (a'[t]/a[t]) hd[t] + (FR'[t]/FR[t]) hd[t] + k^2 h[t]/a[t]^2);

sameEquation = Simplify[cosmicFromConformal - cosmicDirect] === 0;

zed = A[s] Sqrt[G[s]];
confEq[f_] := D[f[s], {s, 2}] + (2 A'[s]/A[s] + G'[s]/G[s]) D[f[s], s] + k^2 f[s];
reducedEq = Simplify[confEq[Function[ss, u[ss]/(A[ss] Sqrt[G[ss]])]] zed];

frictionGone = Simplify[Coefficient[Expand[reducedEq], u'[s]]] === 0;
leftIsWKB = Simplify[reducedEq - (u''[s] + (k^2 - D[zed, {s, 2}]/zed) u[s])] === 0;

Column[{
  Row[{Style["conformal form is the same equation: ", Bold], sameEquation}],
  Row[{Style["h = u/(a Sqrt[fR]) removes the friction: ", Bold], frictionGone}],
  Row[{Style["what is left is u'' + (k^2 - z''/z) u = 0: ", Bold], leftIsWKB}]}]

(*::md::
Well inside the horizon $k^2 \gg z''/z$, so $u$ oscillates with constant amplitude and

$$h \propto \frac{1}{a\sqrt{f_R}}$$

instead of the $1/a$ of general relativity. In general relativity $f_R = 1$ and the two agree,
which is the check that this is a modification and not a change of variables.
::*)

(*::md::
## 7. The same statement, numerically

The reduction above is exact algebra; that $h$ really tracks $1/(a\sqrt{f_R})$ is the WKB
approximation on top of it, and approximations are the part of this series that get tested.

So: integrate the wave equation for a mode well inside the horizon and follow its envelope. For
an oscillator the WKB amplitude is $\sqrt{h^2 + (\dot h/\omega)^2}$ with $\omega = k/a$, no
peak-finding needed. Multiply it by $a$ and it should drift; multiply it by $a\sqrt{f_R}$ and it
should not.

$|f_{R0}| = 0.1$ here on purpose. A viable model would put the whole effect in the fifth decimal
place, where it could not be told from integration error.
::*)

(*::code::*)
Om = 3/10; OL = 7/10; Lam = 3 OL;
R0 = 3 (Om + 4 OL);
aFun = Function[tt, (Om/OL)^(1/3) Sinh[(3/2) Sqrt[OL] tt]^(2/3)];
t0 = (2/(3 Sqrt[OL])) ArcSinh[Sqrt[OL/Om]];
tOfa[av_] := (2/(3 Sqrt[OL])) ArcSinh[Sqrt[OL/Om] av^(3/2)];
Rb[tt_] := 3 (Om/aFun[tt]^3 + 4 OL);

fRat[tt_, xx_] := 1 - xx R0^2/Rb[tt]^2;

kTest = 300; xTest = 1/10; tStart = tOfa[1/10];

waveSol = First @ NDSolve[{
   h''[t] + (3 aFun'[t]/aFun[t] + D[fRat[t, xTest], t]/fRat[t, xTest]) h'[t]
     + (kTest^2/aFun[t]^2) h[t] == 0,
   h[tStart] == 1, h'[tStart] == 0}, h, {t, tStart, t0},
   MaxSteps -> 10^7, AccuracyGoal -> 12, PrecisionGoal -> 12];

wkbAmp[tt_] := Sqrt[(h[t] /. waveSol /. t -> tt)^2
   + ((aFun[tt]/kTest) (h'[t] /. waveSol /. t -> tt))^2];

envelope = Table[
  Module[{tt = tOfa[av], amp}, amp = wkbAmp[tt];
   <|"a" -> N[av], "amplitude" -> amp,
     "x a" -> amp aFun[tt], "x a Sqrt[fR]" -> amp aFun[tt] Sqrt[fRat[tt, xTest]]|>],
  {av, {1/10, 1/4, 1/2, 3/4, 1}}];

spread[key_] := Module[{v = #[key] & /@ envelope}, (Max[v] - Min[v])/Mean[v]];

Column[{Dataset[envelope],
  Row[{Style["spread of a h        : ", Bold], ScientificForm[spread["x a"], 3]}],
  Row[{Style["spread of a Sqrt[fR] h: ", Bold], ScientificForm[spread["x a Sqrt[fR]"], 3]}]}]

(*::md::
## 8. From amplitude to distance

A siren is read by comparing the amplitude that arrives with the amplitude the source produced.
Since $a\sqrt{f_R}\,h$ is what the wave carries unchanged,

$$h_{\rm obs} = h_{\rm emit}\,\frac{a_e\sqrt{f_R(a_e)}}{a_o\sqrt{f_R(a_o)}},$$

while general relativity would have given $h_{\rm emit}\,a_e/a_o$. An observer who assumes
general relativity converts amplitude to distance with $h \propto 1/d_L$, so the distance they
report is off by the ratio of the two:

$$\frac{d_L^{\rm GW}(z)}{d_L^{\rm EM}(z)} = \sqrt{\frac{f_R(0)}{f_R(z)}}.$$

Since $f_R \to 1$ at high curvature and $f_R(0) = 1 - |f_{R0}|$, the ratio is below one, falls
with redshift, and **saturates**: there is nothing left to accumulate once $f_R(z)$ has reached
one. The whole effect is therefore bounded, for every redshift at once, by

$$1 - \sqrt{1 - |f_{R0}|} \;\approx\; \frac{|f_{R0}|}{2}.$$
::*)

(*::code::*)
ClearAll[sirenRatio];
sirenRatio[zz_, xx_] := Sqrt[fRat[t0, xx]/fRat[tOfa[1/(1 + zz)], xx]];

sirenTable = Table[
  <|"z" -> N[zz],
    "|fR0| = 1e-4" -> N[sirenRatio[zz, 1/10000], 12],
    "|fR0| = 1e-6" -> N[sirenRatio[zz, 1/1000000], 12]|>,
  {zz, {1/10, 1/2, 1, 2, 5}}];

bound[xx_] := 1 - Sqrt[1 - xx];

Column[{
  Dataset[sirenTable],
  Row[{Style["saturation value at |fR0| = 1e-4 : ", Bold], N[Sqrt[1 - 1/10000], 12]}],
  Row[{Style["largest deviation, |fR0| = 1e-4  : ", Bold], ScientificForm[N[bound[1/10000]], 4]}],
  Row[{Style["largest deviation, |fR0| = 1e-6  : ", Bold], ScientificForm[N[bound[1/1000000]], 4]}],
  Row[{Style["against |fR0|/2                  : ", Bold],
       ScientificForm[N[1/1000000/2], 4]}]}]

(*::md::
## 9. What this says, and what it says about the other two

Three observables in this series now depend on the same single function. GR-05 found
$G_{\rm eff}/G = (1/f_R)(1+4m)/(1+3m)$ for growth and $\Sigma = 1/f_R$ for lensing, with no
scale dependence in the second at all; this notebook adds
$d_L^{\rm GW}/d_L^{\rm EM} = \sqrt{f_R(0)/f_R(z)}$. At $z=0$, with $f_R = 1 - |f_{R0}|$:

| | deviation from general relativity |
|---|---|
| lensing, $\Sigma - 1$ | $\approx |f_{R0}|$ |
| sirens, $1 - d_L^{\rm GW}/d_L^{\rm EM}$ | $\approx |f_{R0}|/2$ |

so lensing is the more sensitive of the two by a factor of two, on exactly the same parameter.
And that parameter is the one the Solar System already pins: at $|f_{R0}| \lesssim 10^{-6}$ the
siren effect is below $5\times10^{-7}$, against the percent-level distance precision a
third-generation detector hopes for. **Standard sirens cannot constrain viable $f(R)$**, and not
by bad luck — the same $f_R \approx 1$ that lets the theory survive a laboratory forces the wave
to travel as it does in general relativity.

That is worth stating precisely because it does **not** follow for the other two families.
GR-08 computed their tensor kinetic coefficients, and they are $f_\mathcal{T}$ and $f_Q$ — the
same coefficients that appear in $G_{\rm eff}/G = -1/f_\mathcal{T}$ and $1/f_Q$, and those are
not tied to unity by a laboratory bound in the way $f_R$ is. Whether sirens can see them is
therefore an open question rather than a settled one, and it is the obvious next calculation.

What would be needed: the tensor sector of $f(\mathcal{T})$ perturbs a tetrad, not a metric, and
of $f(Q)$ a metric plus a flat connection. The transverse traceless modes do not excite the
Lorentz modes of GR-06 or the Stückelberg modes of GR-07 — that is the simplification which
makes it tractable — but the machinery still has to come from those notebooks rather than from
this one. GR-05 draws the same line for the same reason.
::*)

(*::md::
## 10. Checks
::*)

(*::code::*)
Dataset @ {
 <|"check" -> "a transverse traceless mode does not perturb the Ricci scalar",
   "ok" -> (deltaR === 0)|>,
 <|"check" -> "matter carries no anisotropic stress, so it does not source the tensor mode",
   "ok" -> FreeQ[eqCross, rho | pres]|>,
 <|"check" -> "both polarisations obey the same equation",
   "ok" -> (Simplify[(eqCross /. hc -> hh) - (eqPlus /. hp -> hh)] === 0)|>,
 <|"acceptance test" -> "the wave equation is fR (h'' + 3H h' + k^2 h/a^2) + fR' h'",
   "ok" -> waveEqMatches|>,
 <|"acceptance test" -> "c_GW^2 = 1 identically, for every f",
   "ok" -> (Simplify[cGW2 - 1] === 0)|>,
 <|"acceptance test" -> "general relativity gives back friction 3H and nothing else",
   "ok" -> (Simplify[(eqCross /. abbrev /. {fR -> 1, fRR -> 0, fRRR -> 0})
      - (a[t]^2/2)(hc''[t] + 3 (a'[t]/a[t]) hc'[t] + k^2 hc[t]/a[t]^2)] === 0)|>,
 <|"check" -> "the conformal-time form is the same equation",
   "ok" -> sameEquation|>,
 <|"check" -> "h = u/(a Sqrt[fR]) removes the friction exactly",
   "ok" -> frictionGone|>,
 <|"check" -> "and leaves u'' + (k^2 - z''/z) u = 0",
   "ok" -> leftIsWKB|>,
 <|"acceptance test" -> "numerically, a Sqrt[fR] h is far flatter than a h",
   "ok" -> (spread["x a Sqrt[fR]"] < spread["x a"]/5)|>,
 <|"check" -> "the numerical invariant is flat to better than a percent",
   "ok" -> (spread["x a Sqrt[fR]"] < 1/100)|>,
 <|"check" -> "the siren ratio is below one and falls with redshift",
   "ok" -> (sirenRatio[1/10, 1/10000] > sirenRatio[1, 1/10000] > sirenRatio[5, 1/10000]
            && sirenRatio[1/10, 1/10000] < 1)|>,
 <|"check" -> "it saturates at Sqrt[1 - |fR0|] rather than growing without limit",
   "ok" -> (Abs[sirenRatio[1000, 1/10000] - Sqrt[1 - 1/10000]] < 10^-8)|>,
 <|"check" -> "the largest deviation is |fR0|/2 to within a part in a thousand",
   "ok" -> (Abs[bound[1/1000000]/(1/1000000/2) - 1] < 1/1000)|>,
 <|"acceptance test" -> "at the Solar System bound the siren effect is below 1e-6",
   "ok" -> (bound[1/1000000] < 10^-6)|>,
 <|"check" -> "lensing beats sirens by a factor of two on the same parameter",
   "ok" -> (Abs[(1/fRat[t0, 1/1000000] - 1)/bound[1/1000000] - 2] < 1/100)|>}
