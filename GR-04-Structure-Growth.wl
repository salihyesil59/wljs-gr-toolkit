(* ::Package:: *)

(* ============================================================================
   GR-04 : Growth of Structure and f-sigma-8
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-04-Structure-Growth.wl

   ============================================================================ *)

(*::md::
# Growth of Structure

GR-03 turns $E(z)$ into distances, and distances are where supernovae and BAO live. But a
modified-gravity model can be tuned to reproduce the $\Lambda$CDM expansion history almost
exactly, and then every distance in GR-03 agrees to within the error bars. The expansion
history alone cannot separate those theories.

What separates them is how fast lumps grow inside that expansion. Gravity that is slightly
stronger than Newton's pulls matter together faster, and the difference accumulates over
billions of years into something redshift surveys can measure. That observable is
$f\sigma_8(z)$, and this notebook computes it.

## The equation

For sub-horizon linear perturbations in the matter era, the density contrast
$\delta_m = \delta\rho_m/\rho_m$ obeys

$$
\ddot\delta_m + 2H\dot\delta_m - 4\pi G_{\text{eff}}\rho_m\delta_m = 0 .
$$

Trading time for redshift with $\dot{} = -(1+z)H\,d/dz$ and using
$4\pi G\rho_m/H^2 = \tfrac{3}{2}\Omega_m(z)$ with
$\Omega_m(z) = \Omega_{m0}(1+z)^3/E^2$, this becomes an ODE in $z$ alone:

$$
\delta'' + \left[\frac{E'}{E} - \frac{1}{1+z}\right]\delta' \;=\; \frac{3\,\Omega_{m0}(1+z)}{2\,E(z)^2}\,\frac{G_{\text{eff}}}{G}\,\delta .
$$

Only two things go in: the background $E(z)$, which GR-02 gives you, and the ratio
$G_{\text{eff}}/G$. From the solution come the three things people quote:

$$
D(z) = \frac{\delta(z)}{\delta(0)}, \qquad
f(z) = \frac{d\ln\delta}{d\ln a} = -\frac{(1+z)\,\delta'}{\delta}, \qquad
f\sigma_8(z) = -\,\sigma_{8,0}\,\frac{(1+z)\,\delta'(z)}{\delta(0)} .
$$

## What is derived here and what is not

The growth equation above is derived in the notebook, and everything downstream of it follows
by integration. $G_{\text{eff}}/G$ is **input**, not output, and in general relativity it is
exactly 1.

Where the modified-gravity entries come from splits in two. The $f(R)$ forms are derived from
scratch in **GR-05**, which perturbs the field equations in Newtonian gauge and solves for the
metric potentials; its result drops straight into the library below. The $f(\mathcal{T})$ and
$f(Q)$ forms are quoted from the literature, because their fundamental variables are a tetrad
and a flat non-metric connection, which metric perturbation theory cannot reach. Check the
convention against your own source for those two.

One structural difference is worth keeping in view: in $f(R)$ gravity $G_{\text{eff}}$ depends
on the wavenumber $k$, while in $f(\mathcal{T})$ and $f(Q)$ it does not. That scale dependence
is itself a discriminator, and GR-05 plots it.

Conventions follow GR-02 and GR-03: $1+z = 1/a$ with $a_0 = 1$, and $\sigma_{8,0} = 0.811$ by
default.

## The series

| | |
|---|---|
| **GR-01** | a metric in, curvature tensors out |
| **GR-02** | an action in, field equations, Friedmann and $E(z)$ out |
| **GR-03** | $E(z)$ in, distances, times and the BAO ruler out |
| **GR-04** | $E(z)$ and $G_{\text{eff}}$ in, growth and $f\sigma_8$ out |
| **GR-05** | perturbed field equations in, $G_{\text{eff}}$, slip and lensing out |
| **GR-06** | a tetrad in, torsion, the TEGR identity and $f(\mathcal{T})$ cosmology out |

You are on the growth branch, the one distances cannot substitute for. GR-03 takes the same
$E(z)$ and turns it into what a telescope measures; the models worth arguing about are the ones
those two branches disagree about.
::*)

(*::md::
## 1. Setup

The notebook is self-contained; it needs an $E(z)$ as an expression in `z` and nothing else.
That is exactly the form GR-02 returns in `hz["Branch"]`, so models carry across by copy and
paste.
::*)

(*::code::*)
ClearAll["Global`*"];

$sigma8Default = 0.811;      (* Planck 2018 TT,TE,EE+lowE+lensing *)
$zInitDefault  = 100;        (* deep enough in matter domination that delta ~ a *)

(*::md::
## 2. The solver

The initial conditions are set deep in matter domination, where the growing mode is
$\delta \propto a$, so $\delta(z_i) = 1/(1+z_i)$ and $\delta'(z_i) = -1/(1+z_i)^2$. The overall
normalisation cancels in $D$, $f$ and $f\sigma_8$, so only the shape matters.

**The guard is the important part.** A modified $G_{\text{eff}}$ can develop a pole inside the
integration range — $f(\mathcal{T}) = \mathcal{T} + \beta\mathcal{T}^2$ does exactly that, as
section 8 shows. Without a check, `NDSolve` walks into the singularity, gives up, and then
silently extrapolates the truncated interpolation back to $z = 0$, returning a confident and
completely wrong $f\sigma_8$. So the solver scans $G_{\text{eff}}$ first and refuses the job,
naming the redshift where it goes bad.
::*)

(*::code::*)
ClearAll[GrowthSolution, GrowthFactor, GrowthRate, FSigma8];

Options[GrowthSolution] = {
  "OmegaM"    -> 0.315,
  "Geff"      -> 1,
  "Sigma8"    -> $sigma8Default,
  "zInitial"  -> $zInitDefault,
  "Label"     -> "model"};

GrowthSolution[eExpr_, OptionsPattern[]] :=
 Module[{om0 = OptionValue["OmegaM"], geff = OptionValue["Geff"],
         s8 = OptionValue["Sigma8"], zi = OptionValue["zInitial"],
         probe, bad, del, sol, dfun},

  (* scan Geff and E before integrating, so a pole becomes a message not a wrong number *)
  probe = Table[{zz, N[geff /. z -> zz], N[eExpr /. z -> zz]}, {zz, 0, zi, zi/400.}];
  bad = SelectFirst[probe,
    ! (NumericQ[#[[2]]] && NumericQ[#[[3]]] && #[[2]] > 0 && #[[3]] > 0) &];
  If[! MissingQ[bad],
   Message[GrowthSolution::pole, bad[[1]]]; Return[$Failed]];

  sol = Quiet@NDSolve[{
     del''[z] + (D[eExpr, z]/eExpr - 1/(1 + z)) del'[z]
       - (3/2) om0 (1 + z) geff del[z]/eExpr^2 == 0,
     del[zi] == 1/(1 + zi), del'[zi] == -1/(1 + zi)^2},
    del, {z, 0, zi}, MaxSteps -> 10^6];
  If[sol === {} || Head[sol] =!= List, Message[GrowthSolution::nosol]; Return[$Failed]];
  dfun = del /. First[sol];

  (* the solver must have reached z = 0, not stopped at a singularity on the way *)
  If[First[dfun["Domain"][[1]]] > 10^-6,
   Message[GrowthSolution::short, First[dfun["Domain"][[1]]]]; Return[$Failed]];

  <|"Delta"        -> dfun,
    "GrowthFactor" -> Function[zz, dfun[zz]/dfun[0]],
    "GrowthRate"   -> Function[zz, -(1 + zz) dfun'[zz]/dfun[zz]],
    "FSigma8"      -> Function[zz, -s8 (1 + zz) dfun'[zz]/dfun[0]],
    "OmegaMOf"     -> Function[zz, om0 (1 + zz)^3/(eExpr /. z -> zz)^2],
    "E"            -> eExpr,
    "OmegaM"       -> om0,
    "Sigma8"       -> s8,
    "Label"        -> OptionValue["Label"]|>];

GrowthSolution::pole =
  "Geff or E is non-positive or non-numeric by z = `1` (first bad point on the scan grid, so \
the singularity itself sits at or just below this). Integrating through it would return a \
confidently wrong answer, so nothing was computed.";
GrowthSolution::nosol = "NDSolve returned no solution.";
GrowthSolution::short =
  "The integration only reached z = `1`, not z = 0. Suspect a singularity in Geff.";

(*::md::
## 3. Library of backgrounds and $G_{\text{eff}}$

Backgrounds are expressions in `z`; $G_{\text{eff}}/G$ likewise. They are chosen
independently, which is the whole point: you can hold the background fixed and vary only the
gravity.
::*)

(*::code::*)
(* ============================================================================
   LIBRARY -- copy a background and a Geff into the cell below.
   ============================================================================

   -- Backgrounds E(z), all expressions in z --------------------------------

   flat LCDM            eBg = Sqrt[om (1 + z)^3 + (1 - om)];
   with radiation       eBg = Sqrt[orad (1 + z)^4 + om (1 + z)^3 + (1 - om - orad)];
   Einstein-de Sitter   eBg = (1 + z)^(3/2);
   wCDM                 eBg = Sqrt[om (1 + z)^3 + (1 - om) (1 + z)^(3 (1 + w))];
   CPL                  eBg = Sqrt[om (1 + z)^3
                                 + (1 - om) (1 + z)^(3 (1 + w0 + wa)) Exp[-3 wa z/(1 + z)]];
   from GR-02           eBg = <paste hz["Branch"] here>;

   -- Geff/G, quasi-static and sub-horizon ----------------------------------
   (these are quoted from the literature, not derived in this notebook)

   general relativity   geff = 1;

   phenomenological mu, the form used in DESI and Euclid forecasts:
                        geff = 1 + mu0 (1 - om)/(om (1 + z)^3 + (1 - om));

   f(T) = T + b T^2, with T = -6 H0^2 E^2 so f_T = 1 - 12 bt E^2:
                        geff = 1/(1 - 12 bt eBg^2);      (* bt = b H0^2 *)

   f(Q) = Q + b Q^2, same structure with Q = -6 H0^2 E^2:
                        geff = 1/(1 - 12 bt eBg^2);

   f(R), large-scale limit (k -> 0), with fR = df/dR evaluated on the background:
                        geff = 1/fR;

   f(R), small-scale limit (k -> infinity), the famous 4/3 enhancement:
                        geff = 4/(3 fR);

   f(R), general scale dependence, x = (k^2/a^2)(fRR/fR):
                        geff = (1/fR) (1 + 4 x)/(1 + 3 x);

   ============================================================================ *)

(*::md::
## 4. Your model
::*)

(*::code::*)
ClearAll[om, mu0, eBg, geff, gr, mg];

om   = 0.315;
eBg  = Sqrt[om (1 + z)^3 + (1 - om)];

gr = GrowthSolution[eBg, "OmegaM" -> om, "Geff" -> 1,
   "Label" -> "general relativity"];

Dataset @ <|
  "label"        -> gr["Label"],
  "delta(0)"     -> gr["Delta"][0],
  "f(0)"         -> gr["GrowthRate"][0],
  "f sigma8 (0)" -> gr["FSigma8"][0],
  "D(1)"         -> gr["GrowthFactor"][1],
  "D(2)"         -> gr["GrowthFactor"][2]|>

(*::md::
## 5. Growth factor, growth rate, $f\sigma_8$

$D(z)$ falls monotonically into the past because structure has had less time to assemble.
$f(z)$ rises towards 1 going back, because at high redshift matter dominates completely and
the Einstein-de Sitter answer $\delta \propto a$ takes over exactly. $f\sigma_8$ has a
maximum near $z \approx 0.5$: at higher redshift $\sigma_8(z)$ is still small, at lower
redshift $\Lambda$ has begun shutting growth down.
::*)

(*::code::*)
Row[{
  Plot[Evaluate[{gr["GrowthFactor"][z], gr["GrowthRate"][z]}], {z, 0, 3},
   PlotLegends -> {"D(z)", "f(z)"}, AxesLabel -> {"z", None},
   PlotRange -> All, ImageSize -> 340],
  Plot[Evaluate[gr["FSigma8"][z]], {z, 0, 3},
   AxesLabel -> {"z", "f\[Sigma]8"}, PlotRange -> All, ImageSize -> 340,
   PlotLabel -> "the redshift-survey observable"]}]

(*::md::
## 6. The growth index

For a wide class of dark-energy models the growth rate is captured remarkably well by

$$
f(z) \simeq \Omega_m(z)^{\gamma}, \qquad \gamma \approx 0.55 \ \text{in } \Lambda\text{CDM}.
$$

Modified gravity shifts $\gamma$ — that is why it is quoted as a one-number summary of the
growth history. The cell below reads $\gamma(z) = \ln f/\ln\Omega_m$ straight off the
solution rather than assuming it.
::*)

(*::code::*)
Module[{gamOf},
 gamOf[zz_?NumericQ] := Log[gr["GrowthRate"][zz]]/Log[gr["OmegaMOf"][zz]];
 Column[{
   Row[{"gamma(0) = ", gamOf[0], "     gamma(1) = ", gamOf[1]}],
   Plot[{gamOf[z], 0.55}, {z, 0, 3}, PlotLegends -> {"measured", "0.55"},
    AxesLabel -> {"z", "\[Gamma]"}, PlotRange -> All, ImageSize -> 400]}]]

(*::md::
## 7. Same background, different gravity

This is the figure the whole series has been building towards. Both models below share the
*identical* $E(z)$, so every distance in GR-03 — $D_L$, $D_A$, $\mu(z)$, the BAO ruler — is
exactly the same for the two. Supernovae and BAO cannot tell them apart at any precision.

Only $G_{\text{eff}}$ differs, through the standard phenomenological
$\mu(z) = 1 + \mu_0\,\Omega_\Lambda(z)/\Omega_{\Lambda,0}$, and the growth histories separate
by several percent in $f\sigma_8$ — comfortably inside the reach of a redshift survey.
That gap is the entire observational case for measuring growth.
::*)

(*::code::*)
Module[{mu0v = 0.3},
 mg = GrowthSolution[eBg, "OmegaM" -> om,
    "Geff" -> 1 + mu0v (1 - om)/(om (1 + z)^3 + (1 - om)),
    "Label" -> "mu0 = 0.3"];
 Row[{
   Plot[Evaluate[{gr["FSigma8"][z], mg["FSigma8"][z]}], {z, 0, 3},
    PlotLegends -> {"GR", "\[Mu]0 = 0.3"}, AxesLabel -> {"z", "f\[Sigma]8"},
    PlotRange -> All, ImageSize -> 340],
   Plot[Evaluate[mg["FSigma8"][z]/gr["FSigma8"][z] - 1], {z, 0, 3},
    AxesLabel -> {"z", "fractional shift"}, PlotRange -> All, ImageSize -> 340,
    PlotLabel -> "identical distances, different growth"]}]]

(*::md::
## 8. When $G_{\text{eff}}$ has a pole

Take $f(\mathcal{T}) = \mathcal{T} + \beta\mathcal{T}^2$ from GR-02. With
$\mathcal{T} = -6H^2 = -6H_0^2E^2$ the derivative is
$f_{\mathcal{T}} = 1 - 12\tilde\beta E^2$ where $\tilde\beta = \beta H_0^2$, and the
quasi-static result $G_{\text{eff}}/G = 1/f_{\mathcal{T}}$ therefore blows up wherever
$E^2 = 1/(12\tilde\beta)$.

That $E^2$ scaling is brutal. At $z_i = 100$ a $\Lambda$CDM background already has
$E \approx 561$, so avoiding a pole all the way back to matter domination needs
$\tilde\beta \lesssim 3\times10^{-7}$ — and a coupling that small changes $f\sigma_8$ today by
about a part in a million. The honest conclusion is that a quadratic $f(\mathcal{T})$
correction cannot both stay regular through the matter era and leave a visible growth
signature, at least at the level of this scale-independent $G_{\text{eff}}$.

The cell below shows the pole redshift as a function of $\tilde\beta$, and then asks the
solver for a case that has one. It is meant to fail, and the message it prints instead of a
number is the point.
::*)

(*::code::*)
Module[{zpole, bt = 0.01},
 (* E^2 = 1/(12 bt)  =>  om (1+z)^3 + (1-om) = 1/(12 bt) *)
 zpole[b_] := ((1/(12 b) - (1 - om))/om)^(1/3) - 1;
 Column[{
   Grid[Prepend[
     {#, zpole[#]} & /@ {0.1, 0.03, 0.01, 0.003, 0.001, 10^-5},
     {Style["beta H0^2", Bold], Style["pole at z =", Bold]}],
    Alignment -> Left, Spacings -> {2, 0.7}],
   Style["asking the solver for beta H0^2 = 0.01, which has its pole at z = "
     <> ToString[zpole[bt]] <> ":", Italic],
   GrowthSolution[eBg, "OmegaM" -> om, "Geff" -> 1/(1 - 12 bt eBg^2)]}]]

(*::md::
## 9. Checks

Seven tests, each against an identity or a number that does not come from this notebook:

1. Einstein-de Sitter has $f = 1$ exactly and $\delta \propto a$, because in a
   matter-only universe the growing mode is the scale factor itself.
2. $\Lambda$CDM growth matches the exact quadrature for the growing mode,
   $D(z) \propto E(z)\int_z^\infty (1+z')\,E(z')^{-3}dz'$.
3. $f(0) \simeq \Omega_{m0}^{0.55}$, and the growth-index fit holds to better than 2 percent
   out to $z = 2$.
4. $f\sigma_8(0) \approx 0.43$ for Planck parameters.
5. $G_{\text{eff}} = 1$ written as an expression reproduces the plain GR run exactly.
6. A stronger $G_{\text{eff}}$ grows more structure.
7. The pole guard fires on a $G_{\text{eff}}$ that is singular in range, instead of returning
   a wrong number.
::*)

(*::code::*)
Module[{eds, lcdm, dExact, ratios, devs, omz, strong, guarded},
 eds  = GrowthSolution[(1 + z)^(3/2), "OmegaM" -> 1, "Geff" -> 1];
 lcdm = gr;

 dExact[zz_?NumericQ] :=
  (eBg /. z -> zz) NIntegrate[(1 + zp)/(eBg /. z -> zp)^3, {zp, zz, Infinity}];
 ratios = Table[
   (lcdm["GrowthFactor"][zz])/(dExact[zz]/dExact[0]), {zz, 0.5, 3, 0.5}];

 omz[zz_] := om (1 + zz)^3/(eBg /. z -> zz)^2;
 devs = Table[
   Abs[lcdm["GrowthRate"][zz] - omz[zz]^0.55]/lcdm["GrowthRate"][zz], {zz, 0, 2, 0.25}];

 strong  = GrowthSolution[eBg, "OmegaM" -> om, "Geff" -> 1.2];
 guarded = Quiet@GrowthSolution[eBg, "OmegaM" -> om, "Geff" -> 1/(1 - 12 (0.01) eBg^2)];

 Dataset @ {
  <|"check" -> "Einstein-de Sitter has f = 1 and delta proportional to a",
    "ok" -> (Abs[eds["GrowthRate"][0] - 1] < 10^-3 &&
             Abs[eds["GrowthRate"][2] - 1] < 10^-3 &&
             Abs[eds["GrowthFactor"][1] - 1/2] < 10^-3)|>,
  <|"check" -> "LCDM matches the exact growing-mode quadrature",
    "ok" -> (Max[Abs[ratios - 1]] < 0.002)|>,
  <|"check" -> "f(0) agrees with OmegaM^0.55",
    "ok" -> (Abs[lcdm["GrowthRate"][0] - om^0.55] < 0.02 om^0.55)|>,
  <|"check" -> "growth index fit holds to 2 percent out to z = 2",
    "ok" -> (Max[devs] < 0.02)|>,
  <|"check" -> "f sigma8 (0) is about 0.43 for Planck parameters",
    "ok" -> (Abs[lcdm["FSigma8"][0] - 0.43] < 0.02)|>,
  <|"check" -> "Geff = 1 as an expression reproduces plain GR",
    "ok" -> (Abs[GrowthSolution[eBg, "OmegaM" -> om, "Geff" -> 1 + 0 z]["FSigma8"][0]
              - lcdm["FSigma8"][0]] < 10^-8)|>,
  <|"check" -> "stronger gravity grows more structure",
    "ok" -> (strong["FSigma8"][0] > lcdm["FSigma8"][0])|>,
  <|"check" -> "the pole guard refuses a singular Geff",
    "ok" -> (guarded === $Failed)|>
 }]
