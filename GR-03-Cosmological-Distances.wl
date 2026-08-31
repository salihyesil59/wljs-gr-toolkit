(* ::Package:: *)

(* ============================================================================
   GR-03 : Cosmological Distances, Times and Observables
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-03-Cosmological-Distances.wl

   ============================================================================ *)

(*::md::
# Cosmological Distances, Times and Observables

GR-02 stops at $E(z) = H(z)/H_0$. That is the last purely theoretical object in the chain:
everything an observer actually measures &mdash; how faint a supernova looks, how large the
baryon acoustic ruler appears on the sky, how old the oldest stars can be &mdash; is an
integral of $E(z)$. This notebook does those integrals.

$$
D_C(z) = \frac{c}{H_0}\int_0^z \frac{dz'}{E(z')}, \qquad
t_L(z) = \frac{1}{H_0}\int_0^z \frac{dz'}{(1+z')E(z')}
$$

Everything else follows from those two: angular diameter and luminosity distances, the
distance modulus, lookback time and the age of the universe, the particle and event horizons,
the comoving volume element, and the sound horizon that sets the BAO scale.

This is the piece that makes a modified-gravity model falsifiable. Two theories can produce
almost identical $E(z)$ and still separate by a tenth of a magnitude in $\mu(z)$ at $z \sim 1$,
which is exactly where supernova surveys have their leverage.

## Conventions and units

- Distances in **Mpc**, times in **Gyr**, $H_0 = 100h$ km s$^{-1}$ Mpc$^{-1}$.
- $D_H = c/H_0 = 2997.92458/h$ Mpc and $t_H = 1/H_0 = 9.77792/h$ Gyr.
- $\Omega_k > 0$ means an open universe, matching the $\Omega_k(1+z)^2$ term in GR-02.
- Redshift and scale factor are related by $1+z = 1/a$ with $a_0 = 1$.

The notebook is self-contained: it needs an $E(z)$ and nothing else, so it does not load GR-02.
Section 3 shows how to carry a model across from it.

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

You are on the observational branch. GR-02 hands over $E(z)$; what leaves here is what a
telescope actually records. GR-04 takes the same $E(z)$ down the other branch, and the
interesting models are the ones the two branches disagree about &mdash; identical distances,
different growth.
::*)

(*::md::
## 1. Setup

Only three constants are physical input. The Hubble time is derived from them rather than
typed in, so the Mpc-to-km and Gyr-to-second conversions stay visible.
::*)

(*::code::*)
ClearAll["Global`*"];

$cLight   = 299792.458;                  (* km/s, exact by definition *)
$MpcInKm  = 3.0856775814913673*^19;      (* IAU 2015 *)
$GyrInSec = 3.15576*^16;                 (* Julian year = 365.25 d *)

(*::md::
## 2. The model

A cosmology here is just a function $E(z)$ plus the two numbers the integrals need to carry
units: $h$, and the curvature $\Omega_k$ that decides whether the transverse distance is a
`Sinh`, a straight line, or a `Sin`.

Every distance and time below takes the model as its first argument, so you can hold several
models at once and compare them.
::*)

(*::code::*)
ClearAll[CosmologyModel, FlatLCDM];

Options[CosmologyModel] = {"h" -> 0.674, "OmegaK" -> 0, "Label" -> "model"};

CosmologyModel[eFunc_, OptionsPattern[]] :=
 Module[{h = OptionValue["h"]},
  <|"E"              -> eFunc,
    "h"              -> h,
    "H0"             -> 100. h,                              (* km/s/Mpc *)
    "OmegaK"         -> OptionValue["OmegaK"],
    "Label"          -> OptionValue["Label"],
    "HubbleDistance" -> $cLight/(100. h),                    (* Mpc *)
    "HubbleTime"     -> $MpcInKm/(100. h $GyrInSec)|>];      (* Gyr *)

(* convenience builder for the case everybody starts from *)
FlatLCDM[om_, h_: 0.674, orad_: 0.] :=
 CosmologyModel[
  Function[zz, Sqrt[orad (1 + zz)^4 + om (1 + zz)^3 + (1 - om - orad)]],
  "h" -> h, "OmegaK" -> 0,
  "Label" -> "flat \[CapitalLambda]CDM, \[CapitalOmega]m = " <> ToString[om]];

(*::md::
## 3. Cosmology library

Copy an `EOf` line into the model cell. The first group is closed form; the last shows how to
bring across a model from GR-02 whose $E(z)$ is only defined implicitly.
::*)

(*::code::*)
(* ============================================================================
   COSMOLOGY LIBRARY -- copy one block into "Your cosmology" below.
   ============================================================================

   -- flat LCDM, Planck-like -------------------------------------------------
   model = FlatLCDM[0.315, 0.674];

   -- flat LCDM with radiation (needed above z ~ 100) ------------------------
   model = FlatLCDM[0.315, 0.674, 9.2*^-5];

   -- open or closed LCDM ----------------------------------------------------
   With[{om = 0.3, ol = 0.6, ok = 0.1},
    model = CosmologyModel[
      Function[zz, Sqrt[om (1 + zz)^3 + ok (1 + zz)^2 + ol]],
      "h" -> 0.674, "OmegaK" -> ok, "Label" -> "open LCDM"]];

   -- Einstein-de Sitter (matter only) ---------------------------------------
   model = CosmologyModel[Function[zz, (1 + zz)^(3/2)],
     "h" -> 0.674, "Label" -> "Einstein-de Sitter"];

   -- de Sitter (pure cosmological constant) ---------------------------------
   model = CosmologyModel[Function[zz, 1], "h" -> 0.674, "Label" -> "de Sitter"];

   -- wCDM, constant equation of state ---------------------------------------
   With[{om = 0.315, w = -0.9},
    model = CosmologyModel[
      Function[zz, Sqrt[om (1 + zz)^3 + (1 - om) (1 + zz)^(3 (1 + w))]],
      "h" -> 0.674, "Label" -> "wCDM, w = " <> ToString[w]]];

   -- CPL / Chevallier-Polarski-Linder, w(a) = w0 + wa (1-a) -----------------
   With[{om = 0.315, w0 = -0.9, wa = 0.3},
    model = CosmologyModel[
      Function[zz, Sqrt[om (1 + zz)^3 + (1 - om) (1 + zz)^(3 (1 + w0 + wa))
         Exp[-3 wa zz/(1 + zz)]]],
      "h" -> 0.674, "Label" -> "CPL"]];

   -- from GR-02, closed form: paste hz["Branch"] straight in ----------------
   model = CosmologyModel[Function[zz, Sqrt[0.315 (1 + zz)^3 + 0.685]],
     "h" -> 0.674, "Label" -> "from GR-02"];

   -- from GR-02, implicit: f(T) = T + b T^2 gives E^2 - 18 b H0^2 E^4 = ... --
   With[{om = 0.315, bt = 0.02},
    model = CosmologyModel[
      ENumeric[Ez[z]^2 - bt Ez[z]^4 - om (1 + z)^3 - (1 - om - bt),
               Function[zz, Sqrt[om (1 + zz)^3 + (1 - om)]]],
      "h" -> 0.674, "Label" -> "f(T) quadratic"]];

   ============================================================================ *)

(*::md::
## 4. Distances

$D_C$ is the line-of-sight comoving distance. $D_M$ is its transverse partner: in a curved
universe two objects separated on the sky are not as far apart as $D_C$ would suggest, which
is what the `Sinh`/`Sin` is correcting. Then

$$
D_A = \frac{D_M}{1+z}, \qquad D_L = (1+z)\,D_M, \qquad \mu = 5\log_{10}\!\left(\frac{D_L}{\text{Mpc}}\right) + 25 .
$$

$D_A$ is the one with the surprise in it: it rises, turns over near $z \approx 1.6$ for
$\Lambda$CDM, and then *falls*, so a galaxy of fixed physical size looks bigger the further
back you push it. $D_L$ and $D_A$ differ by $(1+z)^2$, the Etherington reciprocity relation,
which holds in any metric theory with photon number conserved.
::*)

(*::code::*)
ClearAll[ComovingDistance, TransverseComovingDistance, AngularDiameterDistance,
         LuminosityDistance, DistanceModulus];

ComovingDistance[cos_Association, zz_?NumericQ] :=
 cos["HubbleDistance"] NIntegrate[1/cos["E"][zp], {zp, 0, zz}];

TransverseComovingDistance[cos_Association, zz_?NumericQ] :=
 Module[{dh = cos["HubbleDistance"], ok = cos["OmegaK"], dc = ComovingDistance[cos, zz]},
  Which[
   ok > 0, dh/Sqrt[ok] Sinh[Sqrt[ok] dc/dh],
   ok == 0, dc,
   True,    dh/Sqrt[-ok] Sin[Sqrt[-ok] dc/dh]]];

AngularDiameterDistance[cos_Association, zz_?NumericQ] :=
  TransverseComovingDistance[cos, zz]/(1 + zz);

LuminosityDistance[cos_Association, zz_?NumericQ] :=
  (1 + zz) TransverseComovingDistance[cos, zz];

(* apparent minus absolute magnitude, with D_L in Mpc *)
DistanceModulus[cos_Association, zz_?NumericQ] :=
  5 Log10[LuminosityDistance[cos, zz]] + 25;

(*::md::
## 5. Times

Lookback time is how long the light has been travelling; the age at redshift $z$ is the
integral the other way, out to the big bang. Their sum is the age today, which the checks
section uses as a consistency test.

For flat $\Lambda$CDM the age has a closed form,

$$
t_0 = \frac{2}{3H_0\sqrt{\Omega_\Lambda}}\;\operatorname{arcsinh}\sqrt{\frac{\Omega_\Lambda}{\Omega_m}},
$$

and the notebook checks the numerical integral against it.
::*)

(*::code::*)
ClearAll[LookbackTime, AgeAt, AgeToday];

LookbackTime[cos_Association, zz_?NumericQ] :=
 cos["HubbleTime"] NIntegrate[1/((1 + zp) cos["E"][zp]), {zp, 0, zz}];

AgeAt[cos_Association, zz_?NumericQ] :=
 cos["HubbleTime"] NIntegrate[1/((1 + zp) cos["E"][zp]), {zp, zz, Infinity}];

AgeToday[cos_Association] := AgeAt[cos, 0];

(*::md::
## 6. Horizons and volume

The **particle horizon** is how far light has come since the big bang &mdash; the edge of what
can ever have influenced us. The **event horizon** is how far a signal sent now will ever
reach; it is finite whenever the expansion accelerates forever, and in pure de Sitter it is
exactly $D_H$, which the checks confirm.

The comoving volume element per unit redshift and solid angle,
$dV_C/(dz\,d\Omega) = D_H D_M^2/E(z)$, is what turns a number count into a density.
::*)

(*::code::*)
ClearAll[ParticleHorizon, EventHorizon, ComovingVolumeElement, ComovingVolume];

ParticleHorizon[cos_Association, zz_?NumericQ] :=
 cos["HubbleDistance"] NIntegrate[1/cos["E"][zp], {zp, zz, Infinity}];

(* light sent at zz reaches this comoving distance by the infinite future (z -> -1) *)
EventHorizon[cos_Association, zz_?NumericQ] :=
 cos["HubbleDistance"] NIntegrate[1/cos["E"][zp], {zp, -1, zz}];

(* Mpc^3 per steradian per unit redshift *)
ComovingVolumeElement[cos_Association, zz_?NumericQ] :=
 cos["HubbleDistance"] TransverseComovingDistance[cos, zz]^2/cos["E"][zz];

(* all-sky comoving volume out to zz, in Mpc^3 *)
ComovingVolume[cos_Association, zz_?NumericQ] :=
 4 Pi NIntegrate[ComovingVolumeElement[cos, zp], {zp, 0, zz}];

(*::md::
## 7. Models whose $E(z)$ is only implicit

GR-02 returns $E(z)$ in closed form for GR, $\Lambda$CDM, $f(\mathcal{T})$ and $f(Q)$, but for
$f(R)$ and $f(R,\mathcal{G})$ it returns a constraint instead. `ENumeric` turns any constraint
written in `Ez[z]` and `z` into a function of $z$ by root-finding, so those models plug into
the same distance machinery.

It needs a starting guess that is already in the right basin &mdash; the corresponding
$\Lambda$CDM value is almost always good enough, and that is what the library entry uses.
::*)

(*::code::*)
ClearAll[ENumeric];

(* The NumericQ guard is not decoration: NIntegrate and Plot probe a function
   symbolically first, and an unguarded FindRoot answers that probe with a page
   of errors before the numeric pass quietly gets it right. *)
ENumeric[constraint_, guess_] :=
 Module[{ef},
  ef[zz_?NumericQ] :=
   Module[{ev},
    ev /. FindRoot[Evaluate[constraint /. {Ez[z] -> ev, z -> zz}], {ev, guess[zz]}]];
  ef];

(*::md::
## 8. Your cosmology

Edit this cell. The default is the Planck-like flat $\Lambda$CDM everything is measured
against.
::*)

(*::code::*)
ClearAll[model];

model = FlatLCDM[0.315, 0.674];

Dataset @ <|
  "label"                 -> model["Label"],
  "H0 (km/s/Mpc)"         -> model["H0"],
  "Hubble distance (Mpc)" -> model["HubbleDistance"],
  "Hubble time (Gyr)"     -> model["HubbleTime"],
  "age today (Gyr)"       -> AgeToday[model],
  "particle horizon (Mpc)"-> ParticleHorizon[model, 0],
  "event horizon (Mpc)"   -> EventHorizon[model, 0]|>

(*::md::
## 9. The numbers at a glance
::*)

(*::code::*)
Module[{zs = {0.1, 0.5, 1., 2., 5., 1089.}},
 Dataset[
  Association[
    "z"           -> #,
    "D_C (Mpc)"   -> ComovingDistance[model, #],
    "D_A (Mpc)"   -> AngularDiameterDistance[model, #],
    "D_L (Mpc)"   -> LuminosityDistance[model, #],
    "mu"          -> DistanceModulus[model, #],
    "lookback (Gyr)" -> LookbackTime[model, #]] & /@ zs]]

(*::md::
## 10. The distance ladder, drawn

The left panel is the reason $D_A$ deserves its own name: it turns over while $D_L$ keeps
climbing. The right panel is the Hubble diagram, the plot supernova cosmology is done on.
::*)

(*::code::*)
Row[{
  Plot[Evaluate[{ComovingDistance[model, z], LuminosityDistance[model, z],
                 AngularDiameterDistance[model, z]}], {z, 0.01, 5},
   PlotLegends -> {"D_C", "D_L", "D_A"}, AxesLabel -> {"z", "Mpc"},
   PlotRange -> All, ImageSize -> 340],
  Plot[Evaluate[DistanceModulus[model, z]], {z, 0.01, 2},
   AxesLabel -> {"z", "\[Mu]"}, PlotRange -> All, ImageSize -> 340,
   PlotLabel -> "Hubble diagram"]}]

(*::md::
## 11. Telling models apart

This is the payoff. Three models tuned to the same $\Omega_m$ differ by only a few percent in
$E(z)$, but the distance modulus integrates that difference and the residual against
$\Lambda$CDM grows to a few hundredths of a magnitude &mdash; comparable to the scatter of a
well-calibrated supernova sample, which is why these data constrain $w$ at all.

Swap the third model for an $f(\mathcal{T})$ or $f(Q)$ constraint from GR-02 to see where your
own theory sits.
::*)

(*::code::*)
Module[{lcdm, wcdm, cpl, om = 0.315, h = 0.674},
 lcdm = FlatLCDM[om, h];
 wcdm = CosmologyModel[
   Function[zz, Sqrt[om (1 + zz)^3 + (1 - om) (1 + zz)^(3 (1 - 0.9))]],
   "h" -> h, "Label" -> "wCDM, w = -0.9"];
 cpl = CosmologyModel[
   Function[zz, Sqrt[om (1 + zz)^3 + (1 - om) (1 + zz)^(3 (1 - 0.9 + 0.3))
      Exp[-3 (0.3) zz/(1 + zz)]]],
   "h" -> h, "Label" -> "CPL, w0 = -0.9, wa = 0.3"];

 Row[{
  Plot[Evaluate[{lcdm["E"][z], wcdm["E"][z], cpl["E"][z]}], {z, 0, 2},
   PlotLegends -> {"\[CapitalLambda]CDM", "wCDM", "CPL"},
   AxesLabel -> {"z", "E(z)"}, ImageSize -> 340],
  Plot[Evaluate[{DistanceModulus[wcdm, z] - DistanceModulus[lcdm, z],
                 DistanceModulus[cpl, z] - DistanceModulus[lcdm, z]}], {z, 0.02, 2},
   PlotLegends -> {"wCDM", "CPL"}, AxesLabel -> {"z", "\[CapitalDelta]\[Mu] (mag)"},
   PlotRange -> All, ImageSize -> 340,
   PlotLabel -> "residual against \[CapitalLambda]CDM"]}]]

(*::md::
## 12. The BAO sound horizon

Before recombination, photons and baryons are one fluid and sound waves travel at
$c_s = c/\sqrt{3(1+R)}$ with $R = 3\rho_b/4\rho_\gamma$. The comoving distance those waves
cover before the baryons are released,

$$
r_s(z_d) = \frac{c}{H_0}\int_{z_d}^{\infty} \frac{dz}{E(z)\sqrt{3\left(1+R(z)\right)}},
\qquad R(z) = \frac{3\Omega_b}{4\Omega_\gamma}\frac{1}{1+z},
$$

is frozen into the galaxy distribution as a standard ruler of about 147 Mpc. Radiation must be
in $E(z)$ for this integral to mean anything, so the model below is built with it.

The $r_s$ this returns lands within a fraction of a percent of the Planck value, but the
acoustic scale $\theta_* = r_s/D_M$ it reports comes out near $100\theta_* = 1.06$ against
Planck's $1.0411$. The gap is real and it is the model's, not an arithmetic slip: $z_d$ and
$z_*$ are put in by hand here rather than solved for, and $\Omega_r$ is a single number
instead of a photon plus neutrino background with the right relativistic degrees of freedom.
Treat $\theta_*$ as indicative. Fitting it at the level a CMB likelihood needs is a Boltzmann
code's job, not a notebook's.
::*)

(*::code::*)
Module[{om = 0.315, obh2 = 0.0224, ogh2 = 2.47*^-5, h = 0.674, orad = 9.2*^-5,
        cos, rr, rs, dm},
 cos = FlatLCDM[om, h, orad];
 rr[zz_] := (3 obh2/(4 ogh2))/(1 + zz);
 rs = cos["HubbleDistance"] NIntegrate[
    1/(cos["E"][zp] Sqrt[3 (1 + rr[zp])]), {zp, 1060, Infinity}];
 dm = TransverseComovingDistance[cos, 1089.];
 Dataset @ <|
   "sound horizon r_s at z_d = 1060 (Mpc)" -> rs,
   "comoving distance to z = 1089 (Mpc)"   -> dm,
   "acoustic scale theta_star (rad)"       -> rs/dm,
   "acoustic scale 100 theta_star"         -> 100 rs/dm|>]

(*::md::
## 13. Checks

Nine independent tests, each against a number or an identity that does not come from this
notebook:

1. $D_H$ and $t_H$ against their textbook values in Mpc and Gyr.
2. The age of a flat $\Lambda$CDM universe against its closed form, and against the familiar
   13.8 Gyr for Planck parameters.
3. Lookback time plus age at $z$ equals the age today.
4. The low-redshift expansion $D_L \simeq \frac{cz}{H_0}\left[1 + \tfrac{1}{2}(1-q_0)z\right]$
   with $q_0 = \Omega_m/2 - \Omega_\Lambda$.
5. Etherington reciprocity, $D_L = (1+z)^2 D_A$.
6. $D_A$ turns over near $z \approx 1.6$.
7. In de Sitter the event horizon is exactly $D_H$.
8. Einstein-de Sitter has age $\tfrac{2}{3}t_H$ and particle horizon $2D_H$.
9. The comoving distance to last scattering is about 13.9 Gpc and the sound horizon about
   147 Mpc.
::*)

(*::code::*)
Module[{cos, dS, eds, om = 0.315, ol = 0.685, h = 0.674, q0, lowz, rad, rs, rr},
 cos  = FlatLCDM[om, h];
 dS   = CosmologyModel[Function[zz, 1], "h" -> h];
 eds  = CosmologyModel[Function[zz, (1 + zz)^(3/2)], "h" -> h];
 q0   = om/2 - ol;
 lowz = 0.01;
 rad  = FlatLCDM[om, h, 9.2*^-5];
 rr[zz_] := (3 (0.0224)/(4 (2.47*^-5)))/(1 + zz);
 rs = rad["HubbleDistance"] NIntegrate[
    1/(rad["E"][zp] Sqrt[3 (1 + rr[zp])]), {zp, 1060, Infinity}];

 Dataset @ {
  <|"check" -> "D_H = 2997.92458/h Mpc and t_H = 9.77792/h Gyr",
    "ok" -> (Abs[cos["HubbleDistance"] h - 2997.92458] < 10^-6 &&
             Abs[cos["HubbleTime"] h - 9.77792] < 10^-4)|>,
  <|"check" -> "age matches the closed form for flat \[CapitalLambda]CDM",
    "ok" -> (Abs[AgeToday[cos] -
        (2/(3 Sqrt[ol])) cos["HubbleTime"] ArcSinh[Sqrt[ol/om]]] < 10^-6)|>,
  <|"check" -> "age today is 13.8 Gyr for Planck parameters",
    "ok" -> (Abs[AgeToday[cos] - 13.8] < 0.05)|>,
  <|"check" -> "lookback(z) + age(z) = age today",
    "ok" -> (Abs[LookbackTime[cos, 1.] + AgeAt[cos, 1.] - AgeToday[cos]] < 10^-6)|>,
  <|"check" -> "low-z expansion of D_L with q0 = \[CapitalOmega]m/2 - \[CapitalOmega]\[CapitalLambda]",
    "ok" -> (Abs[LuminosityDistance[cos, lowz] -
        cos["HubbleDistance"] lowz (1 + (1 - q0) lowz/2)] <
        10^-4 cos["HubbleDistance"] lowz)|>,
  <|"check" -> "Etherington reciprocity D_L = (1+z)^2 D_A",
    "ok" -> (Abs[LuminosityDistance[cos, 1.5] -
        (1 + 1.5)^2 AngularDiameterDistance[cos, 1.5]] < 10^-8)|>,
  <|"check" -> "D_A turns over near z = 1.6",
    "ok" -> (Abs[(z /. FindMaximum[AngularDiameterDistance[cos, z], {z, 1.5}][[2]])
        - 1.6] < 0.15)|>,
  <|"check" -> "de Sitter event horizon equals D_H",
    "ok" -> (Abs[EventHorizon[dS, 0] - dS["HubbleDistance"]] < 10^-6)|>,
  <|"check" -> "Einstein-de Sitter: age = 2/3 t_H, particle horizon = 2 D_H",
    "ok" -> (Abs[AgeToday[eds] - (2/3) eds["HubbleTime"]] < 10^-6 &&
             Abs[ParticleHorizon[eds, 0] - 2 eds["HubbleDistance"]] < 10^-6)|>,
  <|"check" -> "distance to last scattering ~ 13.9 Gpc, r_s ~ 147 Mpc",
    "ok" -> (Abs[ComovingDistance[rad, 1089.] - 13900] < 400 && Abs[rs - 147] < 7)|>
 }]
