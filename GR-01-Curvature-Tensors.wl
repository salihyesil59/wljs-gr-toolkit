(* ::Package:: *)

(* ============================================================================
   GR-01 : Curvature Tensors from a Metric
   WLJS Notebook source.  Build the notebook with

       wolframscript -file wl2wln.wls GR-01-Curvature-Tensors.wl

   ============================================================================ *)

(*::md::
# Curvature Tensors in General Relativity

Supply a coordinate chart and a metric tensor $g_{\mu\nu}$; this notebook returns the
**Christoffel symbols**, the **Riemann tensor**, the **Ricci tensor**, the **Ricci scalar**
and the **Einstein tensor** &mdash; listing only the non-zero, algebraically independent components.

The curvature itself is computed by the Wolfram Function Repository resources
`MetricTensor`, `ChristoffelSymbols`, `RiemannTensor`, `RicciTensor` and `EinsteinTensor`.
They are downloaded once and cached locally, so the first evaluation of the setup cell
needs a network connection and later ones do not. What this notebook adds on top is the
part those resources do not do: pulling out the non-zero, independent components and
labelling them by coordinate name.

## Conventions

Signature $(-,+,+,+)$, geometrized units $G = c = 1$, and the Levi-Civita connection
(torsion-free and metric-compatible).

$$
\Gamma^{\lambda}{}_{\mu\nu} \;=\; \tfrac{1}{2}\, g^{\lambda\sigma}\left(\partial_{\mu} g_{\sigma\nu} + \partial_{\nu} g_{\sigma\mu} - \partial_{\sigma} g_{\mu\nu}\right)
$$

$$
R^{\rho}{}_{\sigma\mu\nu} \;=\; \partial_{\mu}\Gamma^{\rho}{}_{\nu\sigma} - \partial_{\nu}\Gamma^{\rho}{}_{\mu\sigma} + \Gamma^{\rho}{}_{\mu\lambda}\Gamma^{\lambda}{}_{\nu\sigma} - \Gamma^{\rho}{}_{\nu\lambda}\Gamma^{\lambda}{}_{\mu\sigma}
$$

$$
R_{\mu\nu} = R^{\lambda}{}_{\mu\lambda\nu}, \qquad R = g^{\mu\nu} R_{\mu\nu}, \qquad G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}\, g_{\mu\nu} R
$$

With these signs a round sphere has $R > 0$, and the vacuum field equations read $G_{\mu\nu} = 0$.

## How to use it

1. Evaluate **Setup** and **Engine** once.
2. Pick a metric from the **Metric library** (one long comment, ready to copy) and paste it into **Your metric**.
3. Evaluate **Run**, then read off the results.

## The series

| | |
|---|---|
| **GR-01** | a metric in, curvature tensors out |
| **GR-02** | an action in, field equations, Friedmann and $E(z)$ out |
| **GR-03** | $E(z)$ in, distances, times and the BAO ruler out |
| **GR-04** | $E(z)$ and $G_{\text{eff}}$ in, growth and $f\sigma_8$ out |
| **GR-05** | perturbed field equations in, $G_{\text{eff}}$, slip and lensing out |

This is where the chain starts, and it is the only notebook that asks for nothing but geometry.
Everything downstream leans on the curvature engine set up here: GR-02 runs it on the FLRW
ansatz to get the Friedmann equations, and GR-05 runs it on a *perturbed* metric to get the
effective gravitational coupling.
::*)

(*::md::
## 1. Setup

The first evaluation pulls the curvature resources from the Function Repository; after
that they come from the local cache. The simplifier is applied on top of the already
reduced output of those resources, and only when you supply assumptions, so leaving
assumptions at True costs nothing.
::*)

(*::code::*)
ClearAll["Global`*"];

$grMetric      = ResourceFunction["MetricTensor"];
$grChristoffel = ResourceFunction["ChristoffelSymbols"];
$grRiemann     = ResourceFunction["RiemannTensor"];
$grRicci       = ResourceFunction["RicciTensor"];
$grEinstein    = ResourceFunction["EinsteinTensor"];

grSimplify = Simplify;   (* Simplify | FullSimplify | Identity *)

(*::md::
## 2. Engine

`GRTensors[metric, coords]` wraps the repository resources into a single Association:

| key | meaning |
|---|---|
| `"Christoffel"` | $\Gamma^{\lambda}{}_{\mu\nu}$, indexed `[[\[Lambda], \[Mu], \[Nu]]]` |
| `"Riemann"` | $R^{\rho}{}_{\sigma\mu\nu}$ |
| `"RiemannLower"` | $R_{\rho\sigma\mu\nu}$ |
| `"Ricci"`, `"RicciScalar"`, `"Einstein"` | as named |
| `"Kretschmann"` | $R_{\rho\sigma\mu\nu}R^{\rho\sigma\mu\nu}$ |
| `"Resources"` | the underlying resource objects, for properties such as `"BianchiIdentities"` or `"ChernPontryaginScalar"` |

Nothing in this section needs editing.
::*)

(*::code::*)
ClearAll[GRTensors];

Options[GRTensors] = {"Simplify" :> grSimplify, "Assumptions" -> True};

GRTensors[metricIn_?MatrixQ, coordsIn_List, OptionsPattern[]] :=
 Module[{n, simp, mt, chr, rie, ric, ein, g, ginv, rs},

  n = Length[coordsIn];
  If[Dimensions[metricIn] =!= {n, n},
    Message[GRTensors::dim, Dimensions[metricIn], n]; Return[$Failed]
  ];

  (* The repository resources already return simplified ("Reduced") output, so a second
     pass is worth its cost only when there are assumptions to feed it. *)
  With[{f = OptionValue["Simplify"], as = OptionValue["Assumptions"]},
    simp = Which[f === Identity, Identity, as === True, Identity, True, Function[expr, f[expr, as]]]
  ];

  mt = $grMetric[metricIn, coordsIn];
  If[simp[Det[metricIn]] === 0, Message[GRTensors::deg]; Return[$Failed]];

  chr = $grChristoffel[mt];
  rie = $grRiemann[mt];
  ric = $grRicci[mt];
  ein = $grEinstein[mt];

  g    = Map[simp, mt["ReducedMatrixRepresentation"], {2}];
  ginv = Map[simp, mt["InverseMetricTensor"]["ReducedMatrixRepresentation"], {2}];
  rs   = simp[Tr[ginv . ric["ReducedMatrixRepresentation"]]];

  <|
   "Coordinates"   -> coordsIn,
   "Dimension"     -> n,
   "Metric"        -> g,
   "MetricInverse" -> ginv,
   "Christoffel"   -> Map[simp, chr["ReducedTensorRepresentation"], {3}],
   "Riemann"       -> Map[simp, rie["ReducedTensorRepresentation"], {4}],
   "RiemannLower"  -> Map[simp, rie["CovariantRiemannTensor"]["ReducedTensorRepresentation"], {4}],
   "Ricci"         -> Map[simp, ric["ReducedMatrixRepresentation"], {2}],
   "RicciScalar"   -> rs,
   "Einstein"      -> Map[simp, ein["ReducedMatrixRepresentation"], {2}],
   "Kretschmann"   -> simp[rie["ReducedKretschmannScalar"]],
   "Resources"     -> <|"Metric" -> mt, "Christoffel" -> chr, "Riemann" -> rie,
                        "Ricci" -> ric, "Einstein" -> ein|>,
   "Simplifier"    -> simp
  |>
 ];

GRTensors::dim = "Metric has dimensions `1`, expected `2` x `2`.";
GRTensors::deg = "The metric is degenerate (vanishing determinant).";

(*::code::*)
(* ---- presentation helpers ------------------------------------------------ *)

ClearAll[grGrid, grIndex, ShowMetric, ShowChristoffel, ShowRiemann,
         ShowRicci, ShowRicciScalar, ShowEinstein, ShowSummary];

(* Every loop index below is declared in its Module. Table and Sum scope iterators
   dynamically, so a bare Table[..., {a, n}] would overwrite a symbol named "a" living
   inside the tensor being displayed -- and FLRW calls its scale factor a[t], while Kerr
   calls its spin parameter a. Declaring them makes them lexically renamed locals. *)

grIndex[x_List, is__Integer] := Row[Riffle[x[[{is}]], "\[VeryThinSpace]"]];

grGrid[rows_List] :=
 If[rows === {},
  Style["all components vanish", Italic, GrayLevel[0.5]],
  Grid[{#1, "=", #2} & @@@ rows,
   Alignment -> {{Right, Center, Left}}, Spacings -> {1.2, 0.9}]
 ];

ShowMetric[gr_Association] :=
  MatrixForm[gr["Metric"], TableHeadings -> {gr["Coordinates"], gr["Coordinates"]}];

(* \[CapitalGamma]^a_bc is symmetric in (b,c): keep b <= c *)
ShowChristoffel[gr_Association] :=
 Module[{n = gr["Dimension"], x = gr["Coordinates"], c = gr["Christoffel"], a, b, d},
  grGrid @ Flatten[
    Table[
     If[c[[a, b, d]] === 0, Nothing,
      {Subsuperscript["\[CapitalGamma]", grIndex[x, b, d], x[[a]]], c[[a, b, d]]}],
     {a, n}, {b, n}, {d, b, n}],
    2]
 ];

(* R_abcd is antisymmetric in (ab) and in (cd), and symmetric under (ab) <-> (cd).
   Independent set: index pairs P = (a<b), Q = (c<d) with position[P] <= position[Q]. *)
ShowRiemann[gr_Association] :=
 Module[{n = gr["Dimension"], x = gr["Coordinates"], r = gr["RiemannLower"], pairs,
         i, j, a, b, c, d},
  pairs = Subsets[Range[n], {2}];
  grGrid @ Flatten[
    Table[
     With[{a = pairs[[i, 1]], b = pairs[[i, 2]], c = pairs[[j, 1]], d = pairs[[j, 2]]},
      If[r[[a, b, c, d]] === 0, Nothing,
       {Subscript["R", grIndex[x, a, b, c, d]], r[[a, b, c, d]]}]],
     {i, Length[pairs]}, {j, i, Length[pairs]}],
    1]
 ];

(* R_ab is symmetric *)
ShowRicci[gr_Association] :=
 Module[{n = gr["Dimension"], x = gr["Coordinates"], r = gr["Ricci"], a, b},
  grGrid @ Flatten[
    Table[
     If[r[[a, b]] === 0, Nothing, {Subscript["R", grIndex[x, a, b]], r[[a, b]]}],
     {a, n}, {b, a, n}],
    1]
 ];

ShowRicciScalar[gr_Association] := Row[{Style["R", Italic], " = ", gr["RicciScalar"]}];

ShowEinstein[gr_Association] :=
 Module[{n = gr["Dimension"], x = gr["Coordinates"], e = gr["Einstein"], a, b},
  grGrid @ Flatten[
    Table[
     If[e[[a, b]] === 0, Nothing, {Subscript["G", grIndex[x, a, b]], e[[a, b]]}],
     {a, n}, {b, a, n}],
    1]
 ];

ShowSummary[gr_Association] :=
 Module[{count = Function[tens, Count[Flatten[{tens}], q_ /; q =!= 0]]},
  Dataset @ {
   <|"tensor" -> "Christoffel  \[CapitalGamma]^a_bc", "non-zero" -> count[gr["Christoffel"]]|>,
   <|"tensor" -> "Riemann  R^a_bcd",                  "non-zero" -> count[gr["Riemann"]]|>,
   <|"tensor" -> "Ricci  R_ab",                       "non-zero" -> count[gr["Ricci"]]|>,
   <|"tensor" -> "Einstein  G_ab",                    "non-zero" -> count[gr["Einstein"]]|>
  }
 ];

(*::md::
## 3. Metric library

The next cell is one long comment: nothing in it is evaluated, it is just a shelf
to copy from. Take a block, paste it into **Your metric**, and evaluate.
::*)

(*::code::*)
(* ============================================================================
   METRIC LIBRARY -- copy a block into the "Your metric" cell below.
   ============================================================================

   -- Minkowski, Cartesian ---------------------------------------------------
   coords = {t, x, y, z};
   metric = DiagonalMatrix[{-1, 1, 1, 1}];
   assumptions = True;

   -- Minkowski, spherical ---------------------------------------------------
   coords = {t, r, \[Theta], \[Phi]};
   metric = DiagonalMatrix[{-1, 1, r^2, r^2 Sin[\[Theta]]^2}];
   assumptions = r > 0 && 0 < \[Theta] < Pi;

   -- Schwarzschild (mass M) -------------------------------------------------
   coords = {t, r, \[Theta], \[Phi]};
   metric = DiagonalMatrix[{-(1 - 2 M/r), 1/(1 - 2 M/r), r^2, r^2 Sin[\[Theta]]^2}];
   assumptions = M > 0 && r > 2 M && 0 < \[Theta] < Pi;

   -- Reissner-Nordstrom (mass M, charge Q) ----------------------------------
   coords = {t, r, \[Theta], \[Phi]};
   With[{f = 1 - 2 M/r + Q^2/r^2},
    metric = DiagonalMatrix[{-f, 1/f, r^2, r^2 Sin[\[Theta]]^2}]];
   assumptions = M > 0 && r > 0 && 0 < \[Theta] < Pi;

   -- Schwarzschild-de Sitter (mass M, cosmological constant \[CapitalLambda]) ---
   coords = {t, r, \[Theta], \[Phi]};
   With[{f = 1 - 2 M/r - \[CapitalLambda] r^2/3},
    metric = DiagonalMatrix[{-f, 1/f, r^2, r^2 Sin[\[Theta]]^2}]];
   assumptions = M > 0 && r > 0 && 0 < \[Theta] < Pi;

   -- Kerr, Boyer-Lindquist (mass M, spin a) ---------------------------------
   coords = {t, r, \[Theta], \[Phi]};
   With[{SS = r^2 + a^2 Cos[\[Theta]]^2, DD = r^2 - 2 M r + a^2},
    metric = {
      {-(1 - 2 M r/SS), 0,       0,  -2 M r a Sin[\[Theta]]^2/SS},
      {0,               SS/DD,   0,   0},
      {0,               0,       SS,  0},
      {-2 M r a Sin[\[Theta]]^2/SS, 0, 0,
       (r^2 + a^2 + 2 M r a^2 Sin[\[Theta]]^2/SS) Sin[\[Theta]]^2}
    }];
   assumptions = M > 0 && a > 0 && r > 0 && 0 < \[Theta] < Pi;

   -- FLRW, spatially flat, Cartesian (scale factor a[t]) --------------------
   coords = {t, x, y, z};
   metric = DiagonalMatrix[{-1, a[t]^2, a[t]^2, a[t]^2}];
   assumptions = a[t] > 0;

   -- FLRW with spatial curvature k = -1, 0, +1 ------------------------------
   coords = {t, r, \[Theta], \[Phi]};
   metric = DiagonalMatrix[{-1, a[t]^2/(1 - k r^2), a[t]^2 r^2,
                            a[t]^2 r^2 Sin[\[Theta]]^2}];
   assumptions = a[t] > 0 && r > 0 && 0 < \[Theta] < Pi && 1 - k r^2 > 0;

   -- de Sitter, static patch (cosmological constant \[CapitalLambda]) ----------
   coords = {t, r, \[Theta], \[Phi]};
   With[{f = 1 - \[CapitalLambda] r^2/3},
    metric = DiagonalMatrix[{-f, 1/f, r^2, r^2 Sin[\[Theta]]^2}]];
   assumptions = \[CapitalLambda] > 0 && r > 0 && 0 < \[Theta] < Pi;

   -- Anti-de Sitter, Poincare patch (radius L, boundary at z -> 0) ----------
   coords = {t, x, y, z};
   metric = (L/z)^2 DiagonalMatrix[{-1, 1, 1, 1}];
   assumptions = L > 0 && z > 0;

   -- Round 2-sphere of radius R0 (quick sanity check: R = 2/R0^2) -----------
   coords = {\[Theta], \[Phi]};
   metric = R0^2 DiagonalMatrix[{1, Sin[\[Theta]]^2}];
   assumptions = R0 > 0 && 0 < \[Theta] < Pi;

   ============================================================================ *)

(*::md::
## 4. Your metric

Edit this cell. The assumptions are optional but they pay for themselves: they let
Simplify resolve signs of square roots and trigonometric factors, which usually turns
a page of unreadable output into two terms.
::*)

(*::code::*)
ClearAll[t, r, x, y, z, \[Theta], \[Phi], M, a, Q, k, L, R0, \[CapitalLambda]];

(* ---- Schwarzschild, the default ---- *)
coords = {t, r, \[Theta], \[Phi]};
metric = DiagonalMatrix[{-(1 - 2 M/r), 1/(1 - 2 M/r), r^2, r^2 Sin[\[Theta]]^2}];
assumptions = M > 0 && r > 2 M && 0 < \[Theta] < Pi;

MatrixForm[metric, TableHeadings -> {coords, coords}]

(*::md::
## 5. Run
::*)

(*::code::*)
AbsoluteTiming[gr = GRTensors[metric, coords, "Assumptions" -> assumptions];]

(*::code::*)
ShowSummary[gr]

(*::md::
## 6. Christoffel symbols

$\Gamma^{\lambda}{}_{\mu\nu}$ is symmetric in its lower pair, so each independent symbol
is printed once; the mirrored component $\Gamma^{\lambda}{}_{\nu\mu}$ equals it.
::*)

(*::code::*)
ShowChristoffel[gr]

(*::md::
## 7. Riemann tensor

Printed with all indices down, $R_{\rho\sigma\mu\nu}$, where the symmetries are visible:
antisymmetry in the first pair, antisymmetry in the second, and symmetry under exchange
of the two pairs. Only one representative of each orbit is listed.

The mixed tensor $R^{\rho}{}_{\sigma\mu\nu}$, the one you want for geodesic deviation,
is available as gr["Riemann"].
::*)

(*::code::*)
ShowRiemann[gr]

(*::md::
## 8. Ricci tensor

$R_{\mu\nu} = R^{\lambda}{}_{\mu\lambda\nu}$ is symmetric, so only $\mu \le \nu$ is shown.
::*)

(*::code::*)
ShowRicci[gr]

(*::md::
## 9. Ricci scalar
::*)

(*::code::*)
ShowRicciScalar[gr]

(*::md::
## 10. Einstein tensor

$G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2} g_{\mu\nu} R$. For a vacuum solution this table is empty.
::*)

(*::code::*)
ShowEinstein[gr]

(*::md::
## 11. Checks and invariants

The Ricci scalar vanishes for every vacuum solution, so it cannot tell you whether
spacetime is flat. The Kretschmann scalar $K = R_{\rho\sigma\mu\nu} R^{\rho\sigma\mu\nu}$ can:
for Schwarzschild it gives $48 M^2 / r^6$, finite at the horizon and divergent at $r = 0$.

`gr["Resources"]` carries the underlying repository objects, so the rest of their
properties stay reachable, for instance
`gr["Resources"]["Riemann"]["BianchiIdentities"]` or `["ChernPontryaginScalar"]`.
::*)

(*::code::*)
<|
  "vacuum (G_ab = 0)"  -> (Simplify[Flatten[gr["Einstein"]], assumptions] === Table[0, gr["Dimension"]^2]),
  "flat (R^a_bcd = 0)" -> (Simplify[Flatten[gr["Riemann"]], assumptions] === Table[0, gr["Dimension"]^4]),
  "Kretschmann K"      -> gr["Kretschmann"],
  "curvature singularities" -> gr["Resources"]["Riemann"]["CurvatureSingularities"]
|> // Dataset
