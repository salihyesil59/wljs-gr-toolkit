#!/usr/bin/env wolframscript
(* ============================================================================
   Reference E(z) values for CosmoFit, straight out of GR-02's engine.

       wolframscript -file cosmofit-reference.wls

   The companion library (https://github.com/salihyesil59/CosmoFit) derives the
   same Friedmann constraints independently, in sympy, and solves them with its
   own root finder. Nothing is shared between the two but the physics, so
   agreement is evidence and disagreement is a bug in one of them.

   This prints E(z) at 25 digits for LCDM and for the f(T) power law, by
   loading THIS repository's GR-02 notebook and driving its
   `FriedmannEquations` / `HubbleFunction` directly -- not by re-deriving
   anything here, which would defeat the purpose.

   The numbers are pinned in `tests/test_notebook_agreement.py` on the CosmoFit
   side. Regenerate them here when a derivation in GR-02 changes; that test is
   meant to fail if it changes silently.

   Incidentally, GR-02 gives the f(T) and f(Q) power laws the *identical*
   background constraint, because on flat FLRW both scalars are -6H^2. That is
   correct rather than a bug, and it means background data alone cannot
   separate the two families.
   ============================================================================ *)

Off[General::stop];
Get[FileNameJoin[{DirectoryName[$InputFileName],
     "GR-02-Field-Equations-From-Action.wl"}]];

om = 3/10;
h0 = 70;

(* ---------- LCDM from the action ---------- *)
lcdmBranch = HubbleFunction[FriedmannEquations[Rs - 2 lam],
    "Components" -> {{\[CapitalOmega]m, 0}}]["Branch"];
lamSol = First[Solve[(lcdmBranch /. z -> 0) == 1, lam]];
lcdmE  = Simplify[lcdmBranch /. lamSol];

Print["LCDM-CLOSED ", InputForm[Simplify[lcdmE]]];
Do[Print["WOLFRAM LCDM z=", ToString[N[zz, 3]], " ",
    ToString[N[lcdmE /. {\[CapitalOmega]m -> om, H0 -> h0, z -> zz}, 25], InputForm]],
  {zz, {0, 1/2, 1, 2}}];

(* ---------- f(T) power law ---------- *)
eqTraw = HubbleFunction[FriedmannEquations[Ts + alp (-Ts)^bb],
    "Components" -> {{\[CapitalOmega]m, 0}}]["Equation"];

(* turn the unknown Ez[z] into a plain symbol before any z substitution *)
eqT = eqTraw /. Ez[z] -> ee;

Print["FT-EQUATION ", InputForm[eqT]];

Do[
  Module[{eqn, alpSol, f, root},
    eqn = eqT /. {\[CapitalOmega]m -> om, H0 -> h0, bb -> nn};

    (* closure: E = 1 at z = 0 fixes alp *)
    alpSol = First[Solve[(eqn /. {ee -> 1, z -> 0}) == 0, alp]];

    Do[
      f = eqn /. alpSol /. z -> zz;
      root = ee /. FindRoot[f == 0, {ee, 1, 1/20, 50}, WorkingPrecision -> 30];
      Print["WOLFRAM FT n=", ToString[N[nn, 3]], " z=", ToString[N[zz, 3]],
            " ", ToString[N[root, 25], InputForm]],
      {zz, {0, 1/2, 1, 2}}]],
  {nn, {-1/2, 1/5, 7/10}}];
