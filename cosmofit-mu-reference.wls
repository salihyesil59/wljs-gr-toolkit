#!/usr/bin/env wolframscript
(* ============================================================================
   Reference G_eff/G values for CosmoFit, from GR-06's own derivation.

       wolframscript -file cosmofit-mu-reference.wls

   `cosmofit-reference.wls` compares the *background*. This compares the
   perturbations, which is where a modified-gravity model actually differs --
   and where a sign can hide, because the two codebases disagree about the sign
   of the torsion scalar deliberately.

   GR-06 puts it at T = +6H^2, TEGR at f = -T, and derives G_eff/G = -1/f_T.
   CosmoFit's teleparallel sector puts it at T = -6H^2 and uses +1/f_T. Two
   conventions, two minus signs, and they cancel -- so the same physical model
   must give the same number, and would not if either were flipped alone.

   The same theory in each convention:

       CosmoFit   f(T_CF) = T_CF + A0 (-T_CF)^n,   T_CF = -6E^2
       GR-06      f(T)    = -T   + A0 T^n,         T    = +6E^2

   and the second is TEGR at A0 = 0, which checks the translation before any
   number is compared.

   Pinned on the other side in `tests/test_notebook_agreement.py`. Note n = 0.7
   comes out negative and then through a pole: that region is genuinely
   unphysical, and CosmoFit refuses it rather than returning those numbers.
   ============================================================================ *)

Off[General::stop];
Get[FileNameJoin[{DirectoryName[$InputFileName],
     "GR-06-Teleparallel-Geometry.wl"}]];

Print["G6-GEFF ", InputForm[geffSubHorizon]];

(* CosmoFit's model is  f_CF(T_CF) = T_CF + A0 (-T_CF)^n  with
   T_CF = -6 E^2.  GR-06's T is +6 E^2 = -T_CF, so the same theory
   reads  f = -T + A0 T^n  here, which is TEGR at A0 = 0. *)
fOfT[T_, A0_, n_] := -T + A0 T^n;

om = 3/10;

(* A0 from the closure E(0) = 1 on the Friedmann relation
   E^2 + (Om - 1) E^(2n) = Om (1+z)^3, in GR-06's variables. *)
(* geffSubHorizon is a function of the symbol fT; substitute ours *)
muFromG6[nn_, Ez_] := Module[{A0, T, fTval},
  A0 = (1 - om)/(6^(nn - 1) (2 nn - 1));
  T = 6 Ez^2;
  fTval = -1 + A0 nn T^(nn - 1);
  N[geffSubHorizon /. fT -> fTval, 25]];

(* E(z) for the f(T) power law, from the same constraint GR-02 gives *)
Ez[nn_, z_] := Module[{ee, eqn, sol},
  eqn = ee^2 + (om - 1) ee^(2 nn) - om (1 + z)^3;
  ee /. FindRoot[eqn == 0, {ee, 1, 1/20, 50}, WorkingPrecision -> 30]];

Do[
  Module[{e},
    e = Ez[nn, zz];
    Print["G6-MU n=", ToString[N[nn, 3]], " z=", ToString[N[zz, 3]],
          " E=", ToString[N[e, 20], InputForm],
          " mu=", ToString[muFromG6[nn, e], InputForm]]],
  {nn, {-1/2, 1/5, 7/10}}, {zz, {0, 1/2, 1, 2}}];
