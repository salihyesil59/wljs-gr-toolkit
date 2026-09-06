#!/usr/bin/env wolframscript
(* ============================================================================
   cosmofit-export.wls -- turn an action derived here into a runnable CosmoFit
   model, with the Wolfram answer carried along as an assertion.

       wolframscript -file cosmofit-export.wls
       wolframscript -file cosmofit-export.wls "Ts + A0 (-Ts)^bb" A0 bb=0.2

   The first form writes the built-in examples. The second takes an action in
   GR-02's own symbols, the name of the parameter to close on, and values for
   the rest.

   Why this exists
       `cosmofit-reference.wls` and `cosmofit-mu-reference.wls` compare the two
       codebases on models someone has already written twice. This goes the
       other way: derive something new here, and get the Python that fits it.

   What it will not do
       Guess. Three things make it refuse rather than emit a model:

         - the action is fourth order, so GR-02 hands back an ODE rather than a
           closed form and there is no reference value to assert;
         - the closure has no solution, or more than one;
         - E(z) leaves the reals somewhere in the sample range.

       The last one is not hypothetical. f = Q + alpha Q^2 looks perfectly
       reasonable and, once the closure fixes alpha, its discriminant goes
       negative at z = 0.06: past that the model has no real expansion history
       at all. Root-finding through it returns a confident number near 0.845
       that barely moves with redshift. An exporter that did not check would
       hand over Python for a theory that does not exist.

   The conventions, which are the whole difficulty
       GR-02 uses Ts = Qs = -6H^2, so that a linear f is GR in both sectors.
       CosmoFit's compiler uses T = -6E^2 but Q = +6E^2, and *also* calls a
       linear f the GR limit in both -- checked, not assumed: f = T - 2Lam and
       f = Q - 2Lam each reproduce LCDM to 0.0e0 there, while f = -T - 2Lam and
       f = -Q - 2Lam give NaN.

       Torsion therefore needs no translation at all. Non-metricity does, and it
       is not a rename: the two STEGR Lagrangians differ by an overall sign,
       which the E(0) = 1 closure hides for a linear f and hides for nothing
       else. The map that survives is

           f_CosmoFit(Q) = - f_GR02(-Q)

       sending STEGR to +Q and sending Qs + A0 (-Qs)^b to Q - A0 Q^b. That was
       measured rather than reasoned out, and two plausible alternatives were
       measured first: the literal Qs -> Q gives a complex closure, and Qs -> -Q
       on its own is rejected as having no real solution. Both failures were
       loud only because the emitted script asserts numbers -- which is the
       argument for emitting them.
   ============================================================================ *)

Off[General::stop];

$here = If[StringQ[$InputFileName] && $InputFileName =!= "",
   DirectoryName[AbsoluteFileName[$InputFileName]], Directory[]];

Get[FileNameJoin[{$here, "GR-02-Field-Equations-From-Action.wl"}]];

$sampleZ = {0, 1/2, 1, 2};
$Om = 3/10;
$H0 = 70;

(* Wolfram writes a 20-digit real as 1.234`20., and the precision marker is not
   Python. Strip from the backtick to the end of the literal. *)
pyFloat[x_] := StringReplace[ToString[N[x, 20], InputForm], "`" ~~ ___ -> ""];

(* ------------------------------------------------------------ translation *)

$scalarOf = {Rs -> "R", Ts -> "T", QCF -> "Q"};
$geometryOf = {Rs -> None, Ts -> "teleparallel", Qs -> "symmetric"};

(* The symmetric sector needs more than a rename, and getting it wrong is the
   whole reason the emitted script asserts numbers.

   GR-02 writes actions in Qs = -6H^2; CosmoFit's compiler uses Q = +6E^2. Both
   call a *linear* f the GR limit, so the two STEGR Lagrangians differ by an
   overall sign, which the E(0) = 1 closure hides for linear f and does not hide
   for anything else. The map that survives is

       f_CosmoFit(Q) = - f_GR02(-Q)

   which sends STEGR to +Q, as CosmoFit wants, and sends GR-02's
   Qs + A0 (-Qs)^b to Q - A0 Q^b. Measured, not reasoned: that expression
   reproduces this repository's E(z) to 0.0e0, while the literal Qs -> Q gives a
   complex closure and Qs -> -Q alone is rejected as having no real solution.

   Torsion needs none of this. GR-02 and CosmoFit both put Ts = T = -6H^2, so
   the map is the identity and the f(T) power law already agreed to 0.0e0. *)

forExport[action_, scalar_] := If[scalar === Qs,
   Expand[-(action /. Qs -> -QCF)],
   action];

geometryFor[expr_] := Module[{present},
  present = Select[{Rs, Ts, Qs}, ! FreeQ[expr, #] &];
  Which[
    Length[present] == 0, $Failed,
    Length[present] > 1, $Failed,
    True, {First[present], First[present] /. $geometryOf}]];

(* A recursive converter rather than string surgery on InputForm: "^" -> "**"
   by replacement also rewrites the ^ inside any symbol name a user picks, and
   bracket matching for Exp[...] is exactly the kind of thing that works on
   every example until it does not. *)

(* Dispatch on x_Plus and x_Times, never on Plus[a__].

   Plus[a__] does not survive being written down: Plus of a single argument
   returns that argument, so the pattern evaluates to a__ and the definition is
   stored as py[a__] -- a sequence catch-all sitting *behind* py[x_], which then
   swallows every sum and product and reports them as untranslatable. The
   symptom is that atoms convert fine and nothing else does. *)

ClearAll[py];
py[x_Integer] := ToString[x];
py[x_Rational] := "(" <> ToString[Numerator[x]] <> "/" <> ToString[Denominator[x]] <> ")";
py[x_Real] := ToString[CForm[x]];
py[Pi] := "pi";
py[s_Symbol] := With[{n = s /. $scalarOf}, If[StringQ[n], n, SymbolName[s]]];
py[x_Plus] := "(" <> StringRiffle[py /@ (List @@ x), " + "] <> ")";
py[x_Times] := "(" <> StringRiffle[py /@ (List @@ x), "*"] <> ")";
py[Power[E, b_]] := "exp(" <> py[b] <> ")";
py[Power[a_, Rational[1, 2]]] := "sqrt(" <> py[a] <> ")";
py[Power[a_, -1]] := "(1/" <> py[a] <> ")";
py[Power[a_, b_]] := "(" <> py[a] <> "**" <> py[b] <> ")";
py[Log[a_]] := "log(" <> py[a] <> ")";
py[x_] := (Print["error: cannot translate ", InputForm[x]]; Abort[]);

(* --------------------------------------------------------- the reference *)


referenceFor[action_, closureName_Symbol, values_List] := Module[
  {geom, hf, branch, eqn, closureSolutions, closureValue, samples, bad, fixed},

  geom = geometryFor[action];
  If[geom === $Failed,
    Return[<|"ok" -> False, "why" -> "the action must use exactly one of Rs, Ts, Qs"|>]];

  hf = HubbleFunction[FriedmannEquations[action],
     "Components" -> {{\[CapitalOmega]m, 0}}];

  (* GR-02's own flag, and the right one to ask. False means the constraint
     carries Ez' and Ez'' -- a fourth-order f(R). There is still an expansion
     history, solved numerically, but no closed form to put in an assertion. *)
  If[hf["Algebraic"] =!= True,
    Return[<|"ok" -> False,
      "why" -> "the constraint is a differential equation rather than an "
               <> "algebraic one, which means a fourth-order f(R). GR-02 will "
               <> "integrate it, but there is no closed form to assert against"|>]];

  (* GR-02 always supplies a "Branch" key, but when E(z) is only implicit the
     value is a failed replacement that still mentions Ez. Testing KeyExistsQ
     therefore says "closed form" for every model there is; test the contents. *)
  branch = hf["Branch"] /. values;

  fixed = {\[CapitalOmega]m -> $Om, H0 -> $H0};

  If[FreeQ[branch, Ez],

    (* closed form: solve the closure on the branch itself *)
    eqn = branch /. fixed;
    closureSolutions = Solve[(eqn /. z -> 0) == 1, closureName];
    If[Length[closureSolutions] =!= 1,
      Return[<|"ok" -> False,
        "why" -> "the closure E(0) = 1 has " <> ToString[Length[closureSolutions]]
                 <> " solutions for " <> SymbolName[closureName]|>]];
    closureValue = closureName /. First[closureSolutions];
    samples = Table[N[eqn /. First[closureSolutions] /. z -> zz, 25], {zz, $sampleZ}],

    (* implicit: root-find on the constraint, the way cosmofit-reference.wls does *)
    eqn = (hf["Equation"] /. Ez[z] -> ee) /. values /. fixed;
    closureSolutions = Solve[(eqn /. {ee -> 1, z -> 0}) == 0, closureName];
    If[Length[closureSolutions] =!= 1,
      Return[<|"ok" -> False,
        "why" -> "the closure E(0) = 1 has " <> ToString[Length[closureSolutions]]
                 <> " solutions for " <> SymbolName[closureName]|>]];
    closureValue = closureName /. First[closureSolutions];
    samples = Table[
       Quiet @ Check[
         ee /. FindRoot[(eqn /. First[closureSolutions] /. z -> zz) == 0,
            {ee, 1, 1/50, 50}, WorkingPrecision -> 30], Indeterminate],
       {zz, $sampleZ}]
  ];

  bad = Select[Transpose[{$sampleZ, samples}],
     (! FreeQ[Last[#], Complex | Indeterminate]) || ! Positive[Re[Last[#]]] &];

  If[bad =!= {},
    Return[<|"ok" -> False,
      "why" -> "E(z) is not real and positive at z = "
               <> StringRiffle[ToString[N[First[#]]] & /@ bad, ", "]
               <> "; the closure puts this model outside its own domain"|>]];

  <|"ok" -> True, "geometry" -> geom[[2]], "scalar" -> geom[[1]],
    "closure" -> N[closureValue, 20], "E" -> samples|>];

(* ------------------------------------------------------------- emission *)

emit[action_, closureName_Symbol, values_List, name_String] := Module[
  {ref, geomArg, paramLines, code},

  ref = referenceFor[action, closureName, values];

  If[! ref["ok"],
    Print["REFUSED  ", name, ": ", ref["why"]];
    Print[""];
    Return[$Failed]];

  geomArg = If[ref["geometry"] === None, "",
     "    geometry=\"" <> ref["geometry"] <> "\",\n"];

  (* the closure parameter is declared too. CosmoFit wants it present even
     though it derives the value itself, and carrying the Wolfram value as the
     default puts the two closures side by side the moment they disagree. *)
  paramLines = StringJoin[
    ("        \"" <> SymbolName[First[#]] <> "\": {\"default\": "
       <> ToString[CForm[N[Last[#]]]] <> "},\n") & /@
      Append[values, closureName -> ref["closure"]]];

  code = StringJoin[
    "# Generated by cosmofit-export.wls from the GR-02 action\n",
    "#     ", ToString[action, InputForm], "\n",
    "# The assertion below is this repository's own answer, at 20 digits.\n\n",
    "import numpy as np\n",
    "from CosmoFit.theory import Action\n\n",
    name, " = Action(\n",
    "    \"", py[forExport[action, ref["scalar"]]], "\",\n",
    geomArg,
    "    params={\n", paramLines, "    },\n",
    "    closure=\"", SymbolName[closureName], "\",\n",
    ").build(\"", name, "\")\n\n",
    "m = ", name, "(", name, ".PARAMS_CLASS(H0=", ToString[N[$H0]],
       ", Omega_m=", ToString[N[$Om]], "))\n",
    "z = np.array([", StringRiffle[ToString[N[#]] & /@ $sampleZ, ", "], "])\n",
    "wolfram = np.array([\n",
    StringJoin[("    " <> pyFloat[#] <> ",\n") & /@ ref["E"]],
    "])\n\n",
    "dev = np.max(np.abs(m.E(z) / wolfram - 1.0))\n",
    "print(f\"", name, ": max relative deviation from Wolfram = {dev:.3e}\")\n",
    "assert dev < 1e-12, \"CosmoFit and GR-02 disagree about this model\"\n"];

  Print["EXPORTED ", name, "  (", ref["geometry"] /. None -> "metric",
        ", closure ", SymbolName[closureName], " = ",
        ToString[N[ref["closure"], 8]], ")"];
  Print[code];
  code];

(* ---------------------------------------------------------------- main *)

Module[{args, action, closureName, values},
  args = Rest[$ScriptCommandLine];

  If[Length[args] >= 2,

    action = ToExpression[args[[1]]];
    closureName = Symbol[args[[2]]];
    values = (Module[{parts = StringSplit[#, "="]},
        Symbol[parts[[1]]] -> ToExpression[parts[[2]]]] &) /@ Drop[args, 2];
    emit[action, closureName, values, "Exported"],

    (* the built-in set: one that works, one that is refused for each reason *)
    emit[Ts + A0 (-Ts)^bb, A0, {bb -> 1/5}, "PowerLawFT"];
    emit[Qs + A0 (-Qs)^bb, A0, {bb -> 1/5}, "PowerLawFQ"];
    emit[Qs + A0 Qs^2, A0, {}, "QuadraticFQ"];
    emit[Rs - 2 lam + A0 Rs^2, A0, {lam -> 21/10}, "Starobinsky"];
  ]];
