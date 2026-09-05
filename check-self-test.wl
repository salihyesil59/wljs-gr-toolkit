(* ============================================================================
   check-self-test.wl -- a source that check.wls is *supposed* to fail on.

       wolframscript -file check.wls check-self-test.wl

   Expected: 6 checks, 3 pass, 3 FAILED, 1 error, exit code 1.

   A checker that never reports a failure is indistinguishable from one that
   reports nothing, and every notebook in this series currently passes -- so a
   green run proves nothing on its own. This file is the negative control: it
   contains one of each way a check can go wrong, and the runner has to see all
   of them.

   The filename deliberately does not start with GR-, so the default run
   (wolframscript -file check.wls, which globs GR-*.wl) leaves it alone. It is
   only ever run on purpose.
   ============================================================================ *)

(*::md::
# Negative control

Not physics. Each cell below is a way a check can fail.
::*)

(*::code::*)
(* Two that pass, so the count of passes is checked as well as the count of
   failures -- a runner that reported everything as failing would also "catch"
   the failures below. *)
Dataset @ {
  <|"check" -> "SELF-TEST: this one passes", "ok" -> True|>,
  <|"check" -> "SELF-TEST: so does this one", "ok" -> (Simplify[Sin[x]^2 + Cos[x]^2 - 1] === 0)|>
}

(*::code::*)
(* The ordinary failure: ok -> False. *)
Dataset @ {
  <|"check" -> "SELF-TEST: must be reported as FAIL", "ok" -> False|>
}

(*::code::*)
(* A check that never evaluated to a boolean. This is the dangerous one: an
   association whose "ok" is a leftover symbolic expression is not False, so
   anything testing `=== False` would call it a pass. *)
Dataset @ {
  <|"check" -> "SELF-TEST: ok is not a boolean at all", "ok" -> undefinedThing[3]|>
}

(*::code::*)
(* A row outside a Dataset, under the GR-05 spelling of the label. Both have to
   be found, or two notebooks' worth of checks go unseen. *)
{<|"statement" -> "SELF-TEST: bare row, statement key, fails", "ok" -> (1 === 2)|>}

(*::code::*)
(* A cell that dies. The failure has to be attributed to a cell rather than
   silently ending the file. Abort rather than a syntax error, because the whole
   point of the .wl format is that it stays directly evaluatable with Get, and a
   deliberate syntax error here would break that. *)
Abort[]

(*::code::*)
(* Evaluation must continue past the cell that died, or a single broken cell
   would hide every check after it. *)
Dataset @ {
  <|"check" -> "SELF-TEST: reached after the Abort", "ok" -> True|>
}
