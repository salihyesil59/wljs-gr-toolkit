#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ============================================================================
   check.wls -- run the acceptance checks in the GR notebook sources, headless.

   Usage
       wolframscript -file check.wls                 all GR-*.wl, in order
       wolframscript -file check.wls GR-09-....wl    just those files
       wolframscript -file check.wls --verbose       list every check, not only
                                                     the failures
       wolframscript -file check.wls --parse-only    read the sources and count
                                                     cells, evaluate nothing
       wolframscript -file check.wls --timeout=600   seconds allowed per
                                                     expression (0 = no cap)

   The exit code is 0 only when every check passed and no cell errored, so this
   can be wired to anything that reads exit codes.

   What it looks for
       Every notebook ends with verification cells that build rows of the form

           <|"check" -> "...", "ok" -> True|>

       usually inside a Dataset. Opening the notebook is currently the only way
       to see them. This script reads the .wl source with the same cell reader
       wl2wln.wls uses, evaluates the code cells in order, and collects every
       association carrying an "ok" key -- wherever it sits, Dataset or not.

       The label key is not uniform across the series: GR-02, GR-03, GR-04,
       GR-06, GR-07 and GR-09 write "check", GR-05 and GR-08 write "statement",
       and eleven rows write "acceptance test". All three are read here rather
       than edited into agreement, because changing a notebook means rebuilding
       its .wln, and that discards the saved outputs which make the results
       visible on GitHub. A row spelled some fourth way still gets reported,
       under whatever string it carries, so a new spelling shows up as a check
       with an odd label instead of vanishing.

   Why a single process
       The Wolfram licence is single-seat, so this does not spawn a kernel per
       file -- two at once would not both start. Files run one after another in
       this kernel, with Global` cleared between them: the notebooks reuse names
       like f, H and a freely, and a leftover definition would quietly change
       the next file's answers rather than erroring.
   ============================================================================ *)

Begin["GRCheck`"];

$timeout = 1200;
$verbose = False;
$parseOnly = False;

(* Distinct markers, so a cell that legitimately returns $Failed or Null is not
   mistaken for one that died. *)
$timedOut = "GRCheck: timed out";
$wasAborted = "GRCheck: aborted";

(* ------------------------------------------------------------------ parsing *)

(* The same markers and the same Record-based reader as wl2wln.wls, so the two
   agree about where a cell begins and ends. *)

codeCells[path_String] := Module[{lines, cells = {}, mode = None, buf = {}, flush, txt},
  lines = ReadList[path, Record, RecordSeparators -> {"\r\n", "\n", "\r"}, NullRecords -> True];

  flush[] := (
    txt = StringTrim[StringRiffle[buf, "\n"]];
    If[txt =!= "" && mode === "code", AppendTo[cells, txt]];
    buf = {};
  );

  Do[
    Which[
      StringTrim[line] === "(*::md::",     (flush[]; mode = "md"),
      StringTrim[line] === "(*::code::*)", (flush[]; mode = "code"),
      StringTrim[line] === "::*)",         (flush[]; mode = None),
      True,                                 If[mode =!= None, AppendTo[buf, line]]
    ],
    {line, lines}
  ];
  flush[];

  cells
];

(* A cell holds several top-level expressions, and ToExpression would return only
   the last one -- usually, but not always, where the checks are. Read them
   individually instead, so every result gets looked at and a failure can name
   the cell it came from.

   A syntax error makes ReadList stop early and emit a message, leaving a
   perfectly well-formed but truncated list of expressions behind. Silently
   evaluating that prefix is the worst outcome available, so any message during
   the read fails the whole cell. *)

heldExpressions[text_String] := Module[{stream, out, msgs},
  stream = StringToStream[text];
  (* $Messages = {} and not Quiet. Quiet suppresses the *recording* as well as
     the printing -- under it $MessageList stays empty -- so the test below
     would never fire and a truncated cell would be accepted in silence.
     Emptying the output-stream list silences the printing only. *)
  out = Block[{$MessageList = {}, $Messages = {}},
    With[{r = ReadList[stream, Hold[Expression]]}, msgs = $MessageList; r]
  ];
  Close[stream];
  If[MatchQ[out, {___Hold}] && msgs === {}, out, $Failed]
];

(* ------------------------------------------------------- reading check rows *)

$labelKeys = {"check", "statement", "acceptance test", "test", "claim"};

rowLabel[a_Association] := Module[{key, fallback},
  key = SelectFirst[$labelKeys, KeyExistsQ[a, #] &];
  If[StringQ[key],
    ToString[a[key]],
    fallback = SelectFirst[Normal[KeyDrop[a, "ok"]], StringQ[Last[#]] &];
    If[MissingQ[fallback],
      StringTake[ToString[a, InputForm], UpTo[70]],
      ToString[Last[fallback]]
    ]
  ]
];

(* Datasets are unwrapped first, so that a Dataset holding its rows in some
   internal form still yields plain associations.

   The rows are found by KeyExistsQ and not by searching for the string "ok":
   an association's keys are not subexpressions of it, so FreeQ[a, "ok"] is
   True even when a is <|"ok" -> True|>. A guard written that way silently
   reports every notebook as having no checks. *)
checkRows[expr_] := Module[{plain},
  plain = expr //. d_Dataset :> Normal[d];
  If[AtomQ[plain],
    {},
    Cases[plain, a_Association /; KeyExistsQ[a, "ok"], {0, Infinity}]
  ]
];

(* ---------------------------------------------------------------- one file *)

(* Messages are collected and summarised per file, so printing them as well
   would interleave a paragraph of prose with the report -- GR-04 demonstrates
   its own GrowthSolution::pole guard every run, and that message is four lines
   long.

   Emptying $Messages, the list of streams messages are written to, silences the
   printing. Quiet would silence the recording too: under Quiet, $MessageList
   stays empty and the per-file message summary reports nothing, which is how
   this was found. *)

evaluateHeld[h_] := Block[
  {$Context = "Global`", $ContextPath = {"Global`", "System`"}, $Messages = {}},
  CheckAbort[
    If[$timeout > 0,
      TimeConstrained[ReleaseHold[h], $timeout, $timedOut],
      ReleaseHold[h]
    ],
    $wasAborted
  ]
];

runFile[path_String] := Module[
  {name, cells, rows = {}, errors = {}, messages = {}, held, value, seconds, cellNo = 0},

  name = FileNameTake[path];
  cells = codeCells[path];

  If[cells === {},
    Print[name];
    Print["  no (*::code::*) cells found"];
    Print[""];
    Return[<|"file" -> name, "cells" -> 0, "rows" -> {},
             "errors" -> {{0, "no code cells"}}, "messages" -> {}, "seconds" -> 0.|>]
  ];

  If[$parseOnly,
    Print[StringPadRight[name, 44], Length[cells], " cells"];
    Return[<|"file" -> name, "cells" -> Length[cells], "rows" -> {},
             "errors" -> {}, "messages" -> {}, "seconds" -> 0.|>]
  ];

  seconds = First @ AbsoluteTiming[
    Do[
      cellNo++;

      held = Block[{$Context = "Global`", $ContextPath = {"Global`", "System`"}},
        heldExpressions[cell]];

      If[held === $Failed,
        AppendTo[errors, {cellNo, "could not be read"}];
        Continue[]
      ];

      Do[
        value = Block[{$MessageList = {}},
          With[{v = evaluateHeld[h]},
            messages = Join[messages, ToString /@ $MessageList];
            v
          ]
        ];

        Which[
          value === $timedOut,
            AppendTo[errors, {cellNo, "timed out after " <> ToString[$timeout] <> " s"}],
          value === $wasAborted,
            AppendTo[errors, {cellNo, "aborted"}],
          True,
            rows = Join[rows, {cellNo, #} & /@ checkRows[value]]
        ],
        {h, held}
      ],
      {cell, cells}
    ]
  ];

  <|"file" -> name, "cells" -> Length[cells], "rows" -> rows, "errors" -> errors,
    "messages" -> messages, "seconds" -> seconds|>
];

(* ------------------------------------------------------------------ report *)

passedQ[row_] := TrueQ[row[[2]]["ok"]];

plural[n_Integer, word_String] := ToString[n] <> " " <> word <> If[n === 1, "", "s"];

reportFile[r_Association] := Module[{rows, bad, n},
  If[$parseOnly, Return[]];

  rows = r["rows"];
  bad = Select[rows, ! passedQ[#] &];
  n = Length[rows];

  Print[r["file"]];
  Print["  ", r["cells"], " cells, ", ToString[NumberForm[r["seconds"], {6, 1}]], " s, ",
    Which[
      n === 0,     "no checks found",
      bad === {},  plural[n, "check"] <> ", all pass",
      True,        plural[n, "check"] <> ", " <> ToString[Length[bad]] <> " FAILED"
    ]
  ];

  Do[
    Print["    FAIL  cell ", row[[1]], ": ", rowLabel[row[[2]]],
      (* False is the ordinary failure; anything else means the check did not
         even evaluate to a boolean, which is a different kind of broken. *)
      If[row[[2]]["ok"] === False, "",
        "   (ok -> " <> StringTake[ToString[row[[2]]["ok"], InputForm], UpTo[40]] <> ")"]],
    {row, bad}
  ];

  If[$verbose,
    Do[Print["    pass  ", rowLabel[row[[2]]]], {row, Select[rows, passedQ]}]
  ];

  Do[Print["    ERROR cell ", e[[1]], ": ", e[[2]]], {e, r["errors"]}];

  If[r["messages"] =!= {},
    Print["    messages: ", StringRiffle[
      KeyValueMap[If[#2 > 1, #1 <> " x" <> ToString[#2], #1] &, Counts[r["messages"]]], ", "]]
  ];

  Print[""];
];

(* -------------------------------------------------------------------- main *)

Module[{args, files, here, results, allRows, bad, errs, total},
  args = Rest[$ScriptCommandLine];

  $verbose   = MemberQ[args, "--verbose" | "-v"];
  $parseOnly = MemberQ[args, "--parse-only"];

  Do[
    If[StringMatchQ[a, "--timeout=" ~~ __], $timeout = ToExpression[StringDrop[a, 10]]],
    {a, args}
  ];
  If[! (IntegerQ[$timeout] && $timeout >= 0),
    Print["error: --timeout wants a whole number of seconds, 0 or more"];
    Exit[2]
  ];

  here = If[StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[AbsoluteFileName[$InputFileName]],
    Directory[]
  ];

  files = Select[args, ! StringStartsQ[#, "-"] &];
  If[files === {}, files = SortBy[FileNames["GR-*.wl", here], FileNameTake]];

  If[files === {},
    Print["error: no GR-*.wl sources next to check.wls"];
    Exit[2]
  ];

  files = AbsoluteFileName /@ files;
  If[MemberQ[files, $Failed], Print["error: no such file"]; Exit[2]];

  results = Table[
    (* Clearing Global` between files is what keeps file n+1 from inheriting
       file n's definition of f, H or a. *)
    Quiet[Unprotect["Global`*"]; ClearAll["Global`*"]; Remove["Global`*"]];
    (* Reported as each file finishes rather than all at the end: the full
       series takes ten minutes, most of it inside GR-09, and a run that prints
       nothing until it is over is indistinguishable from a run that has hung. *)
    With[{r = runFile[f]}, reportFile[r]; r],
    {f, files}
  ];

  If[$parseOnly,
    Print[""];
    Print[Length[files], " files, ", Total[#["cells"] & /@ results], " code cells"];
    Exit[0]
  ];

  allRows = Join @@ (#["rows"] & /@ results);
  bad = Select[allRows, ! passedQ[#] &];
  errs = Join @@ (#["errors"] & /@ results);
  total = Total[#["seconds"] & /@ results];

  Print[StringRepeat["-", 70]];
  Print[plural[Length[files], "file"], ", ", plural[Length[allRows], "check"], ", ",
    Length[allRows] - Length[bad], " pass, ", Length[bad], " fail, ",
    plural[Length[errs], "error"], ", ", ToString[NumberForm[total/60., {4, 1}]], " min"];

  Exit[If[bad === {} && errs === {}, 0, 1]];
];

End[];
