#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ============================================================================
   patchwln.wls -- replace a passage inside a built .wln notebook, in place.

   Usage
       wolframscript -file patchwln.wls notebook.wln patch.txt [--count n]

   The patch file holds the old and the new text, separated by marker lines:

       %%%OLD%%%
       ...the exact text to find...
       %%%NEW%%%
       ...the text to put there...

   Why this exists
       Rebuilding a notebook with wl2wln.wls regenerates it from the .wl source
       and therefore discards every saved evaluation output. When all you want
       is to reword a paragraph, editing the .wln in place keeps the results.

   What it refuses to do
       A Markdown passage lives in a .wln exactly twice: once in the hidden
       ".md" input cell and once in its rendered "Display: markdown" output
       cell. If the count is anything other than the expected number, the file
       is left alone -- a partial replacement would leave the source and the
       rendered copy disagreeing.

       It also refuses replacement text containing a blank line followed by a
       "%", because "\n\n%" is how the reader recognises a cell separator and
       introducing one mid-cell corrupts the notebook.

   Always edit the matching .wl source too, or the next rebuild silently
   reverts the change.
   ============================================================================ *)

readText[path_String] := FromCharacterCode[Import[path, "Byte"], "UTF8"];

writeText[path_String, text_String] := Module[{stream},
  stream = OpenWrite[path, BinaryFormat -> True];
  BinaryWrite[stream, ToCharacterCode[text, "UTF8"]];
  Close[stream]];

(* verify with WLJS's own reader when the app is installed, otherwise skip *)
parseCheck[path_String] := Module[{ff, stream, res},
  ff = FileNameJoin[{$HomeDirectory, "AppData", "Local", "Programs",
     "WLJS Notebook", "resources", "app", "Kernel", "FileFormat.wl"}];
  If[! FileExistsQ[ff], Return["skipped, WLJS install not found"]];
  Quiet@Get[ff];
  stream = OpenRead[path, DOSTextFormat -> False];
  res = Quiet@CoffeeLiqueur`Notebook`FileFormat`readNotebook[stream, 10];
  Close[stream];
  If[FailureQ[res] || ! AssociationQ[res],
   "FAILED to parse",
   ToString[Length[res["Cells"]]] <> " cells, " <>
    ToString[Count[res["Cells"], c_ /; Lookup[c, "Type", ""] === "Output"]] <>
    " outputs"]];

Module[{args, nbPath, patchPath, want, raw, parts, old, new, src, n, out, before, after},

  args = Rest[$ScriptCommandLine];
  If[Length[args] < 2,
   Print["usage: wolframscript -file patchwln.wls notebook.wln patch.txt [--count n]"];
   Exit[1]];

  nbPath    = args[[1]];
  patchPath = args[[2]];
  want = With[{p = FirstPosition[args, "--count"]},
    If[MissingQ[p], 2, ToExpression[args[[First[p] + 1]]]]];

  If[! FileExistsQ[nbPath],    Print["error: no such notebook: ", nbPath];  Exit[1]];
  If[! FileExistsQ[patchPath], Print["error: no such patch file: ", patchPath]; Exit[1]];

  raw = readText[patchPath];
  parts = StringSplit[raw, {"%%%OLD%%%\n", "%%%NEW%%%\n"}];
  If[Length[parts] < 2,
   Print["error: patch file needs a %%%OLD%%% and a %%%NEW%%% marker line"]; Exit[1]];

  (* markers sit on their own lines, so drop the newline each block ends with *)
  old = StringDrop[parts[[-2]], -1];
  new = StringDrop[parts[[-1]], -1];

  If[StringContainsQ[new, "\n\n%"],
   Print["error: replacement contains a blank line followed by '%', which the "
     <> "reader would take for a cell separator. Nothing written."];
   Exit[1]];

  src = readText[nbPath];
  n = Length[StringPosition[src, old, Overlaps -> False]];
  Print["occurrences found : ", n, "   (expected ", want, ")"];

  If[n =!= want,
   Print["error: refusing to write. A markdown passage should appear twice, "
     <> "once in the hidden .md input cell and once in its rendered output cell. "
     <> "Pass --count to override deliberately."];
   Exit[1]];

  before = parseCheck[nbPath];
  out = StringReplace[src, old -> new];
  writeText[nbPath, out];
  after = parseCheck[nbPath];

  Print["bytes             : ", StringLength[src], " -> ", StringLength[out]];
  Print["before            : ", before];
  Print["after             : ", after];
  If[StringContainsQ[after, "FAILED"],
   Print["the patched file no longer parses -- restore it from git"]; Exit[1]];
  Print["ok"];
];
