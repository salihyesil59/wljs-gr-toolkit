#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ============================================================================
   wl2wln.wls -- convert an annotated Wolfram Language source file into a
   WLJS Notebook document (.wln, current text-based file format).

   Usage
       wolframscript -file wl2wln.wls  input.wl  [output.wln]

   Markup understood in the .wl source
       (*::md::                          start of a Markdown block
       ... markdown text (KaTeX ok) ...
       ::*)                              end of the Markdown block

       (*::code::*)                      start of a Wolfram input cell
       ... Wolfram code ...              (runs until the next marker)

   Everything outside a marked block is ignored, so the .wl file stays a
   valid, directly evaluatable Wolfram source file.

   A Markdown block is emitted as a hidden ".md" input cell plus its rendered
   markdown output cell, which is exactly what WLJS Notebook writes itself.
   ============================================================================ *)

$sepWidth = 75;
$objectFields = {"Objects", "Symbols", "Storage", "ExcalidrawImages", "RuntimeCache"};

(* ---------------------------------------------------------------- separators *)

sepLine[] := "%" <> StringRepeat["-", $sepWidth - 2] <> "%";

sepLine[name_String] := With[{fill = $sepWidth - StringLength[name] - 4},
  "%" <> StringRepeat["-", Floor[fill/2]] <> "%" <> name <> "%" <>
    StringRepeat["-", Ceiling[fill/2]] <> "%"
];

(* -------------------------------------------------------------- mini-YAML out *)

yamlValue[True] := "true";
yamlValue[False] := "false";
yamlValue[s_String] := s;
yamlValue[x_] := ToString[x, InputForm];

cellHeader[meta_Association] := StringRiffle[
  Flatten[{
    If[KeyExistsQ[meta, "Display"], "Display: " <> meta["Display"], Nothing],
    If[KeyExistsQ[meta, "Props"] && Length[meta["Props"]] > 0,
      {"Props:", KeyValueMap["  " <> #1 <> ": " <> yamlValue[#2] &, meta["Props"]]},
      Nothing
    ],
    "Type: " <> Lookup[meta, "Type", "Input"]
  }],
  "\n"
];

(* ------------------------------------------------------------------- parsing *)

parseSource[path_String] := Module[{lines, cells = {}, mode = None, buf = {}, flush, txt},

  (* Record + NullRecords keeps the blank lines that Markdown paragraphs depend on;
     ReadList[..., String] silently drops them. *)
  lines = ReadList[path, Record, RecordSeparators -> {"\r\n", "\n", "\r"}, NullRecords -> True];

  flush[] := (
    txt = StringTrim[StringRiffle[buf, "\n"]];
    If[txt =!= "",
      Which[
        mode === "md",
          AppendTo[cells, {<|"Type" -> "Input", "Props" -> <|"Hidden" -> True|>|>, ".md\n" <> txt}];
          AppendTo[cells, {<|"Type" -> "Output", "Display" -> "markdown"|>, txt}],
        mode === "code",
          AppendTo[cells, {<|"Type" -> "Input", "Props" -> <|"Hidden" -> False|>|>, txt}]
      ]
    ];
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

(* ------------------------------------------------------------------ assembly *)

buildNotebook[cells_List] := Module[{out = {}, n = Length[cells]},

  AppendTo[out, sepLine["Notebook"] <> "\n"];
  AppendTo[out, StringRiffle[Prepend["  - " <> # & /@ $objectFields, "ObjectFields:"], "\n"] <> "\n"];
  AppendTo[out, sepLine[] <> "\n\n"];
  AppendTo[out, "\n"];
  AppendTo[out, sepLine["Cells"] <> "\n"];

  Do[
    AppendTo[out, cellHeader[cells[[i, 1]]] <> "\n"];
    AppendTo[out, sepLine[] <> "\n"];
    AppendTo[out, cells[[i, 2]] <> "\n\n"];
    AppendTo[out,
      If[i < n, sepLine[] <> "\n", sepLine["EndOfCells"] <> "\n\n\n"]
    ],
    {i, n}
  ];

  Do[
    AppendTo[out, "\n" <> sepLine[f] <> "\n" <> sepLine["EndOf" <> f] <> "\n\n\n"],
    {f, $objectFields}
  ];

  AppendTo[out, "\n\n\n" <> sepLine["EndOfNotebook"]];

  StringJoin[out]
];

(* --------------------------------------------------------------------- main *)

Module[{args, in, out, cells, text, stream},
  args = Rest[$ScriptCommandLine];

  If[Length[args] < 1,
    Print["usage: wolframscript -file wl2wln.wls input.wl [output.wln]"];
    Exit[1];
  ];

  in = args[[1]];
  If[!FileExistsQ[in], Print["error: no such file: ", in]; Exit[1]];

  out = If[Length[args] >= 2, args[[2]], FileBaseName[in] <> ".wln"];
  out = If[DirectoryName[out] === "", FileNameJoin[{DirectoryName[AbsoluteFileName[in]], out}], out];

  cells = parseSource[in];
  If[cells === {}, Print["error: no (*::md::*) or (*::code::*) blocks found in ", in]; Exit[1]];

  text = buildNotebook[cells];

  stream = OpenWrite[out, BinaryFormat -> True];
  WriteString[stream, text];
  Close[stream];

  Print["wrote ", out, "  (", Length[cells], " cells, ", StringLength[text], " bytes)"];
];
