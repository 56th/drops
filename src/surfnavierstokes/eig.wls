#!/usr/bin/env wolframscript
(* ::Package:: *)

If[Length@$ScriptCommandLine < 2, Print["Error: please specify testcase path"]; Quit[]];
inDir = $ScriptCommandLine[[2]]<>"/matrices/";


Get["https://raw.githubusercontent.com/dih5/TgBot/master/TgBot/TgBot.m"];
Needs["TgBot`"]
tgChatID = Import["tg.txt", "List"][[1]];
tgToken = Import["tg.txt", "List"][[2]];
BotAPICall["getUpdates",{},{"Token"->tgToken}];
tg[text_, time_]:=BotAPICall["sendMessage", {"chat_id"->tgChatID, "text" -> 
	"eig.wls: " <> text <> " ("<> ToString@TeXForm@ScientificForm[time / 60., 3] <>" min, " <> ToString@TeXForm@ScientificForm[MaxMemoryUsed[] * 10.^-6, 3] <> " Mb)"},
{"Token"->tgToken}];


import[path_] := N@Import[path <> ".mat", "LabeledData"][[1, 2]]


time = First@AbsoluteTiming[
	A = import[inDir<>"A"];
	M = import[inDir<>"M"];
	M = .5(M + Transpose@M);
	Print["size of A      = ", Length@A];
	Print["size of M      = ", Length@M];
	Print["types: ", Head/@{A, M}];
	Print["max memory used:  ", ToString@TeXForm@ScientificForm[MaxMemoryUsed[] * 10.^-6, 3], " MB"];
	Print["memory available: ", ToString@TeXForm@ScientificForm[MemoryAvailable[] * 10.^-6, 3], " MB"];
];
tg["import matrices (" <> inDir <> "*), n = " <> ToString@Length@A, time];


time = First@AbsoluteTiming[
	l = Eigenvalues[{A, M}, -5, Method->{"Arnoldi"}]; 
];
tg["lambda_min = " <> ToString[TeXForm@ScientificForm[#, 3]&/@l], time];
Print[Export[inDir <> "lambda.json", <| "n" -> Length@A, "lambda_min" -> l, "time" -> time |>]];
