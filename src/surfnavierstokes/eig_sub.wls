#!/usr/bin/env wolframscript
(* ::Package:: *)

If[Length@$ScriptCommandLine < 3, Print["Error: please specify testcase path, eps"]; Quit[]];
inDir = $ScriptCommandLine[[2]]<>"/matrices/";
eps = ToExpression[$ScriptCommandLine[[3]]];
If[Length@$ScriptCommandLine >= 4, num = ToExpression[$ScriptCommandLine[[4]]], num = 5];
Print["eps = ", eps];
Print["computing ", num, " smallest eigs"];


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
	A=.5(A+Transpose@A);
	B = import[inDir<>"B"];
	M = import[inDir<>"M"];
	M = .5(M + Transpose@M);
	Cn = import[inDir<>"C"];
	Cn = .5(Cn + Transpose@Cn);
	a = ArrayFlatten[{{A, Transpose@B}, {B, Cn}}];
	b = ArrayFlatten[{{M, 0.}, {0., SparseArray[{{i_,i_} -> eps}, {Length@B, Length@B}]}}];
	Print["size of A      = ", Length@A];
	Print["size of M      = ", Length@M];
	Print["size of B      = ", Length@B, " x ", Length@First@B];
	Print["types: ", Head/@{A, M, B, a, b}];
	Print["max memory used:  ", ToString@TeXForm@ScientificForm[MaxMemoryUsed[] * 10.^-6, 3], " MB"];
	Print["memory available: ", ToString@TeXForm@ScientificForm[MemoryAvailable[] * 10.^-6, 3], " MB"];
];
tg["import matrices (" <> inDir <> "*), n = " <> ToString@Length@A, time];


time = First@AbsoluteTiming[
	l = Eigenvalues[{a, b}, -num, Method->{"Arnoldi"}]; 
];
tg["lambda_min = " <> ToString[TeXForm@ScientificForm[#, 3]&/@l], time];
Print[Export[inDir <> "lambda_eps=" <> ToString@DecimalForm[eps, Infinity] <> "_num=" <> ToString@num <>".json", <| "n" -> Length@a, "lambda_min" -> l, "time" -> time |>]];
