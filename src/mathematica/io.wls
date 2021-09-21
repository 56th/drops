#!/usr/bin/env wolframscript
(* ::Package:: *)

io = <|
	"in" -> "/home/alexander/Documents/drops/output",
	"out" -> "/home/alexander/Documents/Trace-FEM",
	"runs" -> {"evolving_conv_lvl3_gamma=0", "evolving_conv_lvl3_gamma=1", "evolving_conv_lvl4_gamma=0"}
|>;


str2exp[s_] := Quiet@ToExpression@StringReplace[ToString@s, {"e+" :> "*^", "e-" :> "*^-", "\\" :> ""}];
Do[
	SetDirectory[io["in"] <> "/" <> run];
	fileNames = Sort[
		FileNames["t_"~~__~~".json", "stats"], 
		ToExpression@StringDelete[Last@FileNameSplit[#1],{"t_",".json"}] < ToExpression@StringDelete[Last@FileNameSplit[#2], {"t_",".json"}]&
	];
	files[run]= Map[str2exp, #, {-1}]& /@ Import[#, "RawJSON"]& /@ fileNames;
	numSteps[run] = Length@fileNames;
	h[run]=files[run][[1]]["h"];
	(* hSF[testDir] = ScientificForm[h[run], 3]; *)
	input[run] = Map[str2exp, Import["input.json", "RawJSON"], {-1}],
{run, io["runs"]}]


(* ::Section:: *)
(*Num tetras (velocity)*)


ListPlot[Map[#["NumTetras", "Velocity"]&, files /@ io["runs"], {2}], PlotLegends -> io["runs"], PlotRange -> All]


(* ::Section:: *)
(*Residual*)


ListPlot[Map[#["Solver", "Outer", "ResidualNorm", "r_i"]&, files /@ io["runs"], {2}], PlotLegends -> io["runs"]]


(* ::Section:: *)
(*Num iters*)


(*
Do[
	Print@ListPlot[#["Solver", "Outer", "TotalIters"]& /@ files[run], PlotRange\[Rule]All],
{run, io["runs"]}]
*)
ListPlot[Map[#["Solver", "Outer", "TotalIters"]&, files /@ io["runs"], {2}], PlotLegends -> io["runs"]]
