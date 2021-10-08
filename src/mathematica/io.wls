#!/usr/bin/env wolframscript
(* ::Package:: *)

io = <|
	"in" -> "/home/alexander/Documents/drops/output",
	"out" -> "/home/alexander/Documents/Trace-FEM",
	"runs" -> {"evolving_conv_lvl3", "evolving_conv_lvl4", "evolving_conv_lvl5", "evolving_conv_lvl6"}
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


ListLogPlot[Map[{#["t"], #["NumTetras", "Velocity"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"# of active elements", PlotLegends -> io["runs"], PlotRange -> All]


(* ::Section:: *)
(*Residual*)


ListPlot[Map[{#["t"], #["Solver", "Outer", "ResidualNorm", "r_i"]}&, files /@ io["runs"], {2}], Filling -> Bottom, PlotLabel -> "residual norm", PlotLegends -> io["runs"], PlotRange -> All]


(* ::Section:: *)
(*Num iters*)


ListPlot[Map[{#["t"], #["Solver", "Outer", "TotalIters"]}&, files /@ io["runs"], {2}], Filling -> Bottom, PlotLabel->"# iters", PlotLegends -> io["runs"], PlotRange -> All]


(* ::Section:: *)
(*Error*)


ListLogPlot[Map[{#["t"], #["Integral", "Error", "VelocityH1"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"VelocityH1 err", PlotLegends -> io["runs"], PlotRange -> All]
ListLogPlot[Map[{#["t"], #["Integral", "Error", "VelocityL2"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"VelocityL2 err", PlotLegends -> io["runs"], PlotRange -> All]
ListLogPlot[Map[{#["t"], #["Integral", "Error", "PressureL2"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"PressureL2 err", PlotLegends -> io["runs"], PlotRange -> All]


trapezoidalRule[val_, dt_] := dt/2 (First@val + 2 Total[val[[2;;-2]]] + Last@val)
order[table_, i_] := Module[{col = Transpose[table][[i]], stepSize = h /@ io["runs"]},
	Prepend[Table[Log[col[[j]] / col[[j - 1]]]/Log[stepSize[[j]] / stepSize[[j - 1]]], {j, 2, Length@col}], ""]
]
insertCol[table_, col_, i_] := Transpose[Join[Transpose[table][[;;i - 1]], {col}, Transpose[table][[i;;]]]]
applyCol[f_, table_, I_] := Module[{res = Transpose@table}, 
	Do[res[[i]] = f /@ res[[i]], {i, I}];
	Transpose@res
]


header = {"h", "# d.o.f.", "velH1", "order", "velL2", "order", "preL2", "order"};
table = Table[{
	h[run],
	files[run][[1]]["DOF", "Velocity"] + files[run][[1]]["DOF", "Pressure"],
	trapezoidalRule[#["Integral", "Error", "VelocityH1"]& /@ files[run], files[run][[2]]["t"]],
	trapezoidalRule[#["Integral", "Error", "VelocityL2"]& /@ files[run], files[run][[2]]["t"]],
	trapezoidalRule[#["Integral", "Error", "PressureL2"]& /@ files[run], files[run][[2]]["t"]]
}, {run, io["runs"]}];
orderVelH1 = order[table, 3];
orderVelL2 = order[table, 4];
orderPreL2 = order[table, 5];
table = insertCol[table, orderVelH1, 4];
table = insertCol[table, orderVelL2, 6];
table = insertCol[table, orderPreL2, 8];
table = applyCol[ScientificForm[#, 3]&, table, {1, 3, 7, 5}];
table = applyCol[NumberForm[#, 2]&, table, {4, 6, 8}];
Grid[Prepend[table, header], Frame -> All]
