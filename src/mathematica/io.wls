#!/usr/bin/env wolframscript
(* ::Package:: *)

io = <|
	"in" -> "/home/alexander/Documents/drops/output",
	"out" -> "/home/alexander/Documents/Trace-FEM",
	"runs" -> {"spherical_harmonics_u_N"}
|>;
(* modules *)
Get[NotebookDirectory[] <> "TeXTableForm.m"]
ResourceFunction["MaTeXInstall"][]
<<MaTeX`
texStyle = {FontFamily -> "Latin Modern Roman", FontSize -> 12};
(* functions *)
str2exp[s_] := Quiet@ToExpression@StringReplace[ToString@s, {"e+" :> "*^", "e-" :> "*^-", "\\" :> ""}]
exportTable[table_, name_] := Module[{tex, path = io["out"] <> "/tab/" <> name <> ".tex"},
	TeXTableForm[table, 1, path];
	tex = Append[#[[;;-2]], "\\\\ \hline"]& /@ Import[path, "Table"][[2;;-2]];
	Export[path, tex, "Table"]
]
trapezoidalRule[val_, dt_] := dt/2 (First@val + 2 Total[val[[2;;-2]]] + Last@val)
spaceTimeNorm[val_, dt_] := dt Total[val]
order[table_, i_] := Module[{col = Transpose[table][[i]], stepSize = h /@ io["runs"]},
	Prepend[Table[Log[col[[j]] / col[[j - 1]]]/Log[stepSize[[j]] / stepSize[[j - 1]]], {j, 2, Length@col}], ""]
]
insertCol[table_, col_, i_] := Transpose[Join[Transpose[table][[;;i - 1]], {col}, Transpose[table][[i;;]]]]
applyCol[f_, table_, I_] := Module[{res = Transpose@table}, 
	Do[res[[i]] = f /@ res[[i]], {i, I}];
	Transpose@res
]
(* plotting *)
plotmarkers = "OpenMarkers";
markers[data_, plot_:ListPlot, step_:1] := plot[#[[;;;;step]]& /@ data, PlotMarkers -> plotmarkers, Axes -> False]
tickF = Charting`ScaledTicks["Linear"][##] /. {a_, _NumberForm, b___} :> {a, ToString @ PercentForm[a], b} &;


(* ::Section:: *)
(*Import data*)


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
(*Stats*)


kinT = Map[{#["t"], #["Integral", "FESolution", "KineticEnergy"]}&, files /@ io["runs"], {2}];
kinN = Map[{#["t"], #["SurfaceSpeedKineticEnergy"]}&, files /@ io["runs"], {2}];
Show[
	ListPlot[Join[kinN, kinT], 
		BaseStyle -> texStyle, Joined -> True, 
		Frame -> True,
			FrameLabel -> {MaTeX["t", Magnification -> 1.5], None},
		PlotLegends -> PointLegend[
			{
				MaTeX["\\frac12 \|u_N(\cdot, t)\|^2_{L^2(\Gamma_h)}", Magnification -> 1.5], 
				MaTeX["\\frac12 \|{\\bf u}_h(\cdot, t)\|^2_{L^2(\Gamma_h)}", Magnification -> 1.5]
			}, LegendMarkers -> plotmarkers]
	],
	markers[Join[kinN, kinT], ListPlot, 10]
]


run = files@First@io["runs"];
base = (First@run)["Integral", "SurfacaAreaP1"];
err = {#["t"], (base - #["Integral", "SurfacaAreaP1"])/base}& /@ run;
ListPlot[err, 
	BaseStyle -> texStyle, Joined -> True, 
	Frame -> True,
		FrameLabel -> {MaTeX["t", Magnification -> 1.5], None},
		FrameTicks-> {{tickF, None}, {Automatic, None}},
	PlotLegends -> MaTeX["\\frac{|\Gamma_h(0)| - |\Gamma_h(t)|}{|\Gamma_h(0)|}", Magnification -> 1.5]
]


ListPlot[Map[{#["t"], #["Integral", "FESolution", "KineticEnergy"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"kinetic energy", PlotLegends -> io["runs"]]
ListPlot[Map[{#["t"], #["Integral","SurfacaAreaP1"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"surface area", PlotLegends -> io["runs"], AspectRatio -> 1]
ListLogPlot[Map[{#["t"], #["NumTetras", "Velocity"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"# of active elements", PlotLegends -> io["runs"], PlotRange -> All]
ListPlot[Map[{#["t"], #["Solver", "Outer", "ResidualNorm", "r_i"]}&, files /@ io["runs"], {2}], Filling -> Bottom, PlotLabel -> "residual norm", PlotLegends -> io["runs"], PlotRange -> All]
ListPlot[Map[{#["t"], #["MaxSurfaceSpeed"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"max |u_N|", PlotLegends -> io["runs"], PlotRange -> {-0.01, Automatic}]
ListPlot[Map[{#["t"], #["NarrowBandWidth"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"narrow-band width", PlotLegends -> io["runs"], AspectRatio -> 1, AxesOrigin-> {0, 0}, PlotRange -> {-0.01, All}]
ListPlot[Map[{#["t"], #["Solver", "Outer", "TotalIters"]}&, files /@ io["runs"], {2}], Filling -> Bottom, PlotLabel->"# iters", PlotLegends -> io["runs"], PlotRange -> All]


(* ::Section::Closed:: *)
(*Error*)


ListLogPlot[Map[{#["t"], #["Integral", "Error", "VelocityH1"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"VelocityH1 err", PlotLegends -> io["runs"], PlotRange -> All]
ListLogPlot[Map[{#["t"], #["Integral", "Error", "VelocityL2"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"VelocityL2 err", PlotLegends -> io["runs"], PlotRange -> All]
ListLogPlot[Map[{#["t"], #["Integral", "Error", "PressureL2"]}&, files /@ io["runs"], {2}], Joined -> True, PlotLabel->"PressureL2 err", PlotLegends -> io["runs"], PlotRange -> All]


header = {"h", "# d.o.f.", "velH1", "order", "velL2", "order", "preL2", "order"};
table = Table[{
	h[run],
	Mean[#["DOF", "Velocity"] + #["DOF", "Pressure"]& /@ files[run]],
	spaceTimeNorm[#["Integral", "Error", "VelocityH1"]& /@ files[run], files[run][[2]]["t"]],
	spaceTimeNorm[#["Integral", "Error", "VelocityL2"]& /@ files[run], files[run][[2]]["t"]],
	spaceTimeNorm[#["Integral", "Error", "PressureL2"]& /@ files[run], files[run][[2]]["t"]]
}, {run, io["runs"]}];
orderVelH1 = order[table, 3];
orderVelL2 = order[table, 4];
orderPreL2 = order[table, 5];
table = insertCol[table, orderVelH1, 4];
table = insertCol[table, orderVelL2, 6];
table = insertCol[table, orderPreL2, 8];
table = applyCol[ScientificForm[N@#, 3]&, table, {2}];
table = applyCol[ScientificForm[#, 2]&, table, {1, 3, 7, 5}];
table = applyCol[NumberForm[#, 3]&, table, {4, 6, 8}];
table = Prepend[table, header];
Grid[table, Frame -> All]


exportTable[table, "evolving_conv"]
