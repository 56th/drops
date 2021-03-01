#!/usr/bin/env wolframscript
(* ::Package:: *)

surface = "Torus";
solution = "KillingExact";


surfaces = <|
	"Sphere" -> <|
		"dist" -> ({x, y, z, t} |-> Sqrt[x^2 + y^2 + z^2] - r[t]),
		"phi" -> ({x, y, z, t} |-> x^2 + y^2 + z^2 - r[t]r[t]),
		"L" -> 5/3, "T" -> 1,
		"Default" -> {r[t] -> r0 (1 + A Sin[2 Pi t]), r0 -> 1, A -> 1/10}
	|>,
	"Torus" -> <|
		"dist" -> ({x, y, z, t} |-> Sqrt[(Sqrt[x^2 + y^2] - R)^2 + z^2] - r[x, y]),
		"phi" -> ({x, y, z, t} |-> (x^2 + y^2 + z^2 + R^2 - r[x, y]r[x, y])^2 - 4 R^2 (x^2 + y^2)),
		"L" -> 5/3, "T" -> 0,
		"Default" -> {r[x, y] -> r0 + (r1 - r0)(1 - x / Sqrt[x^2 + y^2]) / 2, R -> 1, r1 -> 1/2, r0 -> 1/2}
	|>
|>;


solutions = <|
	"PolynomialExact" -> <|
		"u_T" -> ({x, y, z, t} |-> P[{x, y, z}, t] . {-z^2, y, x}),
		"p" -> ({x, y, z, t} |-> x y^3 + z)
	|>,
	"KillingExact" -> <|
		"u_T" -> ({x, y, z, t} |-> \[Omega] {0, -z, y}),
		"p" -> ({x, y, z, t} |-> 0),
		"Default" -> {\[Omega] -> 1}
	|>,
	"DirectionChangeExact" -> <|
		"u_T" -> ({x, y, z, t} |-> P[{x, y, z}, t] . {1 - 2t, 0, 0}),
		"p" -> ({x, y, z, t} |-> 0)
	|>
|>;


srfDef = Lookup[surfaces[surface], "Default", {}];
slnDef = Lookup[solutions[solution], "Default", {}];


L = surfaces[surface]["L"];
T = surfaces[surface]["T"];
d[{x_, y_, z_}, t_] := surfaces[surface, "dist"][x, y, z, t]
\[Phi][{x_, y_, z_}, t_] := surfaces[surface, "phi"][x, y, z, t]
t0 = RandomReal@T;
x0 = FindRoot[{\[Phi][{x, y, z}, t] //. srfDef /. t -> t0, 0, 0}, {{x, L}, {y, L}, {z, L}}];
Print["x0: ", x0, ", t0 = ", t0]
Print["d(x0, t0) = ", d[{x, y, z}, t] //. srfDef /. x0 /. t -> t0]


replaceCPP[exp_] := StringReplace[ToString@CForm@N[exp /. {Pi -> PI, r' -> rprime, \!\(\*SuperscriptBox[\(r\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, y] -> rgrad0, \!\(\*SuperscriptBox[\(r\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, y] -> rgrad1}], {
	"x" -> "x[0]", "y" -> "x[1]", "z" -> "x[2]", "prime" -> "_prime", "grad0" -> "_grad(x[0], x[1])[0]", "grad1" -> "_grad(x[0], x[1])[1]", 
	"-1.*" -> "-", "- 1.*" -> "- ",
	"PI" -> "M_PI",
	"Power" -> "pow", "Sqrt" -> "std::sqrt", 
	"Cos" -> "cos", "Sin" -> "sin",
	"Arccos" -> "Acos", "Acsin" -> "asin",
	"\[Nu]" -> "nu", "\[Omega]" -> "omega", "\[Delta]0" -> "delta_0"
}]
exportCPP[f_] := Module[{val = f[{x, y, z}, t]}, If[ListQ[val], replaceCPP /@ val, replaceCPP@val]]


Print["d = ", exportCPP@d]
Print["\[Phi] = ", exportCPP@\[Phi]]
uN[{x_,y_,z_},t_] := Evaluate@FullSimplify[-D[d[{x, y, z}, t], t]];
Print["u_N = ", exportCPP@uN]
n[{x_, y_, z_}, t_] := Evaluate@FullSimplify[Grad[d[{x,y,z}, t],{x, y, z}]];
n[{x, y, z}, t]
Print["n = ", exportCPP@n // MatrixForm]


P[{x_, y_, z_}, t_] := Evaluate@FullSimplify[IdentityMatrix[3] - Transpose@{n[{x,y,z},t]} . {n[{x,y,z},t]}];
surfGrad[{x_,y_,z_},t_] := (P[{x,y,z},t] . Grad[#,{x,y,z}])&;
surfVectGrad[{x_,y_,z_},t_] := (P[{x,y,z},t] . Grad[#,{x,y,z}] . P[{x,y,z},t])&;
Es[{x_,y_,z_},t_] := 1/2 (surfVectGrad[{x,y,z},t][#]+Transpose[surfVectGrad[{x,y,z},t][#]])&;
surfVectDiv[{x_,y_,z_},t_] := Tr[surfVectGrad[{x,y,z},t][#]]&;
surfMatDiv[{x_,y_,z_},t_] := {surfVectDiv[{x,y,z},t][#[[1]]],surfVectDiv[{x,y,z},t][#[[2]]],surfVectDiv[{x,y,z},t][#[[3]]]}&
surfCurl[{x_,y_,z_},t_] := surfVectDiv[{x,y,z},t][Cross[#,n[{x,y,z},t]]]&;
surfVectCurl[{x_,y_,z_},t_] := Cross[n[{x,y,z},t],surfGrad[{x,y,z},t][#]]&;
conv[w_,{x_,y_,z_},t_] := {
	w[{x,y,z},t][[1]]D[#[[1]],x]+w[{x,y,z},t][[2]]D[#[[1]],y]+w[{x,y,z},t][[3]]D[#[[1]],z],
	w[{x,y,z},t][[1]]D[#[[2]],x]+w[{x,y,z},t][[2]]D[#[[2]],y]+w[{x,y,z},t][[3]]D[#[[2]],z],
	w[{x,y,z},t][[1]]D[#[[3]],x]+w[{x,y,z},t][[2]]D[#[[3]],y]+w[{x,y,z},t][[3]]D[#[[3]],z]
}&;


uT[{x_, y_, z_}, t_] := solutions[solution, "u_T"][x, y, z, t]
p[{x_, y_, z_}, t_] := solutions[solution, "p"][x, y, z, t]
uT[{x, y, z}, t] - P[{x, y, z}, t] . uT[{x, y, z}, t] // Simplify


uNH[{x_, y_, z_}, t_] := Evaluate@FullSimplify[uN[{x, y, z}, t]Grad[n[{x, y, z}, t], {x, y, z}], \[Phi][{x, y, z}, t] == 0]
Print["uN * H = ", uNH[{x, y, z}, t] // MatrixForm]


mg[{x_, y_, z_}, t_] := -Evaluate@FullSimplify[surfVectDiv[{x, y, z}, t][uT[{x,y,z},t]] + Tr[uNH[{x, y, z}, t]], \[Phi][{x, y, z}, t] == 0 && r[t] > 0]
mg[{x, y, z}, t]
Print["-g = ", exportCPP[mg]]


fStokes[{x_,y_,z_}, t_] := Evaluate@FullSimplify[
	+ D[uT[{x, y, z}, t], t] 
	- \[Nu] P[{x, y, z}, t] . surfMatDiv[{x, y, z}, t][Es[{x, y, z}, t][uT[{x, y, z}, t]] + uNH[{x, y, z}, t]] 
	+ uNH[{x, y, z}, t] . uT[{x, y, z}, t] 
	- uN[{x, y, z}, t]surfGrad[{x, y, z}, t][uN[{x, y, z}, t]] 
	+ surfGrad[{x, y, z}, t][p[{x, y, z}, t]],
\[Phi][{x, y, z}, t] == 0 && r[t] > 0]


fStokes[{x, y, z}, t]
Print["(Stokes) f_T = ", exportCPP[fStokes] // MatrixForm]


w[{x_, y_, z_}, t_] := uT[{x, y, z}, t] + uN[{x, y, z}, t]n[{x, y, z}, t]
fConv[{x_, y_, z_}, t_] := Evaluate@FullSimplify[P[{x, y, z}, t] . conv[w, {x, y, z}, t][uT[{x, y, z}, t]], \[Phi][{x, y, z}, t] == 0 && r[t] > 0]


fConv[{x, y, z}, t]
Print["(Convective term) f_T += ", exportCPP[fConv] // MatrixForm]
