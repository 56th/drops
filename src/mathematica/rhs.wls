#!/usr/bin/env wolframscript
(* ::Package:: *)

surface = "Sphere";
solution = "KillingExact";


surfaces = <|
	"Sphere" -> <|
		"dist" -> ({x, y, z, t} |-> Sqrt[x^2 + y^2 + z^2] - r[t]),
		"phi" -> ({x, y, z, t} |-> x^2 + y^2 + z^2 - r[t]r[t]),
		"Assumptions" -> r[t] > 0,
		"Default" -> {r[t] -> r0 (1 + A Sin[2 Pi t]), r0 -> 1, A -> 1/10},
		"T" -> 1
	|>,
	"Torus" -> <|
		"dist" -> ({x, y, z, t} |-> Sqrt[(Sqrt[x^2 + y^2] - R)^2 + z^2] - r[x, y]),
		"phi" -> ({x, y, z, t} |-> (x^2 + y^2 + z^2 + R^2 - r[x, y]r[x, y])^2 - 4 R^2 (x^2 + y^2)),
		"Assumptions" -> R > r[x, y] && r[x, y] > 0,
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
		"c" -> ({x, y, z, t} |-> 1/2 (1+Tanh[(z Cos[\[Omega] t]-y Sin[\[Omega] t])/(2 Sqrt@2 \[Epsilon])])),
		"Default" -> {\[Omega] -> 1, \[Epsilon] -> .01}
	|>,
	"DirectionChangeExact" -> <|
		"u_T" -> ({x, y, z, t} |-> P[{x, y, z}, t] . {1 - 2t, 0, 0})
	|>
|>;


srfDef = Lookup[surfaces[surface], "Default", {}];
slnDef = Lookup[solutions[solution], "Default", {}];
assumptions = Lookup[surfaces[surface], "Assumptions", True];


complexity[exp_] := StringLength@ToString@InputForm@exp
simplify[exp_, x_:True] := FullSimplify[exp, ComplexityFunction -> complexity, Assumptions -> assumptions && x]


L = Lookup[surfaces[surface], "L", 5/3];
T = Lookup[surfaces[surface], "T", 0];
d[{x_, y_, z_}, t_] := Evaluate[surfaces[surface, "dist"]][x, y, z, t]
\[Phi][{x_, y_, z_}, t_] := Evaluate[surfaces[surface, "phi"]][x, y, z, t]
t0 = RandomReal@T;
x0 = FindRoot[{\[Phi][{x, y, z}, t] //. srfDef /. t -> t0, 0, 0}, {{x, L}, {y, L}, {z, L}}];
Print["x0: ", x0, ", t0 = ", t0]
Print["d(x0, t0) = ", d[{x, y, z}, t] //. srfDef /. x0 /. t -> t0]


Print["d = ", d[{x, y, z}, t]]
Print["\[Phi] = ", \[Phi][{x, y, z}, t]]
uN[{x_,y_,z_},t_] := Evaluate@simplify[-D[d[{x, y, z}, t], t]];
Print["u_N = ", uN[{x, y, z}, t]]
n[{x_, y_, z_}, t_] := Evaluate@simplify[Grad[d[{x,y,z}, t],{x, y, z}]];
Print["n = ", MatrixForm@n[{x, y, z}, t]]


P[{x_, y_, z_}, t_] := Evaluate@simplify[IdentityMatrix[3] - Transpose@{n[{x,y,z},t]} . {n[{x,y,z},t]}, \[Phi][{x, y, z}, t] == 0];
Print["P = ", MatrixForm@P[{x, y, z}, t]]
surfGrad[{x_, y_, z_}, t_] := (P[{x, y, z}, t] . Grad[#, {x, y, z}])&;
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


uT[{x_, y_, z_}, t_] := Evaluate[Lookup[solutions[solution], "u_T", {x, y, z, t} |-> {0, 0, 0}]][x, y, z, t]
p[{x_, y_, z_}, t_] := Evaluate[Lookup[solutions[solution], "p", {x, y, z, t} |-> 0]][x, y, z, t]
c[{x_, y_, z_}, t_] := Evaluate[Lookup[solutions[solution], "c", {x, y, z, t} |-> 0]][x, y, z, t];
Print["|u_T - P u_T| = ", N@Norm@simplify[uT[{x, y, z}, t] - P[{x, y, z}, t] . uT[{x, y, z}, t], \[Phi][{x, y, z}, t] == 0]]


uNH[{x_, y_, z_}, t_] := Evaluate@simplify[uN[{x, y, z}, t]Grad[n[{x, y, z}, t], {x, y, z}], \[Phi][{x, y, z}, t] == 0]
Print["u_N H = ", MatrixForm@uNH[{x, y, z}, t]]


mg[{x_, y_, z_}, t_] := Evaluate@simplify[-(surfVectDiv[{x, y, z}, t][uT[{x,y,z},t]] + Tr[uNH[{x, y, z}, t]]), \[Phi][{x, y, z}, t] == 0]
Print["-g = ", mg[{x, y, z}, t]]


fStokes[{x_,y_,z_}, t_] := Evaluate@simplify[
	+ D[uT[{x, y, z}, t], t] 
	- \[Nu] P[{x, y, z}, t] . surfMatDiv[{x, y, z}, t][Es[{x, y, z}, t][uT[{x, y, z}, t]] + uNH[{x, y, z}, t]] 
	+ uNH[{x, y, z}, t] . uT[{x, y, z}, t] 
	- uN[{x, y, z}, t]surfGrad[{x, y, z}, t][uN[{x, y, z}, t]] 
	+ surfGrad[{x, y, z}, t][p[{x, y, z}, t]],
\[Phi][{x, y, z}, t] == 0]


Print["(Stokes) f_T = ", MatrixForm[fStokes[{x, y, z}, t]]]


w[{x_, y_, z_}, t_] := uT[{x, y, z}, t] + uN[{x, y, z}, t]n[{x, y, z}, t]
fConv[{x_, y_, z_}, t_] := Evaluate@simplify[P[{x, y, z}, t] . conv[w, {x, y, z}, t][uT[{x, y, z}, t]], \[Phi][{x, y, z}, t] == 0]


Print["(Convective term) f_T += ", MatrixForm[fConv[{x, y, z}, t]]]


f0[c_] := 1/4 c^2 (1 - c)^2
M[c_] := c (1 - c)
\[Mu][{x_, y_, z_}, t_] := Evaluate@simplify[f0'[c[{x, y, z}, t]] - \[Epsilon]^2 surfVectDiv[{x, y, z}, t][surfGrad[{x, y, z}, t][c[{x, y, z}, t]]], \[Phi][{x, y, z}, t] == 0]
Print["mu = ", \[Mu][{x, y, z}, t]]


simplify[D[c[{x, y, z}, t], t] + surfGrad[{x, y, z}, t][c[{x, y, z}, t]] . uT[{x, y, z}, t], \[Phi][{x, y, z}, t] == 0]
fC[{x_, y_, z_}, t_] := Evaluate@simplify[-surfVectDiv[{x, y, z}, t][M[c[{x, y, z}, t]]surfGrad[{x, y, z}, t][\[Mu][{x, y, z}, t]]], \[Phi][{x, y, z}, t] == 0]
Print["(CH) f = ", fC[{x, y, z}, t]]


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
	"Arccos" -> "Acos", "Acsin" -> "asin", "Tan" -> "tan", "Sec" -> "sec",
	"\[Nu]" -> "nu", "\[Omega]" -> "omega", "\[Epsilon]" -> "eps", "\[Delta]0" -> "delta_0"
}]
exportCPP[f_] := Module[{val = f[{x, y, z}, t]}, If[ListQ[val], replaceCPP /@ val, replaceCPP@val]]
