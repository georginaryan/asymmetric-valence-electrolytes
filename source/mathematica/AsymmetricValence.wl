(* ::Package:: *)

BeginPackage["AsymmetricValence`"];

AsymmetricValenceResults::usage = "AsymmetricValenceResults[r,J,II, V, \[Epsilon]] returns {\[Phi]Tot, PTot, NTot}, the matched composite solutions for the potential, cation concentration and anion concentration for a binary electrolyte with valence ratio r.";

Begin["`Private`"];

(*Solve the left boundary layer solution*)
LBL[nr_?NumericQ, JJ_?NumericQ, r_?NumericQ, II_?NumericQ, V_?NumericQ] := Module[{sol, m, k, \[Alpha]L, mFun, nlVal, zeroJ, NL, J},
(*Relation between nL and nR*)
J = Chop[JJ];
zeroJ = J == 0;
NL[k_] := If[zeroJ, k Exp[V - II], k Exp[V] ((1 - J)/(1 + J))^((II - 2 J)/(2 J))];
(*Set up parameters*)
nlVal = NL[nr];
\[Alpha]L = Sqrt[(1 + J)/nlVal];
(*Solve LBL differential equation*)
sol = Quiet[NDSolve[{m'[x] == Sign[\[Alpha]L - 1]*Sqrt[nlVal/2]*Sqrt[\[Alpha]L^(2 (1 + r))/r*m[x]^(2 - 2 r) + m[x]^4 - (1 + r)/r*\[Alpha]L^2*m[x]^2], m[0] == 1}, m, {x, 0, 40}, MaxSteps -> Infinity, Method -> "StiffnessSwitching", AccuracyGoal -> 8, PrecisionGoal -> 8, WorkingPrecision -> 12, MaxStepSize -> 0.3], NDSolve::precw];mFun = m /. First[sol];
(*return the InterpolatingFunction for m[x]*)
mFun];

(*Solve the right boundary layer solution for the potential*)
RBL[nr_?NumericQ, JJ_?NumericQ, r_?NumericQ, II_?NumericQ, V_?NumericQ] := Module[{sol, m, \[Alpha]R, mFun},
\[Alpha]R = Sqrt[(1 - JJ)/nr];(*force evaluation of nr to a number*)
sol = Quiet[NDSolve[{m'[x] == Sign[\[Alpha]R - 1]*Sqrt[nr/2]*Sqrt[\[Alpha]R^(2 (1 + r))/r*m[x]^(2 - 2 r) + m[x]^4 - (1 + r)/r*\[Alpha]R^2*m[x]^2],m[0] == 1},m, {x, 0, 40},
    MaxSteps -> Infinity,Method -> "StiffnessSwitching",AccuracyGoal -> 8,PrecisionGoal -> 8,WorkingPrecision -> 12,MaxStepSize -> 0.3],NDSolve::precw];
mFun = m /. First[sol];
(*return the InterpolatingFunction*)
mFun];

AsymmetricValenceResults[r_?NumericQ, JJ_?NumericQ,II_?NumericQ, V_?NumericQ, \[Epsilon]_?NumericQ] := Module[{J, Jp, Jn, \[Kappa], x, n, nl, nr, pl, pr, nrSolValue, nlSolValue, plSolValue, prSolValue, c0, K1, \[Phi]0, \[Alpha], nrSolution, IntCondition, nIntData, IntConditionFunction, zeroJ, LBLSol, RBLSol, nrList, LBLFin, RBLFin, PL, PR, NL, PTot, NTot, \[Phi]Tot, t},
J = Chop[JJ];
zeroJ = J == 0;
\[Kappa] = r + 1;
(*Boundary concentrations*)
NL[nr_] := If[zeroJ, nr Exp[V - II], nr Exp[V] ((1 - J)/(1 + J))^((II - 2 J)/(2 J))];
PL[nr_] = ((1 + J)^\[Kappa]/NL[nr]^(r));
PR[nr_] = ((1 - J)^\[Kappa]/nr^(r));
(*Outer Solution*)
c0[x_] = 1 - 2 J*(x - 1/2);
\[Phi]0[x_, nl_] := If[zeroJ, V - II x - Log[nl], II/(2 J) Log[1 + J - 2 J x] + Log[(1 + J)/nl] + V - II/(2 J) Log[1 + J]];
(*Iterate over integral condition Eqn 50*)
nrList = Range[0.05, 3.5, 0.05];
IntCondition = Table[Module[{result}, LBLSol = LBL[nr, J, r, II, V];
RBLSol = RBL[nr, J, r, II, V];
result = Quiet[NIntegrate[Re[PL[nr]*(LBLSol[t])^(-2 r) - NL[nr]*(LBLSol[t])^2],{t, 0, 20}] +NIntegrate[Re[PR[nr]*(RBLSol[t])^(-2 r) - nr*(RBLSol[t])^2],{t, 0, 20}],NIntegrate::precw];{nr, result}], {nr, nrList}];
(*Create interpolation function of integral condition*)
IntConditionFunction = Interpolation[IntCondition, InterpolationOrder -> 5];
(*Solve for nR*)
nrSolution = FindRoot[IntConditionFunction[n] == 0, {n, 0.5}, AccuracyGoal -> 8, PrecisionGoal -> 10];
nrSolValue = n /. nrSolution;
(*Final boundary concentrations*)
prSolValue = PR[nrSolValue];
nlSolValue = NL[nrSolValue];
plSolValue = PL[nrSolValue];
(*Final Boundary layer solutions*)
LBLFin = LBL[nrSolValue, J, r, II, V];
RBLFin = RBL[nrSolValue, J, r, II, V];
(*Final composite matching solutions*)
\[Phi]Tot[t_] = \[Phi]0[t, nlSolValue] + 2*Log[RBLFin[(1 - t)/\[Epsilon]]] - \[Phi]0[1, nlSolValue] + (V + 2*Log[LBLFin[t/\[Epsilon]]] - \[Phi]0[0, nlSolValue]);
PTot[t_] = plSolValue*(LBLFin[t/\[Epsilon]])^(-2 r) - c0[0] + c0[t] + prSolValue*(RBLFin[(1 - t)/\[Epsilon]])^(-2 r) - c0[1];
NTot[t_] = nlSolValue*(LBLFin[t/\[Epsilon]])^(2) - c0[0] + c0[t] + nrSolValue*(RBLFin[(1 - t)/\[Epsilon]])^(2) - c0[1];
{\[Phi]Tot, PTot, NTot}];

End[];
EndPackage[];
