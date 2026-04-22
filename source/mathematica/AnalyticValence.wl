(* ::Package:: *)

BeginPackage["AnalyticValence`"];

r1Results::usage = "r1Results[J,II,V,eps] returns {\[Phi], P, N}, the matched composite solutions for the potential, cation concentration and anion concentration for a binary electrolyte with valence ratio r=1.";
r2Results::usage = "r2Results[J,II,V,eps] returns {\[Phi], P, N}, the matched composite solutions for the potential, cation concentration and anion concentration for a binary electrolyte with valence ratio r=2.";
r12Results::usage = "r12Results[J,II,V,eps] returns {\[Phi], P, N}, the matched composite solutions for the potential, cation concentration and anion concentration for a binary electrolyte with valence ratio r=1/2.";

Begin["`Private`"];
(*r=1 electrolyte*)
symmetricParameters[J_,II_,V_]:=Module[{I,L,B,\[Gamma]2,nR,nL,\[Alpha]l,r,\[Alpha]r,zeroJ},
r=1;
zeroJ = J == 0;
\[Gamma]2:= If[zeroJ, Exp[V - II], Exp[V] ((1 - J)/(1 + J))^((II - 2 J)/(2 J))];
nR=(Sqrt[\[Gamma]2]*(1-J)+J+1)/(\[Gamma]2+Sqrt[\[Gamma]2]);
nL=\[Gamma]2*nR;
\[Alpha]r=Sqrt[(1-J)/nR];
\[Alpha]l=Sqrt[(1+J)/nL];
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}];

phiL[xL_,J_,II_,V_]:=Module[{I,L,B,NN,nR,nL,\[Alpha]l,\[Alpha]r,arg},{NN,nR,nL,\[Alpha]l,\[Alpha]r}=symmetricParameters[J,II,V];
arg=\[Alpha]l*Sqrt[(nL)/2]*xL;
V+(2)*Log[(\[Alpha]l*(1+\[Alpha]l*Tanh[arg]))/(\[Alpha]l+Tanh[arg])]];


phiR[xR_,J_,II_,V_]:=Module[{I,L,B,NN,nR,nL,\[Alpha]l,\[Alpha]r,arg},
{NN,nR,nL,\[Alpha]l,\[Alpha]r}=symmetricParameters[J,II,V];
arg=\[Alpha]r*Sqrt[(nR)/2]*xR;
(2)*Log[(\[Alpha]r*(1+\[Alpha]r*Tanh[arg]))/(\[Alpha]r+Tanh[arg])]];


phiOuter[x_,J_,II_,V_]:=Module[{I,zeroJ,L,B,NN,nR,nL,\[Alpha]l,\[Alpha]r},
zeroJ = J == 0;
{NN,nR,nL,\[Alpha]l,\[Alpha]r}=symmetricParameters[J,II,V];
If[zeroJ, V - II x - Log[nL], II/(2 J) Log[1 + J - 2 J x] + Log[(1 + J)/nL] + V - II/(2 J) Log[1 + J]]];

(*----------Composite----------*)
potentialzz[x_,J_,II_,V_,eps_]:=Module[{phi0,phi1,xL,xR},
phi0=phiOuter[0,J,II,V];
phi1=phiOuter[1,J,II,V];
xL=x/eps;
xR=(1-x)/eps;
phiOuter[x,J,II,V]+phiR[xR,J,II,V]+phiL[xL,J,II,V]-phi0-phi1];

pzz[x_,J_,II_,V_,eps_]:=Module[{I,L,B,\[Gamma]2,nR,p,nL,\[Alpha]l,\[Alpha]r,r,C,c0,c1,xL,xR,pL,pR},{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=symmetricParameters[J,II,V];
r=1;
C[p_]:=1-2*J*(p-0.5);
c0=C[0];
c1=C[1];
xL=x/eps;
xR=(1-x)/eps;
pL=(1+J)^(1+r)/nL^r;
pR=(1-J)^(1+r)/nR^r;
C[x]+pL*Exp[r*(V-phiL[xL,J,II,V])]+pR*Exp[-r*phiR[xR,J,II,V]]-c0-c1];

nzz[x_,J_,II_,V_,eps_]:=Module[{I,L,B,NN,nR,nL,\[Alpha]l,\[Alpha]r,C,c0,c1,xL,xR},{NN,nR,nL,\[Alpha]l,\[Alpha]r}=symmetricParameters[J,II,V];
C[p_]:=1-2*J*(p-0.5);
c0=C[0];
c1=C[1];
xL=x/eps;
xR=(1-x)/eps;
C[x]+nL*Exp[-(V-phiL[xL,J,II,V])]+nR*Exp[phiR[xR,J,II,V]]-c0-c1];

r1Results[J_,II_,V_,eps_]:={potentialzz[#,J,II,V,eps]&,pzz[#,J,II,V,eps]&,nzz[#,J,II,V,eps]&};

(*r=2 electrolyte*)

Parameters21[J_,II_,V_]:=Module[{I,L,B,\[Gamma]2,eq,sol,nR,nL,\[Alpha]l,\[Alpha]r,zeroJ},
zeroJ = J == 0;
\[Gamma]2=If[zeroJ, Exp[V - II], Exp[V] ((1 - J)/(1 + J))^((II - 2 J)/(2 J))];
eq=\[Gamma]2*nR*Sqrt[2*\[Gamma]2*nR+1+J]+\[Gamma]2*nR*Sqrt[2*nR+1-J]-\[Gamma]2*(1-J)*Sqrt[2*nR+1-J]-(1+J)*Sqrt[2*\[Gamma]2*nR+1+J]==0;sol=FindRoot[eq,{nR,0.5},MaxIterations->200];
nR=nR/.sol;(*store numerical value*)
nL=\[Gamma]2*nR;
\[Alpha]r=Sqrt[(1-J)/nR];
\[Alpha]l=Sqrt[(1+J)/nL];
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}];

phiL21[xL_,J_,II_,V_]:=Module[{sqrtTerm,\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r,arg},
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=Parameters21[J,II,V];
arg=\[Alpha]l*Sqrt[3*nL]/2*xL;
sqrtTerm=Sqrt[2/(3*\[Alpha]l^2)+1/3];V+Log[\[Alpha]l^2/2]+Log[3*((Tanh[arg]+sqrtTerm)/(1+sqrtTerm*Tanh[arg]))^2-1]];

phiR21[xR_,J_,II_,V_]:=Module[{I,L,B,\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r,arg,sqrtTermR},
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=Parameters21[J,II,V];
arg=\[Alpha]r*Sqrt[3*nR]/2*xR;sqrtTermR=Sqrt[2/(3*\[Alpha]r^2)+1/3];
Log[\[Alpha]r^2/2]+Log[3*((Tanh[arg]+sqrtTermR)/(1+sqrtTermR*Tanh[arg]))^2-1]];

phiOuter21[x_,J_,II_,V_]:=Module[{I,r,L,B,\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r,zeroJ},
zeroJ = J == 0;
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=Parameters21[J,II,V];
If[zeroJ, V - II x - Log[nL], II/(2 J) Log[1 + J - 2 J x] + Log[(1 + J)/nL] + V - II/(2 J) Log[1 + J]]];

potential21[x_,J_,II_,V_,eps_]:=Module[{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r,phi0,phi1,xL,xR},
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=Parameters21[J,II,V];
phi0=phiOuter21[0,J,II,V];
phi1=phiOuter21[1,J,II,V];
xL=x/eps;
xR=(1-x)/eps;
phiOuter21[x,J,II,V]+phiR21[xR,J,II,V]+phiL21[xL,J,II,V]-phi0-phi1];

p21[x_,J_,II_,V_,eps_]:=Module[{I,L,B,\[Gamma]2,nR,p,nL,\[Alpha]l,r,\[Alpha]r,C,c0,c1,xL,xR,pL,pR},{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=Parameters21[J,II,V];
r=2;
C[p_]:=1-2*J*(p-0.5);
c0=C[0];
c1=C[1];
xL=x/eps;
xR=(1-x)/eps;
pL=(1+J)^(1+r)/nL^r;
pR=(1-J)^(1+r)/nR^r;
C[x]+pL*Exp[r*(V-phiL21[xL,J,II,V])]+pR*Exp[-r*phiR21[xR,J,II,V]]-c0-c1];

n21[x_,J_,II_,V_,eps_]:=Module[{I,L,B,\[Gamma]2,nR,nL,r,\[Alpha]l,\[Alpha]r,C,c0,c1,xL,xR},
{\[Gamma]2,nR,nL,\[Alpha]l,\[Alpha]r}=Parameters21[J,II,V];
r=2;
C[p_]:=1-2*J*(p-0.5);
c0=C[0];
c1=C[1];
xL=x/eps;
xR=(1-x)/eps;
C[x]+nL*Exp[-(V-phiL21[xL,J,II,V])]+nR*Exp[phiR21[xR,J,II,V]]-c0-c1];

r2Results[J_,II_,V_,eps_]:={potential21[#,J,II,V,eps]&,p21[#,J,II,V,eps]&,n21[#,J,II,V,eps]&};

(*r=1/2 electrolyte*)

Parameters12[J_,II_,V_]:=Module[{H,\[Gamma],eq,sol,nR,nL,\[Alpha]l,\[Alpha]r,zeroJ},
zeroJ = J == 0;
\[Gamma]=If[zeroJ,Exp[(V - II)/2],Exp[V/2]*((1-J)/(1+J))^(II/(4*J)-1/2)]; (*Note this is \[Gamma], not \[Gamma]^2*)
\[Alpha]r=Sqrt[(1-J)/nR];
\[Alpha]l=Sqrt[(1+J)/(\[Gamma]^2*nR)];
H[a_]:=Sqrt[1+2 a];
eq=Sqrt[3]*(-1-\[Alpha]r+2 \[Alpha]r^2+(-1-\[Alpha]l+2 \[Alpha]l^2)*\[Gamma])+(-1-\[Alpha]r+2 \[Alpha]r^2-3 \[Gamma]+3 \[Alpha]l \[Gamma])*H[\[Alpha]l]+(-3+3 \[Alpha]r+(-1-\[Alpha]l+2 \[Alpha]l^2)*\[Gamma])*H[\[Alpha]r]+Sqrt[3]*(-1+\[Alpha]r+(-1+\[Alpha]l)*\[Gamma])*H[\[Alpha]l]*H[\[Alpha]r]==0;
sol=FindRoot[eq,{nR,0.5},MaxIterations->200];
nR=nR/.sol;(*store numerical value*)
nL=\[Gamma]^2*nR;
{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}];

phiL12[xL_,J_,II_,V_]:=Module[{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r,arg,sqrtTerm,numerator,denominator,frac},{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}=Parameters12[J,II,V];
arg=(\[Alpha]l*Sqrt[6 nL]*xL)/4;
sqrtTerm=Sqrt[3/(1+2 \[Alpha]l)];
numerator=2 \[Alpha]l*(Tanh[arg]+sqrtTerm)^2;
denominator=3*(1+sqrtTerm*Tanh[arg])^2-(Tanh[arg]+sqrtTerm)^2;
frac=numerator/denominator;
V+2*Log[frac]];

phiR12[xR_,J_,II_,V_]:=Module[{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r,arg,sqrtTerm,numerator,denominator,frac},{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}=Parameters12[J,II,V];
arg=(\[Alpha]r*Sqrt[6 nR]*xR)/4;
sqrtTerm=Sqrt[3/(1+2 \[Alpha]r)];
numerator=2 \[Alpha]r*(Tanh[arg]+sqrtTerm)^2;
denominator=3*(1+sqrtTerm*Tanh[arg])^2-(Tanh[arg]+sqrtTerm)^2;
frac=numerator/denominator;
2*Log[frac]];

phiOuter12[x_,J_,II_,V_]:=Module[{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r,zeroJ},
zeroJ = J == 0;
{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}=Parameters12[J,II,V];
If[zeroJ, V - II x - Log[nL], II/(2 J) Log[1 + J - 2 J x] + Log[(1 + J)/nL] + V - II/(2 J) Log[1 + J]]];

potential12[x_,J_,II_,V_,eps_]:=Module[{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r,phi0,phi1,xL,xR},
{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}=Parameters12[J,II,V];
phi0=phiOuter12[0,J,II,V];
phi1=phiOuter12[1,J,II,V];
xL=x/eps;
xR=(1-x)/eps;
phiOuter12[x,J,II,V]+phiR12[xR,J,II,V]+phiL12[xL,J,II,V]-phi0-phi1];

p12[x_,J_,II_,V_,eps_]:=Module[{\[Gamma],nR,s,nL,\[Alpha]l,r,\[Alpha]r,C,c0,c1,xL,xR,pL,pR},
{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}=Parameters12[J,II,V];
r=1/2;
C[s_]:=1-2*J*(s-0.5);
c0=C[0];
c1=C[1];
xL=x/eps;
xR=(1-x)/eps;
pL=(1+J)^(1+r)/nL^r;
pR=(1-J)^(1+r)/nR^r;
C[x]+pL*Exp[r*(V-phiL12[xL,J,II,V])]+pR*Exp[-r*phiR12[xR,J,II,V]]-c0-c1];

n12[x_,J_,II_,V_,eps_]:=Module[{I,L,B,\[Gamma],nR,nL,r,\[Alpha]l,\[Alpha]r,C,c0,c1,xL,xR,s},
{\[Gamma],nR,nL,\[Alpha]l,\[Alpha]r}=Parameters12[J,II,V];
r=1/2;
C[s_]:=1-2*J*(s-0.5);
c0=C[0];
c1=C[1];
xL=x/eps;
xR=(1-x)/eps;
C[x]+nL*Exp[-(V-phiL12[xL,J,II,V])]+nR*Exp[phiR12[xR,J,II,V]]-c0-c1];

r12Results[J_,II_,V_,eps_]:={potential12[#,J,II,V,eps]&,p12[#,J,II,V,eps]&,n12[#,J,II,V,eps]&};

End[];
EndPackage[];
