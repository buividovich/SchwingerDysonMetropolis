(* ::Package:: *)

(* Systematization of the low-order results *)
ExactCorrelators[LT_,LS_,\[Lambda]_]:=Module[{P2,G0,\[CapitalSigma]0,\[CapitalSigma]1,G1,m0,m1,\[Phi]x\[Phi]y0,\[Phi]x\[Phi]y1,Gx0,Gx1,TimeCorrelator0,TimeCorrelator1,SpaceCorrelator0,SpaceCorrelator1},
	P2[p0_,p1_]:=4*Sin[p0/2]^2 + 4*Sin[p1/2]^2;
	G0=Table[1/(\[Lambda]/4+P2[(2*\[Pi]*m0)/LT,(2*\[Pi]*m1)/LS]),{m0,0,LT-1},{m1,0,LS-1}];
	\[CapitalSigma]0=1/(LT LS) Total[G0,2];
	\[CapitalSigma]1[p0_,p1_]:=-2*(1 + 0.25*P2[p0,p1]*(1+(4 - \[Lambda]/4)\[CapitalSigma]0));
	G1=Table[G0[[m0+1,m1+1]]^2 \[CapitalSigma]1[(2*\[Pi]*m0)/LT,(2*\[Pi]*m1)/LS],{m0,0,LT-1},{m1,0,LS-1}];
	\[Phi]x\[Phi]y0=1/(LT LS) Re[Fourier[G0,FourierParameters->{1,1}]];
	\[Phi]x\[Phi]y1=1/(LT LS) Re[Fourier[G1,FourierParameters->{1,1}]];
	Gx0=1+2(-(\[Lambda]/8))\[CapitalSigma]0;
	Gx1= 1+2*(-\[Lambda]/8)*(\[CapitalSigma]0-\[Lambda]/8*\[Phi]x\[Phi]y1[[0+1,0+1]])+2*((-\[Lambda]/8)^2)*2*\[CapitalSigma]0^2;
	TimeCorrelator0=-1 + 2*Gx0 +  4*(-\[Lambda]/8)*(-1)*Table[\[Phi]x\[Phi]y0[[t+1,0+1]],{t,0,LT-1}];
	SpaceCorrelator0=-1 + 2*Gx0 +  4*(-\[Lambda]/8)*(-1)*Table[\[Phi]x\[Phi]y0[[0+1,x+1]],{x,0,LS-1}];
	TimeCorrelator1=-1 + 2*Gx1 + Table[ 4*(-\[Lambda]/8)*(-1)*(\[Phi]x\[Phi]y0[[t+1,0+1]]+(-\[Lambda]/8)*\[Phi]x\[Phi]y1[[t+1,0+1]])+4*(-\[Lambda]/8)^2 (-4*\[Phi]x\[Phi]y0[[t+1,0+1]]*\[CapitalSigma]0 + \[Phi]x\[Phi]y0[[t+1,0+1]]^2 + \[CapitalSigma]0^2),{t,0,LT-1}];
	SpaceCorrelator1=-1 + 2*Gx1 + Table[4*(-\[Lambda]/8)*(-1)*(\[Phi]x\[Phi]y0[[0+1,x+1]]+(-\[Lambda]/8)*\[Phi]x\[Phi]y1[[0+1,x+1]])+4*(-\[Lambda]/8)^2 (-4*\[Phi]x\[Phi]y0[[0+1,x+1]]*\[CapitalSigma]0 + \[Phi]x\[Phi]y0[[0+1,x+1]]^2 + \[CapitalSigma]0^2),{x,0,LS-1}];
	{{Gx0,Gx1},
	 {TimeCorrelator0[[1+1]],TimeCorrelator1[[1+1]]},
	 {SpaceCorrelator0[[1+1]],SpaceCorrelator1[[1+1]]},
	 {TimeCorrelator0,TimeCorrelator1},
	 {SpaceCorrelator0,SpaceCorrelator1}
	}
];
