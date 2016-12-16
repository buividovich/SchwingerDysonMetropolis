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
(* Numerical results of Vicari and Rossi from hep-lat/9307014 *)
VicariRossi\[Lambda]s={5.0,4.0,3.57,3.45,3.33,3.23,3.10,3.0120};
VicariRossiFiniteNData={
{{1/6,0.77592},{1/9,0.78087},{1/15,0.781395}},
{{1/6,0.68234},{1/9,0.70339},{1/15,0.70846}},
{{1/6,0.58690},{1/9,0.62920},{1/15,0.64990}},
{{1/6,0.54730},{1/9,0.58801},{1/15,0.62134}},
{{1/6,0.51139},{1/9,0.53847},{1/15,0.56809}},
{{1/6,0.48186},{1/9,0.50035},{1/15,0.51202}},
{{1/6,0.45209},{1/9,0.46607}},
{{1/6,0.43327},{1/9,0.44558}}
};
VicariRossiN6Data =  {{5.`,0.22407999999999995`},{4.`,0.31766000000000005`},{3.57`,0.4131`},{3.45`,0.4527`},{3.33`,0.48861`},{3.23`,0.51814`},{3.012`,0.56673`}};
VicariRossiN9Data =  {{5.`,0.21913000000000005`},{4.`,0.29661000000000004`},{3.57`,0.3708`},{3.45`,0.41198999999999997`},{3.33`,0.46153`},{3.23`,0.49965000000000004`},{3.012`,0.55442`}};
VicariRossiN15Data = {{5.`,0.21860500000000005`},{4.`,0.29154`},{3.57`,0.35009999999999997`},{3.45`,0.37866`},{3.33`,0.43191`},{3.23`,0.48797999999999997`}};
VicariRossiExtrapolatedData = {{5.0,0.21781750000000022},{4.0,0.28393499999999994}, {3.57,0.31904999999999983},{3.45,0.3286650000000003}, {3.33,0.38748000000000005}, {3.23,0.470475},{3.012,0.5298}};
(* Importing Semen's numerical results for finite temperature *)
Get["G:\\LAT\\sd_metropolis\\convergence_analysis.m"];
N6Data=Import["G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace\\N6_normal_L108_correlator.txt","Table"];
N9Data=Import["G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace\\N9_normal_L108_correlator.txt","Table"];
MCLTS=DeleteDuplicates[Join[Transpose[N6Data][[1]],Transpose[N9Data][[1]]]];
MCLTS=Sort[MCLTS];
N6Data=Partition[N6Data,108];
N9Data=Partition[N9Data,108];
MCCorrelatorPlots=Table[{},{LT,MCLTS}];
GeneratePlots[data_,color_]:=Module[{MyData,Gr1,Gr2,LT,ilt},
For[it=1,it<=Length[data],it++,
LT=data[[it]][[1,1]];
ilt=Position[MCLTS,LT][[1,1]];
MyData=Drop[data[[it]],0,1];
PrependTo[MyData,{0,1,0}];
(*Gr1=PlotWithErrorsTranspose[MyData,PlotRange->All,PlotStyle->{color}];*)
Gr2=ListPlot[Drop[MyData,0,-1 ],PlotRange->All,Joined->True,PlotStyle->{color,Opacity[0.3],Thickness[0.01]}];
AppendTo[MCCorrelatorPlots[[ilt]],Gr2];
(*AppendTo[MCCorrelatorPlots[[ilt]],Gr1];*)
];
];
GeneratePlots[N6Data,Green];
GeneratePlots[N9Data,Blue];
