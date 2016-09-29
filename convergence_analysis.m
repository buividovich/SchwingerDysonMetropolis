(* ::Package:: *)

$HistoryLength=0;
(* Useful function for formatting - double number with n digits after the comma *)
fsndig[x_,n_]:=ToString[PaddedForm[x,{n+2,n},NumberPadding->{"","0"}]];
(* Assuming that X are integers from 0 to xmax-1, calculate vevs and autocorrs *)
XYCorrelatedDataStatistics::wrongdata = "Wrong number of data points provided: `1`";
XYCorrelatedDataStatistics[data_,ycol_,xmax_]:=Module[{np,aX,cXX,ip,x,x1,x2},
	np=Quotient[Length[data],xmax]; 
	If[Mod[Length[data],xmax]!=0,Message[XYCorrelatedDataStatistics::wrongdata,Length[data]]];
	aX = 1/np Sum[Table[data[[ip xmax + x+1,ycol]],{x,0,xmax-1}],{ip,0,np-1}];
    cXX= 1/np Sum[Table[data[[ip xmax + x1+1,ycol]]data[[ip xmax + x2+1,ycol]],{x1,0,xmax-1},{x2,0,xmax-1}],{ip,0,np-1}] - KroneckerProduct[aX,aX];
	{aX,cXX/(np-1)}
];
(* Fitting with autocorrelation matrix taken into account *)
FitCorrelatedData[x_,y_,cyy_,funcs_,arg_]:=Module[{nf,nd,pars,cyyi,y0,i,j,m,m1,m2,mFit,pcm},
    nf = Length[funcs];
    nd = Min[Length[x],Length[y],Min[Dimensions[cyy]]];
	pars = Table[Symbol["p"<>ToString[i]],{i,1,nf}];
    cyyi = Inverse[cyy];
	y0   = Table[Sum[pars[[i]](funcs[[i]]/.{arg->x[[m]]}),{i,1,nf}],{m,1,nd}];
	mFit = NMinimize[(y - y0).cyyi.(y-y0),pars];
    pcm  = Table[Sum[(funcs[[i]]/.{arg->x[[m1]]})cyyi[[m1,m2]](funcs[[j]]/.{arg->x[[m2]]}),{m1,1,nd},{m2,1,nd}],{i,1,nf},{j,1,nf}];
	{pars/.mFit[[2]],Inverse[pcm]}
];
(* Useful functions for plotting *)
Needs["ErrorBarPlots`"];
PlotWithErrors::wrongdata = "Lengths of x, y, and e columns not equal `1` `2` `3`";
PlotWithErrors[x_,y_,e_,opts:OptionsPattern[]]:=Module[{DataToPlot,MinLength,i},
	If[!(Length[x]==Length[y] && Length[y]==Length[e]),
		Message[PlotWithErrors::wrongdata, Length[x],Length[y],Length[e]];
	];
	MinLength = Min[Min[Length[x],Length[y]],Length[e]];
	DataToPlot=Table[{{x[[i]],y[[i]]},ErrorBar[e[[i]]]},{i,1,Length[x]}];
	ErrorListPlot[DataToPlot,Evaluate[FilterRules[{opts},Options[ErrorListPlot]]]]
];
(* Column statistics *)
PrintStatistics[data_,label_]:=Module[{mean,error},
mean=Mean[data];
error=StandardDeviation[data]/Sqrt[Length[data]-1];
Print[label," = ",mean," +/- ",error];
];
