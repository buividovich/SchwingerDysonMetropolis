(* ::Package:: *)

Get["G:\\LAT\\sd_metropolis\\convergence_analysis.m"];
NoisyScalarDataExtrapolation[data_,col_,nord_,{exact0_,exact1_},fitfuncs_,arg_]:=Module[{dataMean,dataCorrelationMatrix,Res0,NFactor,xs,dataExtrapolation,dataExtrapolationCovariance,x,GrData,GrExtrapolation,GrExact,res,err},
	{dataMean,dataCorrelationMatrix}=XYCorrelatedDataStatistics[data,col,nord];
	NFactor=(exact0-1)/dataMean[[1]];
	dataMean = 1 + NFactor dataMean;
	dataCorrelationMatrix *= NFactor^2; 
	xs = Table[1/m,{m,1,nord}];
	Clear[arg];
	{dataExtrapolation,dataExtrapolationCovariance}=FitCorrelatedData[xs,dataMean,dataCorrelationMatrix,fitfuncs,arg];
	GrData = PlotWithErrors[xs,dataMean,Sqrt[Diagonal[dataCorrelationMatrix]],PlotRange->{{0.0,All},{0.0,All}},PlotStyle->{Red,PointSize[0.01]}];
	GrExtrapolation = Plot[Sum[dataExtrapolation[[i]]fitfuncs[[i]],{i,1,Length[fitfuncs]}],{arg,0,1},PlotRange->All,PlotStyle->{Green}];
	GrExact = ListPlot[{{1, exact0},{0.5,exact1}},PlotStyle->{RGBColor[1.0,0.7,0.7],PointSize[0.02]}, PlotRange->{{0.0,All},{0.0,All}}];
	{{dataExtrapolation[[1]],Sqrt[dataExtrapolationCovariance[[1,1]]]},NFactor,Show[GrExact,GrData,GrExtrapolation]}
];
PrintScalarDataExtrapolation[exdata_,label_]:=Module[{},
	Print[label," = ",exdata[[1,1]]," +/- ",exdata[[1,2]],", NFactor = ",exdata[[2]]];
];
ImportCorrelators[Prefix_,LS_,nord_,NFactor_,{exact0_,exact1_},fitfuncs_,arg_]:=Module[{xcs,GrS,m,r,GxyData,np,aGxy,dGxy,cGxy,ip,midpoints,newmidpoints},
	xcs=Table[x,{x,0,LS-1}];
	GrS={};
	midpoints={};
	For[m=0,m<nord,m++,
		r=m/(nord-1);
		GxyData=Import[Prefix<>"_o"<>ToString[m]<>".dat","Real32"];
		GxyData=Partition[GxyData,LS];
		np=Length[GxyData];
		{aGxy,dGxy}=1/np Sum[
			cGxy=Re[Fourier[GxyData[[ip]],FourierParameters->{1,1}]];
			AppendTo[midpoints,{m, cGxy[[Quotient[LS,2]+1]]}];
			{cGxy,cGxy^2},{ip,1,Length[GxyData]}];
		dGxy=Sqrt[(dGxy-aGxy^2)/(np-1)];
		aGxy=1 + NFactor aGxy;
		dGxy = NFactor dGxy;
		AppendTo[GrS,PlotWithErrors[xcs,aGxy,dGxy,PlotStyle->{RGBColor[r,0,1-r]},PlotRange->{0.0,1.1}]];
		If[m==0,AppendTo[GrS,ListPlot[Transpose[{xcs,exact0}],PlotStyle->{RGBColor[0.7+0.3 r,0.7,0.7+0.3(1-r)],Thickness[0.01]},PlotRange->{0.0,1.0},Joined->True]]];
		If[m==1,AppendTo[GrS,ListPlot[Transpose[{xcs,exact1}],PlotStyle->{RGBColor[0.7+0.3 r,0.7,0.7+0.3(1-r)],Thickness[0.01]},PlotRange->{0.0,1.0},Joined->True]]];
	];
    midpoints = Partition[midpoints,np];
    newmidpoints = Table[midpoints[[m,i]],{i,1,np},{m,1,nord}];
	newmidpoints=Flatten[newmidpoints,1];
	{Show[Reverse[GrS]],NoisyScalarDataExtrapolation[newmidpoints,2,nord,{exact0[[Quotient[LS,2]+1]],exact1[[Quotient[LS,2]+1]]},fitfuncs,arg]}
];
