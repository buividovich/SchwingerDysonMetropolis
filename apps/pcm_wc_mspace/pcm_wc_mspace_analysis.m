(* ::Package:: *)

Get["G:\\LAT\\sd_metropolis\\convergence_analysis.m"];
NoisyScalarDataExtrapolation[data_,col_,nord_,{exact0_,exact1_},fitfuncs_,arg_,nptoomit_]:=Module[
{dataMean,dataCorrelationMatrix,dataMeanToFit,dataCorrelationMatrixToFit,Res0,NFactor,xs,xsToFit,dataExtrapolation,dataExtrapolationCovariance,x,GrData,GrExtrapolation,GrExact,res,err,ExactDataToPlot,mcpoint,i,m},
	{dataMean,dataCorrelationMatrix}=XYCorrelatedDataStatistics[data,col,nord];
	NFactor=(exact0-1)/dataMean[[1]];
	dataMean = 1 + NFactor dataMean;
	dataCorrelationMatrix *= NFactor^2; 
	dataMeanToFit = Drop[dataMean,nptoomit];
	dataCorrelationMatrixToFit = Drop[dataCorrelationMatrix,nptoomit,nptoomit];
	xs = Table[1/m,{m,1,nord}];
	xsToFit=Drop[xs,nptoomit];
	Clear[arg];
	{dataExtrapolation,dataExtrapolationCovariance}=FitCorrelatedData[xsToFit,dataMeanToFit,dataCorrelationMatrixToFit,fitfuncs,arg];
	GrData = PlotWithErrors[xs,dataMean,Sqrt[Diagonal[dataCorrelationMatrix]],PlotRange->{{0.0,All},{0.0,All}},PlotStyle->{Red,PointSize[0.01]}];
	GrExtrapolation = Plot[Sum[dataExtrapolation[[i]]fitfuncs[[i]],{i,1,Length[fitfuncs]}],{arg,0,1},PlotRange->All,PlotStyle->{Green}];
	ExactDataToPlot ={{1, exact0},{0.5,exact1}};
	GrExact = ListPlot[ExactDataToPlot,PlotStyle->{RGBColor[1.0,0.7,0.7],PointSize[0.02]}, PlotRange->{{0.0,All},{0.0,All}}];
	{{dataExtrapolation[[1]],Sqrt[dataExtrapolationCovariance[[1,1]]]},NFactor,Show[GrExact,GrData,GrExtrapolation]}
];
PrintScalarDataExtrapolation[exdata_,label_]:=Module[{},
	Print[label," = ",exdata[[1,1]]," +/- ",exdata[[1,2]],", NFactor = ",exdata[[2]]];
];
ImportCorrelators[Prefixes_,LS_,nord_,NFactor_,{exact0_,exact1_},fitfuncs_,arg_]:=Module[{x,xcs,GrS,m,i,r,GxyData,np,aGxy,dGxy,cGxy,ip,midpoints,newmidpoints,Prefix},
	xcs=Table[x,{x,0,LS-1}];
	GrS={};
	midpoints={};
	For[m=0,m<nord,m++,
		r=m/(nord-1);
		GxyData={};
		Do[
			If[FileExistsQ[Prefix<>"_o"<>ToString[m]<>".dat"],GxyData = Join[GxyData,Import[Prefix<>"_o"<>ToString[m]<>".dat","Real32"]];];,
		{Prefix,Prefixes}];
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
		If[m==0,AppendTo[GrS,ListPlot[Transpose[{xcs,exact0}],PlotStyle->{RGBColor[r,0,1-r],Opacity[0.5],Thickness[0.01]},PlotRange->{0.0,1.0},Joined->True]]];
		If[m==1,AppendTo[GrS,ListPlot[Transpose[{xcs,exact1}],PlotStyle->{RGBColor[r,0,1-r],Opacity[0.5],Thickness[0.01]},PlotRange->{0.0,1.0},Joined->True]]];
	];
    midpoints = Partition[midpoints,np];
    newmidpoints = Table[midpoints[[m,i]],{i,1,np},{m,1,nord}];
	newmidpoints=Flatten[newmidpoints,1];
	{Show[Reverse[GrS]],NoisyScalarDataExtrapolation[newmidpoints,2,nord,{exact0[[Quotient[LS,2]+1]],exact1[[Quotient[LS,2]+1]]},fitfuncs,arg,0]}
];
