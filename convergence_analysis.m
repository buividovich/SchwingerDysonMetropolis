(* ::Package:: *)

$HistoryLength=0;
GetConvergenceData[fname_,ocol_,dcol_,add_,norm_,maxord_]:=Module[{oX,oX2,np,FData,i, io,Res,mX,eX},
(* Exporting the data *)
oX=Table[0,{io,1,maxord}];
oX2=Table[0, {i0,1,maxord}];
np=Table[0, {i0,1,maxord}];
FData=Import[fname,"Table"];
(*oXh=Table[{},{io,1,MaxOrder}];*)
For[i=1,i<=Length[FData],i++,
io=FData[[i,ocol]]+1;
np[[io]]++;
oX[[io]]+=FData[[i,dcol]];
oX2[[io]]+=FData[[i,dcol]]^2;
];
Res={};
For[io=1,io<=maxord,io++,
If[np[[io]]>1,
mX=oX[[io]]/np[[io]];
eX=Sqrt[(oX2[[io]]/np[[io]]-mX^2)/(np[[io]]-1)];
AppendTo[Res,{N[1/io],add +norm mX,norm eX}];
];
];
Res
];
LinearTransformColumn[t_,c_,a_,b_]:=Module[{st},
st=Transpose[t];
st[[c]]=a+b  st[[c]];
Transpose[st]
];
(* Useful functions for plotting *)
Needs["ErrorBarPlots`"];
PlotWithErrors[data_,col_]:=Module[{DataToPlot,i},
DataToPlot=Table[{{data[[i,1]],data[[i,2]]},ErrorBar[data[[i,3]]]},{i,1,Length[data]}];
ErrorListPlot[DataToPlot,PlotStyle->{col,PointSize[0.01]},PlotRange->{{0.0,All},All}]
];
(* Useful functions for formatting *)
fsndig[x_,n_]:=ToString[PaddedForm[x,{n+2,n},NumberPadding->{"","0"}]];
(* Column statistics *)
PrintStatistics[data_,label_]:=Module[{mean,error},
mean=Mean[data];
error=StandardDeviation[data]/Sqrt[Length[data]-1];
Print[label," = ",mean," +/- ",error];
];
