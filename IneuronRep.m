(* ::Package:: *)

IneuronRate[lambda_,mu_,p_,size_]:=Module[{theta,b,deltax,deltay,idxdiag,e,mulnormal,w,v,iemis,iemisy,iefull,rmis,rmisy,rfull,s,sign,stimVec,J0,J0p,side,chiral,idp,a0,thre,a,hImis,hImisy,hIfull,kappa,interp},
theta = 0;
b=150;
deltax[b_,x_,y_]:=((1+b/2) x-b/2 Sqrt[mu] y)/(1+b+(1-mu) b^2/4);
deltay[b_,x_,y_]:=((1+b/2) y-b/2 Sqrt[mu] x)/(1+b+(1-mu) b^2/4);
idxdiag = Table[{i,i},{i,1,size}];
e =Table[1,{i,1,size}];
(*Calculate the weight*)
mulnormal = RandomVariate[\[ScriptD]=MultinormalDistribution[{3,3},{{1,Sqrt[mu]},{Sqrt[mu],1}}],{p,size}];
w = Ramp[mulnormal[[;;,;;,1]]];
v = Ramp[mulnormal[[;;,;;,2]]];
(*Calculate E neurons' response*)
iemis=(w[[1]]-Mean[w[[1]]])*deltax[b,1,0]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,0];
iemisy=(w[[1]]-Mean[w[[1]]])*deltax[b,0,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,0,1];
iefull=(w[[1]]-Mean[w[[1]]])*deltax[b,1,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,1];
rmis=b*Ramp[iemis-theta];
rmisy=b*Ramp[iemisy-theta];
rfull=b*Ramp[iefull-theta];
s=(rmis-Mean[rmis])*(rfull-Mean[rfull]);

(*Calculate stimulus input vector and put in the matrix*)
stimVec = Join[Table[w[[i]],{i,1,p}],Table[v[[i]],{i,1,p}]];
J0 = IdentityMatrix[size];
side=stimVec[[1;;2 p,2p+1;;]];
J0p=stimVec[[1;;2 p,1;;2 p]];
(*The order*)
chiral = Det[J0p];
idp=IdentityMatrix[2 p];
If[chiral<= 0,
idp=ReplacePart[idp,{1->idp[[2]],2->idp[[1]]}]
];

(*Diagonal elements*);
thre=0.8-Log[0.5]/4;
a0=HeavisideTheta[-s](1.4+12Exp[1.5s])+0.2HeavisideTheta[s](0.01+HeavisideTheta[s-thre]);
kappa =lambda;
a =kappa*a0+(1-kappa)*e;
interp=kappa*J0p/size+(1-kappa)*idp;
J0=ReplacePart[J0,{i_,i_}:>a[[i]]];
J0[[1;;2 p,1;;2 p]] = interp;
J0[[1;;2 p,2p+1;;]] = kappa*side/size;
hImis = J0 . rmis;
hImisy = J0 . rmisy;
hIfull = J0 . rfull;
{hImis,hImisy,hIfull}
]


IneuronJie[lambda_,mu_,p_,size_]:=Module[{theta,b,deltax,deltay,idxdiag,e,mulnormal,w,v,Jout,iemisx,iemisy,iefull,rmisx,rmisy,rfull,s,sign,stimVec,J0p,Jie,Jei,side,chiral,idp,a0,thre,a,hImisx,hImisy,hIfull,alambda,kappa,interp},
theta = 0;
b=150;
deltax[b_,x_,y_]:=((1+b/2) x-b/2 Sqrt[mu] y)/(1+b+(1-mu) b^2/4);
deltay[b_,x_,y_]:=((1+b/2) y-b/2 Sqrt[mu] x)/(1+b+(1-mu) b^2/4);
idxdiag = Table[{i,i},{i,1,size}];
e =Table[1,{i,1,size}];
(*Calculate the weight*)
mulnormal = RandomVariate[\[ScriptD]=MultinormalDistribution[{3,3},{{1,Sqrt[mu]},{Sqrt[mu],1}}],{p,size}];
w = Ramp[mulnormal[[;;,;;,1]]];
v = Ramp[mulnormal[[;;,;;,2]]];
Jout = (w\[Transpose] . w+v\[Transpose] . v)/size;
(*Calculate E neurons' response*)
iemisx=(w[[1]]-Mean[w[[1]]])*deltax[b,1,0]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,0];
iemisy=(w[[1]]-Mean[w[[1]]])*deltax[b,0,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,0,1];
iefull=(w[[1]]-Mean[w[[1]]])*deltax[b,1,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,1];
rmisx=b*Ramp[iemisx-theta];
rmisy=b*Ramp[iemisy-theta];
rfull=b*Ramp[iefull-theta];
s=(rmisx-Mean[rmisx])*(rfull-Mean[rfull]);

(*Calculate stimulus input vector and put in the matrix*)
stimVec = Join[Table[w[[i]],{i,1,p}],Table[v[[i]],{i,1,p}]];
side=stimVec[[1;;2 p,2p+1;;]];
J0p=stimVec[[1;;2 p,1;;2 p]];
(*The order*)
chiral = Det[J0p];
idp=IdentityMatrix[2 p];
If[chiral<= 0,
idp=ReplacePart[idp,{1->idp[[2]],2->idp[[1]]}]
];

(*Diagonal elements*);
thre=0.8-Log[0.5]/4;
a0=HeavisideTheta[-s](1.4+12Exp[1.5s])+0.2HeavisideTheta[s](0.01+HeavisideTheta[s-thre]);
Jie = IdentityMatrix[size];
Jei = ConstantArray[0,{size,size}];
Jei[[All,1;;2p]]=stimVec\[Transpose];

kappa =lambda;
a =kappa*a0+(1-kappa)*e;
interp=kappa*J0p/size+(1-kappa)*idp;
Jie=ReplacePart[Jie,{i_,i_}:>a[[i]]];
Jie[[1;;2 p,1;;2 p]]=interp;
Jie[[1;;2 p,2p+1;;]] = kappa*side/size;
hImisx =Jie . rmisx;
hImisy =Jie . rmisy;
hIfull =Jie . rfull;
alambda = DiagonalMatrix[(1-kappa)/a];
Jei[[All,2p+1;;]]=Jout[[All,2p+1;;]] . alambda[[2p+1;;,2p+1;;]];
{{iemisx,iemisy,iefull,hImisx,hImisy,hIfull},Jie,Jei}
]


EImeanRate[mu_,p_,size_]:=Module[{theta,b,deltax,deltay,e,mulnormal,w,v,iemis,iemisy,iefull,rmis,rmisy,rfull,Jee,Ji,Jfull},
theta = 0;
b=150;
deltax[b_,x_,y_]:=((1+b/2) x-b/2 Sqrt[mu] y)/(1+b+(1-mu) b^2/4);
deltay[b_,x_,y_]:=((1+b/2) y-b/2 Sqrt[mu] x)/(1+b+(1-mu) b^2/4);
e =Table[1,{i,1,size}];
(*Calculate the weight*)
mulnormal = RandomVariate[\[ScriptD]=MultinormalDistribution[{3,3},{{1,Sqrt[mu]},{Sqrt[mu],1}}],{p,size}];
w = Ramp[mulnormal[[;;,;;,1]]];
v = Ramp[mulnormal[[;;,;;,2]]];
(*Calculate E neurons' response*)
iemis=(w[[1]]-Mean[w[[1]]])*deltax[b,1,0]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,0];
iemisy=(w[[1]]-Mean[w[[1]]])*deltax[b,0,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,0,1];
iefull=(w[[1]]-Mean[w[[1]]])*deltax[b,1,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,1];
rmis=b*Ramp[iemis-theta];
rmisy=b*Ramp[iemisy-theta];
rfull=b*Ramp[iefull-theta];
Jee=((w\[Transpose] . w+v\[Transpose] . v)+Mean[w[[1]]] (w\[Transpose] . {e}+{e}\[Transpose] . w)+Mean[v[[1]]](v\[Transpose] . {e}+{e}\[Transpose] . v))/size;
Ji=(2(w\[Transpose] . w+v\[Transpose] . v)+Mean[w[[1]]]^2+Mean[v[[1]]]^2)/size;
{Jee . rmis,Jee . rmisy,Jee . rfull,Ji . rmis,Ji . rmisy,Ji . rfull,iemis,iemisy,iefull}
]


IneuronJieEvolve1[lambda_,mu_,p_,size_,x0_,y0_,z0_]:=Module[{theta,b,deltax,deltay,idxdiag,e,mulnormal,w,v,wv,Jout,iemisx,iemisy,iefull,rmisx,rmisy,rfull,s,sign,stimVec,J0p,Jie,Jei,side,chiral,idp,a0,thre,a,hImisx,hImisy,hIfull,alambda,kappa,interp},
theta = 0;
b=150;
deltax[b_,x_,y_]:=((1+b/2) x-b/2 Sqrt[mu] y)/(1+b+(1-mu) b^2/4);
deltay[b_,x_,y_]:=((1+b/2) y-b/2 Sqrt[mu] x)/(1+b+(1-mu) b^2/4);
idxdiag = Table[{i,i},{i,1,size}];
e =Table[1,{i,1,size}];
(*Calculate the weight*)
wv=Table[{3+Sqrt[Sqrt[mu]]x0[[i]]+Sqrt[1-Sqrt[mu]]y0[[i]],3+Sqrt[Sqrt[mu]]x0[[i]]+Sqrt[1-Sqrt[mu]]z0[[i]]},{i,1,Length[x0]}];
w = {Ramp[wv[[All,1]]]};
v = {Ramp[wv[[All,2]]]};
Jout = (w\[Transpose] . w+v\[Transpose] . v)/size;
(*Calculate E neurons' response*)
iemisx=(w[[1]]-Mean[w[[1]]])*deltax[b,1,0]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,0];
iemisy=(w[[1]]-Mean[w[[1]]])*deltax[b,0,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,0,1];
iefull=(w[[1]]-Mean[w[[1]]])*deltax[b,1,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,1];
rmisx=b*Ramp[iemisx-theta];
rmisy=b*Ramp[iemisy-theta];
rfull=b*Ramp[iefull-theta];
s=(rmisx-Mean[rmisx])*(rfull-Mean[rfull]);

(*Calculate stimulus input vector and put in the matrix*)
stimVec = Join[Table[w[[i]],{i,1,p}],Table[v[[i]],{i,1,p}]];
side=stimVec[[1;;2 p,2p+1;;]];
J0p=stimVec[[1;;2 p,1;;2 p]];
(*The order*)
chiral = Det[J0p];
idp=IdentityMatrix[2 p];
If[chiral<= 0,
idp=ReplacePart[idp,{1->idp[[2]],2->idp[[1]]}]
];

(*Diagonal elements*);
thre=0.8-Log[0.5]/4;
a0=HeavisideTheta[-s](1.4+12Exp[1.5s])+0.2HeavisideTheta[s](0.01+HeavisideTheta[s-thre]);
Jie = IdentityMatrix[size];
Jei = ConstantArray[0,{size,size}];
Jei[[All,1;;2p]]=stimVec\[Transpose];

kappa =lambda;
a =kappa*a0+(1-kappa)*e;
interp=kappa*J0p/size+(1-kappa)*idp;
Jie=ReplacePart[Jie,{i_,i_}:>a[[i]]];
Jie[[1;;2 p,1;;2 p]]=interp;
Jie[[1;;2 p,2p+1;;]] = kappa*side/size;
hImisx =Jie . rmisx;
hImisy =Jie . rmisy;
hIfull =Jie . rfull;
alambda = DiagonalMatrix[(1-kappa)/a];
Jei[[All,2p+1;;]]=Jout[[All,2p+1;;]] . alambda[[2p+1;;,2p+1;;]];
{{iemisx,iemisy,iefull,hImisx,hImisy,hIfull},Jie,Jei}
]


IneuronJieEvolve2[lambda_,mu_,p_,size_,x0_,y0_]:=Module[{theta,b,deltax,deltay,idxdiag,e,mulnormal,w,v,wv,Jout,iemisx,iemisy,iefull,rmisx,rmisy,rfull,s,sign,stimVec,J0p,Jie,Jei,side,chiral,idp,a0,thre,a,hImisx,hImisy,hIfull,alambda,kappa,interp},
theta = 0;
b=150;
deltax[b_,x_,y_]:=((1+b/2) x-b/2 Sqrt[mu] y)/(1+b+(1-mu) b^2/4);
deltay[b_,x_,y_]:=((1+b/2) y-b/2 Sqrt[mu] x)/(1+b+(1-mu) b^2/4);
idxdiag = Table[{i,i},{i,1,size}];
e =Table[1,{i,1,size}];
(*Calculate the weight*)
wv=Table[{3+x0[[i]],3+Sqrt[mu]x0[[i]]+Sqrt[1-mu]y0[[i]]},{i,1,Length[x0]}];
w = {Ramp[wv[[All,1]]]};
v = {Ramp[wv[[All,2]]]};
Jout = (w\[Transpose] . w+v\[Transpose] . v)/size;
(*Calculate E neurons' response*)
iemisx=(w[[1]]-Mean[w[[1]]])*deltax[b,1,0]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,0];
iemisy=(w[[1]]-Mean[w[[1]]])*deltax[b,0,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,0,1];
iefull=(w[[1]]-Mean[w[[1]]])*deltax[b,1,1]+(v[[1]]-Mean[v[[1]]])*deltay[b,1,1];
rmisx=b*Ramp[iemisx-theta];
rmisy=b*Ramp[iemisy-theta];
rfull=b*Ramp[iefull-theta];
s=(rmisx-Mean[rmisx])*(rfull-Mean[rfull]);

(*Calculate stimulus input vector and put in the matrix*)
stimVec = Join[Table[w[[i]],{i,1,p}],Table[v[[i]],{i,1,p}]];
side=stimVec[[1;;2 p,2p+1;;]];
J0p=stimVec[[1;;2 p,1;;2 p]];
(*The order*)
chiral = Det[J0p];
idp=IdentityMatrix[2 p];
If[chiral<= 0,
idp=ReplacePart[idp,{1->idp[[2]],2->idp[[1]]}]
];

(*Diagonal elements*);
thre=0.8-Log[0.5]/4;
a0=HeavisideTheta[-s](1.4+12Exp[1.5s])+0.2HeavisideTheta[s](0.01+HeavisideTheta[s-thre]);
Jie = IdentityMatrix[size];
Jei = ConstantArray[0,{size,size}];
Jei[[All,1;;2p]]=stimVec\[Transpose];

kappa =lambda;
a =kappa*a0+(1-kappa)*e;
interp=kappa*J0p/size+(1-kappa)*idp;
Jie=ReplacePart[Jie,{i_,i_}:>a[[i]]];
Jie[[1;;2 p,1;;2 p]]=interp;
Jie[[1;;2 p,2p+1;;]] = kappa*side/size;
hImisx =Jie . rmisx;
hImisy =Jie . rmisy;
hIfull =Jie . rfull;
alambda = DiagonalMatrix[(1-kappa)/a];
Jei[[All,2p+1;;]]=Jout[[All,2p+1;;]] . alambda[[2p+1;;,2p+1;;]];
{{iemisx,iemisy,iefull,hImisx,hImisy,hIfull},Jie,Jei}
]
