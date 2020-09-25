function [dL,dE]=gradCost_myProblem(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,t0,tf,data)

%GRADCOST - Return the gradient of the cost in analytic form
%
% Syntax:  [dL,dE]=gradCost(L,X,Xr,U,Ur,P,t,E,x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%   dL - Gradient of the stage cost L  wrt. t,p,x,u 
%   dE - Gradient of E with respect to tf,p,x0,u0,uf,xf
%  
%  Gradient dL:
%  If the analytic form for L is supplied, set  dL.flag=1 otherwise set
%  dL.flag=0. If tf is a variable of the problem  and dL.flag=1 it is necessary to
%  specify  the derivative of L with respect to time t. 
%  You need to specify dL.dx and dL.du whenever dL.flag=1. 
%  dL.dx, dL.du,  dL.dp and dL.dt must be vectorized.
%  The derivative of L with respect to the i-th state variable, evaluated along all the horizon,
%  corresponds to the i-th column of dL.dx. The same rule holds for dL.du, dL.dp, dL.dt 
%  The i-th state and input
%  xi and ui, evaluated at the time instants t=[t0,...tk,...tf], are column vectors 
%  taken as X(:,i) and U(:,i)
%  The i-th state and input references  xri and uri for the i-th variable are column vectors taken 
%   as Xr(:,i) and Ur(:,i)  
%    
%  Gradient dE:
%  If the analytic form for E is supplied set  dE.flag=1 otherwise set
%  dE.flag=0. If the derivative of E does not depend on some variable it is possible to set
%  the respective variable as the empty matrix. For instance if it does not depend on x0 set
%  dE.dx0=[];
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
[ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScaleBack(X,Xr,U,Ur,P,x0,xf,u0,uf,p,data);   % It scales back system variables

% Stage cost gradients
[n,m]=deal(data.sizes{3:4});
Lt=ones(size(t));


Price_carbon=data.data.Price_carbon;
Carbon_em=data.data.Carbon_em;

c_em=interp1(data.data.T_mex,Carbon_em,t,'linear');
ct=interp1(data.data.T_mex,data.data.Pr_elec,t,'previous');  


 if size(t,1)>1
  delt=t(2:end)-t(1:end-1); 
  delt=[delt;delt(end)];
 else
  delt=data.data.T_mex(2)-data.data.T_mex(1); 
 end 
 
%Computing derivative with resect to u  
% Derivative with respect to U(:,6) (Bouught power)
e5=zeros(1,m);
e5(5)=1;
du5=kron(ct(1:size(t,1))+Price_carbon*c_em(1:size(t,1)).*delt,e5);
%du6=kron(ct,e6);

% Derivative with respect to U(:,7) (Sold power)
e6=zeros(1,m);
e6(6)=1;
du6=kron(-0.98*ct(1:size(t,1)),e6);


dL.flag=1; 
dL.dp=[];
dL.dt=[];
dL.dx=kron(zeros(1,n),Lt);
dL.du=du5+du6;



% Terminal cost gradients
dE.flag=1;
dE.dt0=[];
dE.dtf=[];
dE.dp=[];
dE.dx0=[];
dE.du0=[];
dE.dxf=[];
dE.duf=[];

[ dL,dE ] = batchScaleGradCost(dL,dE,data);
%------------- END CODE --------------



