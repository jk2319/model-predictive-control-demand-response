    
function [df,dg,db]=jacConst_term_MPC(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data)

%  jacConst - Return the gradient of the plant model, path constraints and  boundary constraints 
%
% Syntax: [df,dg,db]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data)
%
% Inputs:
%    see directCollocation.m
%    Notice that the i-th state and input xi and ui, evaluated at the time instants t=[t0,...tk,...tf], 
%    are column  vectors taken as X(:,i) and U(:,i)
%
% Outputs:
%     df - Gradient of plant model with respect to  t, p, x, u
%     dg - Gradient of the path constraints with respect to t, p, x, u
%     db - Gradient of the boundary constraints with respect to  tf,p,x0,u0,uf,xf 
%    
%     Gradient df: 
%     df.flag - Set  df.flag=1 if the analytic form for df is supplied otherwise set df.flag=0; 
%     df.dx - Gradient of f(x,u,p,t)  wrt. x
%     df.du - Gradient of f(x,u,p,t)  wrt. u
%     df.dt - Gradient of f(x,u,p,t)  wrt. t
%     df.dp - Gradient of f(x,u,p,t)  wrt. p
%     If df.flag=1, the gradients have to be set considering the following rules:
%       1) If tf is a variable of the problem it is necessary to
%         specify  the derivative of f with respect to time t otherwise it possible to set df.dt=[].
%       2) If p=[p1, p2, ....,pn] is a variable of the problem it is necessary to
%         specify  the derivatives of f with respect to parameters pi otherwise it possible to set df.dp=[].  
%       3) df.dx and df.du have to be specified 
%       4) df.dx, df.du,  df.dp and df.dt are cell arrays and  must be vectorized.
%          For instance df.dx{i} contains the derivative  of f(x,u,p,t) with respect to the i-th
%          state variable xi. df.dx{i} is a matrix with n column where n is the dimension of the system.        
%          Each row of df.dx{i} stores the derivative of f(x,u,p,t)  (composed by  f1, f2, ....fn) 
%          evaluated at a time instant. The derivative of fj with respect to the i-th state variable, 
%          evaluated along all the horizon, corresponds to the i-th column of df.dx{i}(:,j);
%          Here, if no you find df.dt difficult to derive, can leave df.dt=[];
%          
%     Gradient dg: 
%     dg.flag - Set  dg.flag=1 if the analytic form for dg is supplied otherwise set dg.flag=0; 
%     dg.dx - Gradient of g(x,u,p,t)  wrt. x
%     dg.du - Gradient of g(x,u,p,t)  wrt. u
%     dg.dt - Gradient of g(x,u,p,t)  wrt. t
%     dg.dp - Gradient of g(x,u,p,t)  wrt. p
%     If dg.flag=1, the gradients have to be set following the same rules explained for df.
%     Here the cell arrays  df.dx{i} is a matrix with ng column where ng is the number of path 
%     constraints. 
%
%     Gradient db;
%     db.flag - Set  db.flag=1 if the analytic form for db is supplied otherwise set db.flag=0; 
%     db.dtf  - Gradient of b(x0,xf,u0,uf,p,tf)  wrt. tf
%     db.dp   - Gradient of b(x0,xf,u0,uf,p,tf)  wrt. p
%     db.dx0  - Gradient of b(x0,xf,u0,uf,p,tf)  wrt. x0
%     db.du0  - Gradient of b(x0,xf,u0,uf,p,tf)  wrt. u0
%     db.dxf  - Gradient of b(x0,xf,u0,uf,p,tf)  wrt. xf
%     db.duf  - Gradient of b(x0,xf,u0,uf,p,tf)  wrt. uf
%     If the derivative of b does not depend on some variable it is possible to set
%     the respective variable as the empty matrix. For instance if it does not depend on tf set
%     db.dtf=[]; The gradient of b with respect to a vector variable  has to be specified in a matrix
%     with dimension nb x s where s is the size of the vector variable. For instance db.dxf has to be a matrix 
%     with dimension nb x n and the entry (i,j) correspond to the derivative of th i-th constraints 
%     with respect of the j-th state variable in xf.
%     
% 
%
% Other m-files required: none
%         Subfunctions: none
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

%[ X,~,U,~,P,x0,xf,u0,uf,p,data] = batchScaleBack(X,[],U,[],P,x0,xf,u0,uf,p,data);
[ X,Xr,U,Ur,P,x0,xf,u0,uf,p,data ] = batchScaleBack(X,[],U,[],P,x0,xf,u0,uf,p,data);   % It scales back system variables

%[nt,np,n,m]=deal(data.sizes{1:4});   % extract variable dimension
% if nt is 1 the vector df.dt can not be empty

% x1 - dwelling Temperature
% x2 - state of charge of seasonal storage
% X3 - state of charge of LiFePO4 battery 
% x4 - battery degradation

% u1 - Heating  HP power consumption 
% u2 - Seasonal storage power electricity consumption
% u3 -Cooling HP power consumpyion 
% ud - Battery discharge rate
% uc - Battery charge rate
% u_B_grid - Power bought from the grid
% u_s_grid - Power inject into the grid

UV= data.data.UV;
A= data.data.A;
rhoAir= data.data.rhoAir;
B_V= data.data.B_V;
C_air= data.data.C_air;
n_ac=data.data.n_ac;
C_Build = data.data.C_Build;
m_COP=data.data.m_COP;
COP_init=data.data.COP_init;
T_init=data.data.T_init;
%Conv_S_S=data.data.Conv_S_S;
K_sat=data.data.K_sat;
m_cop_cool=data.data.m_cop_cool;
T_ex=data.data.T_ex ;
T_mex=data.data.T_mex;
alpha1=data.data.alpha1;


a1=interp1(T_mex ,alpha1,min(max(t,0),T_mex(end)),'previous',T_ex(end)); %,vdat.alpha1(end));
Ta=interp1(T_mex,T_ex,min(max(t,0),T_mex(end)),'linear',T_ex(end));
% x1 room temperature
COPTA=m_COP*(Ta-T_init)+COP_init;  % COP for edHP
% Thermal balance
%dx(:,1)=(UV*A+rhoAir* B_V*C_air*n_ac)*(Ta-x1)./C_Build+COPTA.*u1./C_Build-m_cop_cool*u2./C_Build;

% battery storage
  %dx(:,2)=vdat.ch_eff*uc-ud./vdat.dis_eff;


Lt=ones(size(t));

df.flag=1;               % df.flag=1 if the analytic is supplied otherwise set
                         % df.flag=0;
                         
df.dp{1}=[];              % Derivatives of f(x,u,p,t) wrt. p
df.dt{1}=[];              % Derivatives of f(x,u,p,t) wrt. t, if not sure, can set df.dt=[];
df.dx{1}=[-(UV*A+rhoAir* B_V*C_air*n_ac)*Lt./C_Build, 0*t]; %0*t % Derivatives of f(x,u,p,t) wrt. x1
df.dx{2}=[0*t, 0*t]; % Derivative of f(x,u,p,t) wrt. x2
%df.dx{3}=[0*t, 0*t, 0*t]; % Derivative of f(x,u,p,t) wrt. x3 

%dx1(x,u,p,t) dx2(x,u,p,t) dx3(x,u,p,t)
df.du{1}=[Lt.*COPTA./C_Build, 0*t]; %, 0*t];  % Derivative of f(x,u,p,t) wrt. u1
%df.du{2}=[Lt*Conv_S_S./C_Build, -Conv_S_S*Lt,    0*t]; 
df.du{2}=[-m_cop_cool*Lt./C_Build, 0*t];  % Derivative of f(x,u,p,t) wrt. u2
df.du{3}=[0*t, -Lt./data.data.dis_eff];  % Derivative of f(x,u,p,t) wrt. u3 (ud)
df.du{4}=[0*t,  data.data.ch_eff*Lt];  % Derivative of f(x,u,p,t) wrt. u4 (uc)
df.du{5}=[0*t,  0*t];  % Derivative of f(x,u,p,t) wrt. u5 (u_B_grid)
df.du{6}=[0*t,  0*t];  % Derivative of f(x,u,p,t) wrt. u6 (u_S_grid)
df.du{7}=[0*t,  0*t];  % Derivative of f(x,u,p,t) wrt. u7 (u_ADR)
df.du{8}=[0*t,  0*t];  % Derivative of f(x,u,p,t) wrt. u8 (s_t)
df.du{9}=[0*t,  0*t];  % Derivative of f(x,u,p,t) wrt. u9 (u_DR)

dg.flag =1;
%     Gradient dg: 
%     dg.flag - Set  dg.flag=1 if the analytic form for dg is supplied otherwise set dg.flag=0; 
%     dg.dx - Gradient of g(x,u,p,t)  wrt. x
%     dg.du - Gradient of g(x,u,p,t)  wrt. u
%     dg.dt - Gradient of g(x,u,p,t)  wrt. t
%     dg.dp - Gradient of g(x,u,p,t)  wrt. p

%g_eq1=u_B_grid-u_s_grid+ud-uc+PV_el-u1-u2-u3-u_ADR;
%g_eq2=u_DR-(a1).*(u_b0-u_B_grid)-(a1).*u_ADR;

% neq1=x2-(1-g1).*10-g1.*30+g1.*(a1).*(30-10);
% 
% neq2=x1-(1-g1).*18-g1.*19+g1.*(a1).*(19-18);
% 
% neq3=x1-(1-g1).*22-g1*21-g1.*(a1).*(22-21);
% 
% neq4=u_DR+(a1).*s_t-(a1).*1; %u_min= 1kW
% 
% neq5=s_t;
% 
% neq6=u_ADR;
% 
% neq7=u_ADR-(a1)*5; %DR_max= 5kW
% 
% neq8=(a1).*u_b0+(1-(a1)).*80-u_B_grid; %u_B_grid uup=80kW

dg.dp{1}=[];
dg.dt{1}=[];
dg.dx{1}=[0*Lt, 0*Lt, 0*Lt, Lt, Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];
dg.dx{2}=[0*Lt, 0*Lt, Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];
%dg.dx{3}=[0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];

dg.du{1}=[-Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];
dg.du{2}=[-Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];
dg.du{3}=[Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];
dg.du{4}=[-Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt];
dg.du{5}=[Lt, (a1).*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, -Lt]; %wrt u_B_grid
dg.du{6}=[-Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt]; %wrt u_S_grid
dg.du{7}=[-Lt, -(a1).*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, Lt, Lt, 0*Lt]; %wrt u_ADR
dg.du{8}=[0*Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt, (a1).*Lt, Lt, 0*Lt, 0*Lt, 0*Lt]; %wrt s_t
dg.du{9}=[0*Lt, Lt, 0*Lt, 0*Lt, 0*Lt, Lt, 0*Lt, 0*Lt, 0*Lt, 0*Lt]; %wrt u_DR

db.flag=1;
db.dtf=[];
db.dt0=[];
db.dp=[];
db.dx0=[];
db.du0=[];
db.duf=[];
db.dxf=[];  %  example [1./sqrt(xf(1))./xf(1)/2, 0, 1; 0, 1, 0];     % derivatives with respect to [x1(tf) ,x2(tf), x3(tf)]


[ df,dg,db ] = batchScalejacConst(df,dg,db,data);
%------------- END CODE --------------
