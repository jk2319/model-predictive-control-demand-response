function dx = TermalOpt_Dynamics_Sim(x,u,p,t,vdat)
% Supersonic Aircraft Minimum Fuel Climb Problem - Dynamics - simulation
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)
% 
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
%
%------------- BEGIN CODE --------------

x1 = x(:,1); u1 = u(:,1); u2 = u(:,2); 
ud=u(:,3);   uc=u(:,4);  s_t=u(:,8);


UV= vdat.UV;
A= vdat.A;
rhoAir= vdat.rhoAir;
B_V= vdat.B_V;
C_air= vdat.C_air;
n_ac=vdat.n_ac;
C_Build = vdat.C_Build;
T_ex=  vdat.T_ex ;
T_ex2=vdat.T_ex2;
T_mex=vdat.T_mex;
m_COP=vdat.m_COP;
COP_init=vdat.COP_init;
T_init=vdat.T_init;
Ir=vdat.Ir;
%T_hour=vdat.T_hour;

%Conv_S_S=vdat.Conv_S_S;
eff_T_S=vdat.eff_T_S;
Area_T_S=vdat.Area_T_S;
%eff_SS=vdat.eff_SS;
K_sat=vdat.K_sat;
%xupSS=vdat.xupSS;
m_cop_cool=vdat.m_cop_cool;
AT_S=vdat.TSo_Ap*vdat.AFl;


% Solar Radiance
IRd=interp1(T_mex,Ir,t,'linear'); 
% Ambient Temperature
Ta=interp1(T_mex,T_ex,t,'linear');
Ta2=interp1(T_mex,T_ex2,t,'linear'); %inexact

%Thermal solar term
P_T_S=eff_T_S*IRd.*AT_S;

Ta=interp1(T_mex,T_ex,t,'linear',T_ex(end));

Ta2=interp1(T_mex,T_ex2,t,'linear',T_ex(end)); %inexact
% x1 room temperature

COPTA=m_COP*(Ta-T_init)+COP_init;  % COP for edHP
  
COPTA=m_COP*(Ta2-T_init)+COP_init;  % COP for edHP inexact

% Thermal balance
dx(:,1)=(UV*A+rhoAir* B_V*C_air*n_ac)*(Ta2-x1)./C_Build+COPTA.*u1./C_Build-m_cop_cool*u2./C_Build;

% LiFePO4 Battery dis_eff ch_eff
dx(:,2)=vdat.ch_eff*uc-ud./vdat.dis_eff;
% 
% % Slack variable for penalised energy
% dx(:,3)=s_t;
