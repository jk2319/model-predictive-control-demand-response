function [dx,g_eq] = TermalOpt_Dynamics_Internal(x,u,p,t,vdat)
%Supersonic Aircraft Minimum Fuel Climb Problem - Dynamics - Internal
%
% Syntax:  
%          [dx] = Dynamics(x,u,p,t,vdat)	(Dynamics Only)
%          [dx,g_eq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Eqaulity Path Constraints)
%          [dx,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics and Inqaulity Path Constraints)
%          [dx,g_eq,g_neq] = Dynamics(x,u,p,t,vdat)   (Dynamics, Equality and Ineqaulity Path Constraints)
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
%    g_eq - constraint function for equality constraints
%    g_neq - constraint function for inequality constraints
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

x1 = x(:,1); x2 = x(:,2);
u1 = u(:,1); u2 = u(:,2);  ud=u(:,3);   uc=u(:,4);  
u_B_grid=u(:,5); u_s_grid=u(:,6);

%u_b0=u(:-1,5);
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

UV= vdat.UV;
A= vdat.A;
rhoAir= vdat.rhoAir;
B_V= vdat.B_V;
C_air= vdat.C_air;
n_ac=vdat.n_ac;
C_Build = vdat.C_Build;
T_ex=  vdat.T_ex ;
T_mex=vdat.T_mex;
m_COP=vdat.m_COP;
COP_init=vdat.COP_init;
T_init=vdat.T_init;
Ir=vdat.Ir;
AT_S=vdat.TSo_Ap*vdat.AFl;
PV_S=vdat.PV_s_Ap*vdat.AFl;
%Conv_S_S=vdat.Conv_S_S;
eff_T_S=vdat.eff_T_S;
Area_T_S=vdat.Area_T_S;
%eff_SS=vdat.eff_SS;
K_sat=vdat.K_sat;
%xupSS=vdat.xupSS;
m_cop_cool=vdat.m_cop_cool;



% Solar Radiance

IRd=interp1(T_mex,Ir,t,'linear'); 
% Ambient Temperature
Ta=interp1(T_mex,T_ex,t,'linear'); 

%Thermal solar term
P_T_S=eff_T_S*IRd.*AT_S;

% Power produced by PVs
PV_el=(vdat.t1*IRd+vdat.t2*IRd.^2+vdat.t3*IRd.*Ta).*PV_S;

% x1 room temperature

COPTA=m_COP*(Ta-T_init)+COP_init;  % COP for edHP
  
% Thermal balance
dx(:,1)=(UV*A+rhoAir* B_V*C_air*n_ac)*(Ta-x1)./C_Build+COPTA.*u1./C_Build-m_cop_cool*u2./C_Build;


% battery storage
dx(:,2)=vdat.ch_eff*uc-ud./vdat.dis_eff;


% LiFePO4 Battery dis_eff ch_eff

g_eq=u_B_grid-u_s_grid+ud-uc+PV_el-u1-u2;

%g_neq=u_s_grid;



