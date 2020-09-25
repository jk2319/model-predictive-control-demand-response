function [problem,guess] = TermalOpt_mpc(x0,Tf,T_ref)
% Syntax:  [problem,guess] = TermalOpt_mpc
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: ICLOCS
% MAT-files required: case_study_sim.mat
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
% Lookup Table 1: U.S. 1976 Standard Atmosphere (altitude, density and pressure)
load case_study_sim

%%

% Plant model name, used for Adigator
InternalDynamics=@TermalOpt_Dynamics_Internal;
SimDynamics=@TermalOpt_Dynamics_Sim;

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=@gradCost_term_MPC;
problem.analyticDeriv.hessianLagrangian=@hessianLagrangian_term_MPC;
problem.analyticDeriv.jacConst=@jacConst_term_MPC;

% Settings file
problem.settings=@settings_TermalOpt;

%Initial Time. t0<tf
problem.time.t0_min=0;
problem.time.t0_max=0;
guess.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=Tf;     
problem.time.tf_max=Tf; 
guess.tf=Tf;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system
problem.states.x0=x0;

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=problem.states.x0; 
problem.states.x0u=problem.states.x0; 

% State bounds. xl=< x <=xu
problem.states.xl=[xlow];
problem.states.xu=[xup];

% State error bounds
problem.states.xErrorTol_local=[0.1 0.1];
problem.states.xErrorTol_integral=[0.1 0.1];

% State constraint error bounds
problem.states.xConstraintTol=[0.1 0.1];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[xlow]; 
problem.states.xfu=[xup];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[x0(1) x0(1)];
guess.states(:,2)=[x0(2) x0(2)];
%guess.states(:,3)=[x0(3) x0(3)];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[ulow];
problem.inputs.uu=[uup];

% Bounds on first control action
problem.inputs.u0l=[ulow];
problem.inputs.u0u=[uup]; 

% Input constraint error bounds
problem.inputs.uConstraintTol=[1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1];

% Guess the input sequences with [u0 uf]
uge= (ulow+uup)/2;
guess.inputs(:,1)=[uge(1) uge(1)];
guess.inputs(:,2)=[uge(2) uge(2)];
guess.inputs(:,3)=[uge(3) uge(3)];
guess.inputs(:,4)=[uge(4) uge(4)];
guess.inputs(:,5)=[uge(5) uge(5)];
guess.inputs(:,6)=[uge(6) uge(6)];
guess.inputs(:,7)=[uge(7) uge(7)];
guess.inputs(:,8)=[uge(8) uge(8)];
guess.inputs(:,9)=[uge(9) uge(9)];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.ng_eq=2;
problem.constraints.gTol_eq=[0.05 0.05];

problem.constraints.gl=[0 0 -inf 0 0 0 -inf 0];
problem.constraints.gu=[60 inf 0 inf inf inf 0 inf];
problem.constraints.gTol_neq=[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01];

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];
problem.constraints.bTol=[];

t_index=(0<=T_mex)&(T_mex<=Tf);
t_index_p=(T_ref<=T_mex)&(T_mex<=T_ref+Tf);

% store the necessary problem parameters used in the functions
problem.data.u_b0=u_b0;
problem.data.gamma1=gamma1(t_index_p);
problem.data.alpha1=alpha1(t_index_p);

problem.data.Price_carbon=Price_carbon;
problem.data.Carbon_em=Carbon_em(t_index_p);
problem.data.UV=UV;
problem.data.A=A;
problem.data.rhoAir=rhoAir;
problem.data.B_V=B_V;
problem.data.C_air=C_air ;
problem.data.n_ac=n_ac  ;
problem.data.C_Build=C_Build ;
problem.data.T_ex=T_ex(t_index_p) ;
problem.data.T_ex2=T_ex2(t_index_p) ;
problem.data.Ir=Ir(t_index_p);
problem.data.Pr_elec=Pr_elec(t_index_p);
problem.data.T_mex=T_mex(t_index);
problem.data.m_COP=m_COP;
problem.data.COP_init=COP_init;
problem.data.T_init=T_init;
problem.data.K_sat=K_sat;
%problem.data.Conv_S_S=Conv_S_S;
problem.data.eff_T_S=eff_T_S;
problem.data.Area_T_S=Area_T_S;
problem.data.TSo_Ap=TSo_Ap;
%problem.data.eff_SS=eff_SS;
%problem.data.xupSS=xupSS;
problem.data.m_cop_cool=m_cop_cool;
% Data on PV power output
problem.data.t1=t1;
problem.data.t2=t2;
problem.data.t3=t3;
problem.data.PV_s_Ap=PV_s_Ap;
% Total floor area
problem.data.AFl=AFl;
%Batteries
problem.data.dis_eff=dis_eff;
problem.data.ch_eff=ch_eff;
problem.data.uB_up=uB_up;
problem.data.ADR_max=ADR_max;
problem.data.DR_min=DR_min;

% Get function handles and return to Main.m
problem.data.InternalDynamics=InternalDynamics;
problem.data.functionfg=@fg;
problem.data.plantmodel = func2str(InternalDynamics);
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.sim.functions=SimDynamics;
problem.sim.inputX=[];
problem.sim.inputU=1:length(problem.inputs.ul);
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
problem.data.functions_unscaled=problem.functions_unscaled;
problem.data.ng_eq=problem.constraints.ng_eq;
problem.constraintErrorTol=[problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,vdat)

% L_unscaled - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    stageCost - Scalar or vectorized stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimization. 
%          Example: stageCost = 0*t;

%------------- BEGIN CODE --------------
%u1 = u(:,1); u2 = u(:,2);  u3 = u(:,3);
u_B_grid=u(:,5); 
u_s_grid=u(:,6);
s_t=u(:,8);
u_DR=u(:,9);
Price_carbon=vdat.Price_carbon;
Carbon_em=vdat.Carbon_em;
gamma1=vdat.gamma1; %DR Window 1000-1300 & 1800-2100
alpha1=vdat.alpha1; %Call 1100-1200, Call 1900-2000


c_em=interp1(vdat.T_mex,Carbon_em,t,'linear'); 
ct=interp1(vdat.T_mex,vdat.Pr_elec,t,'previous');  
g1=interp1(vdat.T_mex ,gamma1, t,'previous'); %,vdat.gamma1(end));
a1=interp1(vdat.T_mex ,alpha1, t,'previous'); %,vdat.alpha1(end));



 if size(t,1)>1
  delt=t(2:end)-t(1:end-1); 
  delt=[delt;delt(end)];
 else
   delt=vdat.T_mex(2)-vdat.T_mex(1); 
 end 
stageCost=ct(1:size(t,1)).*u_B_grid-0.98*ct(1:size(t,1)).*u_s_grid+Price_carbon*c_em(1:size(t,1)).*u_B_grid.*delt-a1(1:size(t,1)).*(0.2).*u_DR+(0.1).*s_t; 
 

%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 

% E_unscaled - Returns the boundary value cost
%
% Syntax:  boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

boundaryCost=0;




%------------- END OF CODE --------------


function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)

% b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    bc - column vector containing the evaluation of the boundary function 
%
%------------- BEGIN CODE --------------
varargin=varargin{1};
bc=[];  % [xf(1)-x0(1);xf(2)-x0(2);xf(3)-x0(3)];
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
%------------- BEGIN CODE --------------
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if ((strcmp(options.discretization,'hpLGR')) || (strcmp(options.discretization,'globalLGR')))  && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

