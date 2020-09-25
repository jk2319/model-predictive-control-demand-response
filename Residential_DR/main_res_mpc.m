% Main script to  define and solve solve the Optimal Control Problem using
% a Receding Horizon formulation (MPC)
%

% Define the problem

clear all; %close all;  
format compact;

load External_data   % read forecast of external data
 
addpath('MPC_vanilla')

% The used data have a temporal resolution of 15 minutes and consequnetly some of the parameters can be defined as a multiple of a quartar 
% of an hour only at the moment. 

tstep=1/4;          % Time step of simulation (a quarter of an hour). At the moment can not be decreased.
t_end=168;      % End of the simulation time in hours
T_ref=0;  %5040;    % At the moment it can be set a multiple of 0.25 only.    
Tf=24*2;            % prediction horizons in hours (open loop run)
N_N_node=24;        % It must contain a  positive value such that N_N_node/4 is an integer. 
                    % It is used to define the number of discretization
                    % points in each interval of a 1/4 of hour. If
                    % N_N_node=4 the discretization points are only
                    % extremal to the interval (i.e in 0 and 1/4) if
                    % N_N_node is 8, the number of discretization points in
                    % [0,1/4] is 3 (see N_input) and they are  0, 1/8, 1/4
x0=[20,0];          % Initial condition of the state.  


N_node=Tf*4+1;          % Number of nodes of the external signals for the chosen prediction horizon Tf. 
                        % This is because the data have temporal resolution of 15 minutes  
N_input= N_N_node/4+1;  % number of discretization point in a single interval of a 1/4 of an hour
Ist=find(T_mex>=T_ref);
T_start_data=Ist(1)-1;


options= settings_TermalOpt(Tf*N_N_node+1);          % Get options and solver settings
[problem,guess]=TermalOpt_mpc(x0,Tf,T_ref);          % Fetch the problem definition


par=problem.data;
[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});     % Retrieve parameters describing the size of the problem

[L]=deal(data.functions{1});                                                                                                                                                                            
Ldata=data.data;
Ldata=rmfield(Ldata,'Xscale');  

% Define simulation parameters
odesolver='ode113';
problem.sim.functions=@TermalOpt_Dynamics_Sim_exact;   
% problem.sim.functions=@TermalOpt_Dynamics_Sim_Inexact;

% Initialise vectors
%update_steps=[0:tstep:t_end]';
%n_steps=length(update_steps);
T=[];               % Initialize vector of times for the closed loop formulation
X=[];               % Initialize vector of states for the closed loop formulation
U=[];               % Initialize vector of inputs for the closed loop formulation
status_ref=[];      % Initialize vector storing the status of the optimization
cost_cl=[];         % Initialize vector storing the cost of the optimization
cost_cl2=[];        % Initialize vector storing the cost of the modified simulation
simtime=0;          % initialize simulation time (Since the time is not a decision variable and the prediction horizon is constant simstart starts 
iter_step=1;        % Intialize iteration step
tic;
for simtime=0:tstep:t_end
       
  % Solve optimisation problem
  [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
  
   %simulate closed loop system 
   [ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.25 );
   
   if (status.status>1)||(status.status<0)  
      break  % if the status of the optimisation is not acceptable interrupt the loop 
   end  

   % Update closed loop
   U=[U;mean(solution.U(1:N_input-1,:))];             % Store input solution at the current simulation time; in case you simulate the closed loop use uv
   X=[X;mean(solution.X(1:N_input-1,:))];             % Store state solution at the current simulation time; in case you simulate the closed loop use xv          
   T=[T;simtime+mean(solution.T(1:N_input-1))];       % Store  current evaluation times; in case you simulate the closed loop use tv
   status_ref=[status_ref;status.status];       % Store  current status of the optimization problem
   
   cost_cl=[cost_cl;data.map.w(1:N_input-1)'*((solution.tf-solution.t0)*L(solution.X(1:N_input-1,:),[],solution.U(1:N_input-1,:),[],repmat(solution.p,N_input-1,1),solution.T(1:N_input-1,:),Ldata))];

   cost_cl2=[cost_cl2;data.map.w(1:N_input-1)'*((solution.tf-solution.t0)*L(xv(1:N_input-1,:),[],uv(1:N_input-1,:),[],repmat(solution.p,N_input-1,1),tv(1:N_input-1,:),Ldata))];
%       snm=ones(M,1); 
%    sol{phaseNo}.cost= mp.w'*((tf- t0)*L(X,Xr,U,Ur,P,t,vdat).*snm)+E(x0,xf,u0,uf,p,t0,tf,vdat);
%    
 if simtime<=t_end-tstep %simtime is current time solving the opti problem
     
         
%      if(simtime inside, (gamma=1))
%          {
%              alpha3=rand<0.2 ; only if dont have a call (call=0), flag variable
%              (simtime:simtime+3*tstep)
%              for 1 hr
%               }
%  else {0}
% call must be 1 hr and cannot receive another call (call=1) when they already have
% a call.
  % Update parameters of the optimisation problem; If you simulate the closed loop use the corresponding value in xv in place of solution.X(N_input,:)   
         if options.scaling
            %x0=scale_variables(solution.X(N_input,:),data.data.Xscale,data.data.Xshift); 
            x0=scale_variables(xv,data.data.Xscale,data.data.Xshift); 
        else
            %x0=solution.X(N_input,:); 
            x0=xv;
         end     
         
        % Update the initial condition and guesses
        infoNLP.z0(:,1)=[solution.z(1:nt+np);solution.z(nt+np+m+n+1:end);solution.z(end-n-m+1+nt+np:end)];    % Update initial guess 
        infoNLP.z0(nt+np+1:nt+np+n,1)=x0';   % Update the initial condition from new system information
        infoNLP.zl(nt+np+1:nt+np+n,1)=x0';
        infoNLP.zu(nt+np+1:nt+np+n,1)=x0';
        data.x0t(:,1)=x0';
        data.data.u_b0=solution.U(1,5); %use previous value of u_grid, always use 1st value of previous opti problem


  % update forecast of exogenous parameters
  
        data.data.Carbon_em=Carbon_em(T_start_data+iter_step+1:T_start_data+iter_step+N_node);    %we assume that forecast data are available at each optimisation step
        data.data.Pr_elec=Pr_elec(T_start_data+iter_step+1:T_start_data+iter_step+N_node);        %we assume that forecast data are available at each optimisation step
        data.data.T_ex=T_ex(T_start_data+iter_step+1:T_start_data+iter_step+N_node); 
        data.data.T_ex2=T_ex2(T_start_data+iter_step+1:T_start_data+iter_step+N_node); %inexact
        data.data.Ir=Ir(T_start_data+iter_step+1:T_start_data+iter_step+N_node); 
        data.data.gamma1=gamma1(T_start_data+iter_step+1:T_start_data+iter_step+N_node);    %we assume that forecast data are available at each optimisation step
        data.data.alpha1=alpha1(T_start_data+iter_step+1:T_start_data+iter_step+N_node);    %we assume that forecast data are available at each optimisation step
        Ldata.Pr_elec=data.data.Pr_elec;        
        Ldata.Carbon_em=data.data.Carbon_em;  
        Ldata.alpha1=data.data.alpha1;
        
        
        iter_step=iter_step+1; 
 end     
 fprintf('Time: %.2f\n', simtime)
end
toc
fprintf('Ended at: %.2f\n', simtime)