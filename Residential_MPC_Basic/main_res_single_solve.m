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
t_end=1/5;          % End of the simulation time in hours
T_ref=0;            % At the moment it can be set a multiple of 0.25 only.    
Tf=24*2;            % prediction horizons in hours
N_N_node=24;        % It must contain a  positive value such that N_N_node/4 is an integer. 
                    % It is used to define the number of discretization
                    % points in each interval of a 1/4 of hour. If
                    % N_N_node=4 the discretization points are only
                    % extremal to the interval (i.e in 0 and 1/4) if
                    % N_N_node is 8, the number of discretization points in
                    % [0,1/4] is 3 (see N_input) and they are  0, 1/8, 1/4
x0=[18,0];        % Initial condition of the state.  


N_node=Tf*4+1;          % Number of nodes of the external signals for the chosen prediction horizon Tf. 
                        % This is because the data have temporal resolution of 15 minutes  
N_input= N_N_node/4+1;  % number of discretization point in a single interval of a 1/4 of an hour
Ist=find(T_mex>=T_ref);
T_start_data=Ist(1)-1;


options= settings_TermalOpt(Tf*N_N_node+1);          % Get options and solver settings
[problem,guess]=TermalOpt_mpc(x0,Tf,T_ref);          % Fetch the problem definition


par=problem.data;
[solution,MRHistory]=solveMyProblem( problem,guess,options);


%%


figure
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(4), problem.inputs.ul(4)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(4), problem.inputs.uu(4)],'r-' )
plot(solution.T(:,1),speval(solution,'U',4,solution.T),'b-')
xlim([0 solution.tf])
xlabel('Time [hrs]')
ylabel('Control Input 4 (at collocation point)')
grid on

tt=linspace(solution.T(1,1),solution.tf,100000);

figure
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(4), problem.inputs.ul(4)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(4), problem.inputs.uu(4)],'r-' )
plot(tt,speval(solution,'U',4,tt),'b-')
xlim([0 solution.tf])
xlabel('Time [hrs]')
ylabel('Control Input 4 (cont. trajectory)')
grid on
