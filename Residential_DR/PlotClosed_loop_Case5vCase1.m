

%% figure
% xx=linspace(solution.T(1,1),solution.tf,10000);

% loading system data
load case_study_sim
step=4;

par=problem.data;
t_index=T_mex<=T(end);
t_index_p=(T_ref<=T_mex)&(T_mex<=T_ref+T(end));  
Electr_pr=Pr_elec(t_index_p);
T_env=T_ex(t_index_p);
Rad_sol=Ir(t_index_p);

% Note all the data are shifted at the beginning of the time  

xx=linspace(T(1),T(end),10000);
ct=interp1(T_mex(t_index),Electr_pr,xx,'previous');
T_outside=interp1(T_mex(t_index),T_env,xx,'linear');
I_rad=interp1(T_mex(t_index),Rad_sol,xx,'linear');

%% States = Temperature & SoC Battery
% Internal temperature x1 (/w bounds) vs outside temperature
startt=0;
endd=T(end)+step;
figure()
sgtitle('Case 5-A vs Case 1-A') %CASE!!!!!!!!!!
subplot(2,1,1)
hold on
plot(T,Xbase(:,1),':','LineWidth',3)
plot(T,X(:,1),'-','LineWidth',1.5)
plot(xx,T_outside,'k-','LineWidth',2 )
plot([T(1),simtime+T(end)],problem.states.xu(1)*ones(1,2),'r--','LineWidth',1)
plot([T(1),simtime+T(end)],problem.states.xl(1)*ones(1,2),'r--','LineWidth',1)
xlim([startt endd])
xlabel('Time [h]')
ylabel('(x_{1,t}), [\circ C]')
title('Dwelling Temperature')
legend('Case 5-A','Case 1-A','External Temperature') %CASE!!!!!!!!!!!!!
grid minor; box on;
ax = gca;
% Set x and y font sizes
% ax.XAxis.FontSize = 15;
% ax.YAxis.FontSize = 24;
% The below would set everything: title, x axis, y axis, and tick mark label font sizes.
%ax.FontSize = 16;
% Bold all labels.
%ax.FontWeight = 'bold';

% SoC of Battery with charge and discharge statuses
subplot(2,1,2)
hold on
plot(T,Xbase(:,2),':','LineWidth',2)
plot(T,X(:,2),'-','LineWidth',2)
%stairs(T,par.ch_eff*U(:,5)-U(:,4)./par.dis_eff,'k','LineWidth',2);
stairs(T,U(:,3),'Color',[0,0.7,0.9],'LineWidth',2 )
stairs(T,U(:,4),'r','LineWidth',2 )
plot([T(1),T(end)],problem.states.xu(2)*ones(1,2),'b--','LineWidth',1)
plot([T(1),T(end)],problem.states.xl(2)*ones(1,2),'b--','LineWidth',1)
% plot([T(1),T(end)],problem.inputs.ul(4)*ones(1,2),'k--')
plot([T(1),T(end)],-problem.inputs.uu(3)*ones(1,2),'k--','LineWidth',1)
%plot([T(1),T(end)],problem.inputs.ul(5)*ones(1,2),'k--')
plot([T(1),T(end)],problem.inputs.uu(4)*ones(1,2),'k--','LineWidth',1)
xlim([startt endd])
xlabel('Time [h]')
ylabel('Battery SOC (x_{2,t})')
title('Energy of LiFePO4 Battery')
legend( 'Case 5-A (kWh)', 'Case 1-A (kWh)', 'Battery-Discharge (kW)',  'Battery-Charge (kW)')
% CASE!!!!!!!!!!
grid minor; box on;

P_TS=problem.data.eff_T_S*I_rad*problem.data.TSo_Ap*problem.data.AFl;   % Power from thermal solar
P_PV=(problem.data.t1*I_rad+ problem.data.t2*I_rad.^2+problem.data.t3*I_rad.*T_outside)*problem.data.PV_s_Ap*problem.data.AFl;   % Power produced by PVs


%% Price of electricity, Load vs Generation, Bought vs Sold vs DR, Penalty
%Load = HP heating/cooling + battery charge rate
LoadBase=Ubase(:,1)+Ubase(:,2); %+U(:,4); %+U(:,2); SS
Load=U(:,1)+U(:,2); %+U(:,4); %+U(:,2); SS
generationBase=P_PV+interp1(T,Ubase(:,3),xx,'linear',U(end,3));
generation=P_PV+interp1(T,U(:,3),xx,'linear',U(end,3));
%speval(solution,'U',4,xx);
sold=U(:,6);
Bought=U(:,5);
DR=U(:,9);

start1=0;
end1=T(end)+step;
figure() 
sgtitle('Case 5-A vs Case 1-A')
% % Price of electricity
subplot(4,1,1)
plot(xx,ct,'LineWidth',2)
%xlim([105 130]) 
xlim([start1 end1])
xlabel('Time [h]')
ylabel('Price (£/kWh)')
title('Electricity Tariff')
%legend( '£/h')
grid minor; box on;

subplot(4,1,2)
hold on
plot(T,LoadBase,'k:','LineWidth',1.5)
plot(T,Load,'r','LineWidth',1.5)
%xlim([105 130]) 
xlim([start1 end1])
xlabel('Time [h]')
ylabel('Power (kW)')
title('Load Demand, HP')
legend( 'Case 5-A', 'Case 1-A') %CASE!!!!!!!!!!
grid minor; box on;

subplot(4,1,3)
hold on
plot(xx,generationBase,'k:','LineWidth',1.5)
plot(xx,generation,'b--','LineWidth',1.5)
%xlim([105 130]) %T(end)+step
xlim([start1 end1])
xlabel('Time [h]')
ylabel('Power (kW)')
title('Generation, PV + u_d')
legend( 'Case 5-A', 'Case 1-A') %CASE!!!!!!!!!!!
grid minor; box on;

subplot(4,1,4)
hold on
plot(T,Ubase(:,5),'k:','LineWidth',2)
plot(T,Bought,'LineWidth',1.5) %,'m--',
%plot(T,sold,'LineWidth',1.5 ) %,'k--'
plot(T,DR,'LineWidth',1.5) %,'g--',
%xlim([105 130]) %T(end)+step
xlim([start1 end1])
xlabel('Time [h]')
ylabel('Power (kW)')
title('Power traded with grid')
legend('Bought power (Base)','Bought power','DR')%,'Sold power')
grid minor; box on;

%PENALTY
% subplot(4,1,4)
% hold on
% plot(T,U(:,8),'LineWidth',1.5) %,'m--',
% xlim([start1 end1])
% xlabel('Time [h]')
% ylabel('DR Penalty (kW)')
% legend('Penalized Energy (kW)')
% grid minor; box on;
% ylim([0 5]); 



