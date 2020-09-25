

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
caseX=2-A;
figure()
subplot(2,1,1)
hold on
%plot(T,Xbase(:,1),':','LineWidth',3)
plot(T,X(:,1),'-','LineWidth',1.5)
plot(xx,T_outside,'k-','LineWidth',2 )
plot([T(1),simtime+T(end)],problem.states.xu(1)*ones(1,2),'r--','LineWidth',1)
plot([T(1),simtime+T(end)],problem.states.xl(1)*ones(1,2),'r--','LineWidth',1)
%plot([T(1),T(end)],problem.states.xu(1)*ones(1,2),'r--','LineWidth',2)
%plot([T(1),T(end)],problem.states.xl(1)*ones(1,2),'r--','LineWidth',2)
% plot(T(end)+solution.T,solution.X(:,1),'b-' )
% plot(xx,T_outside,'k-' )
% plot([solution.T(1),T(end)+solution.T(end)],problem.states.xu(1)*ones(1,2),'r--')
% plot([solution.T(1),T(end)+solution.T(end)],problem.states.xl(1)*ones(1,2),'r--')
% plot(tv,xv(:,1),'k-.')
% plot(tv,xv(:,2),'k-.')
%xlim([0 solution.tf])
xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('(x_{1,t}), [\circ C]')
title('Dwelling Temperature')
legend('Base','Case 2-A','External Temperature')
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
sgtitle('Base Case - Winter')
subplot(2,1,2)
hold on
%plot(T,Xbase(:,2),':','LineWidth',2)
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
xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('Battery SOC (x_{2,t})')
title('Energy of LiFePO4 Battery')
%legend( 'Base (kWh)', 'Case 1-B (kWh)', 'Battery-Discharge (kW)',  'Battery-Charge (kW)')
legend( 'Energy-Battery (Kwh)',  'Power-Battery-discharge (Kw)','Power-Battery-charge (Kw)')
grid minor; box on;

P_TS=problem.data.eff_T_S*I_rad*problem.data.TSo_Ap*problem.data.AFl;   % Power from thermal solar
P_PV=(problem.data.t1*I_rad+ problem.data.t2*I_rad.^2+problem.data.t3*I_rad.*T_outside)*problem.data.PV_s_Ap*problem.data.AFl;   % Power produced by PVs


%% Heat generated from solar panel (not energy supplied by PV)
% figure()
% hold on
% % %plot(T,X(:,2),'b-','LineWidth',2 )
% %plot([T(1),T(end)],problem.states.xu(2)*ones(1,2),'r--','LineWidth',2)
% %plot([T(1),T(end)],problem.states.xl(2)*ones(1,2),'r--','LineWidth',2)
% plot(xx,P_TS,'k--','LineWidth',1)
% xlim([0 T(end)+step])
% xlabel('Time [h]')
% ylabel('Heat generation from thermal solar panel (kW)')
% legend('Power of Thermal Solar')
% grid minor
% 

% figure(3)
% hold on
% plot(xx,P_TS,'r-')
% xlim([0 T(end)])
% xlabel('Time [h]')
% ylabel('Heat generation from thermal solar panel (kW)')
% grid on
%% Power Consumption of HP
% figure()
% hold on
% plot(T,U(:,1),'r-','LineWidth',2)
% plot(T,U(:,2),'b-','LineWidth',2)
% %plot([T(1),T(end)],problem.inputs.uu(1)*ones(1,2),'r--','LineWidth',2)
% %plot([T(1),T(end)],problem.inputs.ul(1)*ones(1,2),'r--','LineWidth',2)
% %plot(xx,ct*100,'r-' )  % plot scaled prices
% % plot(tv,uv(:,1),'k-.')
% xlim([0 T(end)+step])
% xlabel('Time [h]')
% ylabel('Power Consumption of edHP (heating & cooling) (kW)')
% legend('HP (heating)','HP (cooling)')
% grid minor


% % figure(5)
% % hold on
% % plot(T,U(:,2),'g-','LineWidth',2 )
% % %plot([T(1),T(end)],problem.inputs.uu(2)*ones(1,2),'r--','LineWidth',2)
% % %plot([T(1),T(end)],problem.inputs.ul(2)*ones(1,2),'r--','LineWidth',2)
% % %plot(xx,ct*100,'r-' )  % plot scaled prices
% % % plot(tv,uv(:,1),'k-.')
% % xlim([0 T(end)+step])
% % xlabel('Time [h]')
% % ylabel('Power Consumption - Seasonal storage (KW)')
% % grid on

% figure
% hold on
% stairs(xx,ct,'r-','LineWidth',2)
% xlim([0 solution.tf])
% xlabel('Time [h]')
% ylabel('Flexible tariffs for residential customers £/kW')
%% Bought Power, Sold Power, DR Power
% figure(6)
% hold on
% plot(T,U(:,5),'b','LineWidth',2)
% % plot([T(1),T(end)],problem.inputs.uu(6)*ones(1,2),'b--','LineWidth',2)
% % plot([T(1),T(end)],problem.inputs.ul(6)*ones(1,2),'b--','LineWidth',2)
% plot(T,U(:,6),'c','LineWidth',2)
% % plot([T(1),T(end)],problem.inputs.uu(6)*ones(1,2),'b--','LineWidth',2)
% % plot([T(1),T(end)],problem.inputs.ul(6)*ones(1,2),'b--','LineWidth',2)
% plot(T,U(:,9),'g','LineWidth',2)
% xlim([0 T(end)+step])
% xlabel('Time [h]')
% ylabel('Power kW')
% legend( 'Power bought from Grid (Kw)','Sold power (kw)','DR (kW)')
% grid minor

%% Power provided to satisfy load - PV, Battery Discharge, Bought Power
% figure() %not really needed to show
% hold on 
% plot(xx,P_PV,'r','LineWidth',2)
% stairs(T,U(:,3),'g','LineWidth',2 )
% stairs(T,U(:,5),'r--','LineWidth',2)
% xlim([0 T(end)+step])
% xlabel('Time [h]')
% ylabel('Power provided to satisfy the load (kW)')
% legend( 'PV Power (kW)', 'Battery Discharged Power (kW)','Bought Power (kW)')
% grid minor

%% Price of electricity, Load vs Generation, Bought vs Sold vs DR, Penalty
%Load = HP heating/cooling + battery charge rate
%LoadBase=Ubase(:,1)+Ubase(:,2); %+U(:,4); %+U(:,2); SS
Load=U(:,1)+U(:,2); %+U(:,4); %+U(:,2); SS
%generationBase=P_PV+interp1(T,Ubase(:,3),xx,'linear',U(end,3));
generation=P_PV+interp1(T,U(:,3),xx,'linear',U(end,3));
%speval(solution,'U',4,xx);
sold=U(:,6);
Bought=U(:,5);
DR=U(:,9);

figure() 

% % Price of electricity
subplot(4,1,1)
plot(xx,ct,'LineWidth',2)
xlim([105 130]) 
%xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('Price (£/kWh)')
title('Electricity Tariff')
%legend( '£/h')
grid minor; box on;

subplot(4,1,2)
hold on
%plot(T,LoadBase,'k:','LineWidth',1.5)
plot(T,Load,'r','LineWidth',1.5)
xlim([105 130]) 
%xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('Power (kW)')
title('Load Demand, HP')
legend( 'Base', 'Case 1-B')
grid minor; box on;

subplot(4,1,3)
hold on
%plot(xx,generationBase,'k:','LineWidth',1.5)
plot(xx,generation,'b--','LineWidth',1.5)
xlim([105 130]) %T(end)+step
%xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('Power (kW)')
title('Generation, PV + u_d')
legend( 'Base', 'Case 1-B')
grid minor; box on;

subplot(2,1,1)
hold on
%plot(T,Ubase(:,5),'k:','LineWidth',2)
plot(T,Bought,'LineWidth',1.5) %,'m--',
%plot(T,sold,'LineWidth',1.5 ) %,'k--'
%plot(T,Ubase(:,9),'LineWidth',1.5) %,'g--',
plot(T,DR,'LineWidth',1.5) %,'g--',
xlim([105 130]) %T(end)+step
%xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('Power (kW)')
title('Power traded with grid')
legend('Bought power','DR')%,'Sold power')'Bought power (Base)'
grid minor; box on;

%PENALTY
% subplot(4,1,4)
% hold on
% plot(T,U(:,8),'LineWidth',1.5) %,'m--',
% xlim([0 T(end)+step])
% xlabel('Time [h]')
% ylabel('DR Penalty (kW)')
% legend('Penalized Energy (kW)')
% grid minor; box on;
% ylim([0 5]); 


%% All control input variables (check only)
figure(1)
hold on;
plot(T,U(:,1),'--','LineWidth',2)
figure(2)
hold on;
plot(T,U(:,2),'--','LineWidth',2)
figure(3)
hold on;
plot(T,U(:,3),'LineWidth',2)
figure(4)
hold on;
plot(T,U(:,4),'--','LineWidth',2)
figure(5)
hold on;
plot(T,U(:,5),'LineWidth',2)
figure(6)
hold on;
plot(T,U(:,6),'--','LineWidth',2)
figure(7)
hold on;
plot(T,U(:,7),'--','LineWidth',2)
figure(8)
hold on;
plot(T,U(:,8),':','LineWidth',2)
figure(9)
hold on;
plot(T,U(:,9),'--','LineWidth',2)
xlim([0 T(end)+step])
xlabel('Time [h]')
ylabel('All control input variables (kW)')
legend( 'u_1', 'u_2','u_d','u_c','u_B_{grid}','u_S_{grid}','u_{ADR}','s_t','u_{DR}')
grid minor; box on;
