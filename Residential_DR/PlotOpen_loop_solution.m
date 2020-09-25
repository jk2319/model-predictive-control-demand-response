


%% figure
par=problem.data;

xx=linspace(solution.T(1,1),solution.tf,10000);
ct=interp1(data.data.T_mex,data.data.Pr_elec ,xx,'previous',data.data.Pr_elec(end));
T_outside=interp1(data.data.T_mex,data.data.T_ex,xx,'linear',data.data.T_ex(end));
I_rad=interp1(data.data.T_mex,data.data.Ir,xx,'linear',data.data.Ir(end));




figure(1)
hold on
plot(simtime+solution.T,solution.X(:,1),'Color',[0,0.7,0.9],'LineWidth',2)
plot(simtime+xx,T_outside,'k-','LineWidth',1)
plot([solution.T(1),simtime+solution.T(end)],problem.states.xu(1)*ones(1,2),'r--','LineWidth',1)
plot([solution.T(1),simtime+solution.T(end)],problem.states.xl(1)*ones(1,2),'r--','LineWidth',1)
% plot(tv,xv(:,1),'k-.')
% plot(tv,xv(:,2),'k-.')
xlim([0 simtime+solution.tf])
xlabel('Time [h]')
ylabel('State')
legend('Room Temperature [\circ C]',' External Temperature [\circ C]')
grid on


P_TS=data.data.eff_T_S*data.data.Ir*data.data.TSo_Ap*data.data.AFl;   % Power from thermal solar
% Power produced by PVs
P_PV=(data.data.t1*data.data.Ir+ data.data.t2*data.data.Ir.^2+data.data.t3*data.data.Ir.*data.data.T_ex)*data.data.PV_s_Ap*data.data.AFl;

% figure(2)
% hold on
% pp=plot(simtime+solution.T,solution.X(:,2),'Color',[0,0.7,0.9],'LineWidth',2);
% pp(1).Marker = '*';
% plot(simtime+solution.T,solution.U(:,2),'r','LineWidth',2)
% plot(simtime+data.data.T_mex,P_TS,'k--','LineWidth',1)
% %plot([T(1),T(end)],problem.states.xu(2)*ones(1,2),'r--','LineWidth',2)
% %plot([T(1),T(end)],problem.states.xl(2)*ones(1,2),'r--','LineWidth',2)
% xlim([0 simtime+solution.tf])
% legend( 'State of charge of seasonal storage (Kwh)','Power-discharging (Kw)')
% xlabel('Time [h]')
% ylabel('Seasonal Storage')
% grid on

%modified
figure(3)
hold on
plot(simtime+solution.T,solution.X(:,2),'Color',[0,0.7,0.9],'LineWidth',2)
plot(simtime+solution.T,solution.U(:,4),'Color',[0.7,0.0,0.9],'LineWidth',2);
plot(simtime+solution.T,-solution.U(:,3),'Color',[0.7,0.9,0],'LineWidth',2);
%plot([solution.T(1),solution.T(end)],problem.states.xu(3)*ones(1,2),'b--','LineWidth',2)
plot([solution.T(1),solution.T(end)],-problem.inputs.uu(3)*ones(1,2),'k--','LineWidth',2)
% %plot([T(1),T(end)],problem.inputs.ul(5)*ones(1,2),'k--')
% plot([T(1),T(end)],problem.inputs.uu(5)*ones(1,2),'k--','LineWidth',2)
xlim([0 simtime+solution.tf])
xlabel('Time [h]')
ylabel('State of charge of the Battery')
%legend( 'Energy-Battery (Kwh)',  'Power-Battery-discharge/charge (Kw)')
%legend( 'Energy-Battery (Kwh)',  'Power-Battery-discharge (Kw)','Power-Battery-charge (Kw)')
grid on

%modified
figure(4)
hold on
plot(simtime+solution.T,solution.U(:,1),'Color',[0,0.7,0.9],'LineWidth',2)
plot(simtime+solution.T,solution.U(:,2),'r','LineWidth',2)
plot([solution.T(1),solution.T(end)],problem.inputs.uu(1)*ones(1,2),'r--','LineWidth',2)
plot([solution.T(1),solution.T(end)],problem.inputs.ul(1)*ones(1,2),'r--','LineWidth',2)
xlim([0 simtime+solution.tf])
xlabel('Time [h]')
ylabel('Power Consumption - edHP and Cooling (KW)')
grid on


% figure(5)
% hold on
% plot(simtime+solution.T,solution.U(:,2),'Color',[0,0.7,0.9],'LineWidth',2 )
% plot([solution.T(1),solution.T(end)],problem.inputs.uu(2)*ones(1,2),'r--','LineWidth',2)
% plot([solution.T(1),solution.T(end)],problem.inputs.ul(2)*ones(1,2),'r--','LineWidth',2)
% xlim([0 simtime+solution.tf])
% xlabel('Time [h]')
% ylabel('Power Consumption - Seasonal storage (KW)')
% grid on

figure(9)
hold on
stairs(simtime+xx,ct,'r-','LineWidth',2)
xlim([0 simtime+solution.tf])
xlabel('Time [h]')
ylabel('Flexible tariffs for residential customers £/kW')

%modified
figure(6)
hold on
plot(solution.T,solution.U(:,5),'b','LineWidth',2)
plot(solution.T,solution.U(:,6),'r','LineWidth',2)
plot(solution.T,solution.U(:,9),'g','LineWidth',2)
plot([solution.T(1),solution.T(end)],problem.inputs.uu(5)*ones(1,2),'b--','LineWidth',2) %uup
plot([solution.T(1),solution.T(end)],problem.inputs.uu(6)*ones(1,2),'r--','LineWidth',2)
%plot([solution.T(1),solution.T(end)],problem.inputs.uu(9)*ones(1,2),'g--','LineWidth',2)
xlim([0 solution.T(end)])
xlabel('Time [h]')
ylabel('Power kW')
legend('Power bought from Grid (Kw)','Sold power (kw)','DR (kW)')

%modified
figure(7)
hold on 
plot(simtime+data.data.T_mex,P_PV,'r','LineWidth',2)
plot(simtime+solution.T,solution.U(:,5),'b','LineWidth',2 )
plot(simtime+solution.T,solution.U(:,3),'k','LineWidth',2 )
plot(simtime+solution.T,solution.U(:,6),'r--','LineWidth',2)
plot(simtime+solution.T,solution.U(:,9),'g','LineWidth',2)
xlim([0 simtime+solution.tf])
xlabel('Time [h]')
ylabel('Power provided  to satisfy the load kW')
legend( 'Power-PV (Kw)','Power from Grid (Kw)',  'Power-Battery-discharge (Kw)','Sold power (kw)','DR (kW)')

%Load = HP heating & cooling + battery charging
Load=solution.U(:,1)+solution.U(:,2)+solution.U(:,4); %+solution.U(:,2)
sold=solution.U(:,6);
Bought=solution.U(:,5);
DR=solution.U(:,9);

figure(8) 
hold on
plot(simtime+solution.T,Load,'r','LineWidth',2)
plot(simtime+solution.T,sold,'k--','LineWidth',2 )
plot(simtime+solution.T,Bought,'m','LineWidth',2 )
plot(simtime+solution.T,DR,'g','LineWidth',2 )
xlim([0 simtime+solution.tf])
xlabel('Time [h]')
ylabel('Power, Load, Sold power, DR (kW)')
legend( 'Load Demand (Kw)','Sold power (kw)','Bought power (kw)','DR Scheme (kw)')
