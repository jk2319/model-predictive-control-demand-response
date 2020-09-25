    % Read data relative to the European Case study


clear all


% File containing data

filename='data_problem.xlsx';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define problem data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two dimensional system: 1 thermal equation,  1 battery storage 

% State and input bounds for the battery
Cap_b=60;  %kwh
p_rate=10;  % Kw 


% State upper bound of the seasonal storage 
xupSS=12500;  % The solar wall installed in the ABC office can generate up to 125000kwh/year
 
u_b0=5;
xlow=[18,0];
xup=[22,Cap_b];



% 9 inputs 
% 1) Electricty (Power) used by edHP, 2) cooling HP, 3) Battery
% discharge rate,  4) Battery charge rate, 5) Power bought from the grid
% 6) Power inject into the grid, 7) Active DR, 8) S_t slack, 9) total DR
 
ulow=[0, 0, 0, 0, 0, 0, 0, 0, 0];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carbon Emission half hourly data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Carbon_e=xlsread(filename,5,'B2:B17522')';  
T_carb=[0:1/2:8760];  % Vector of times

Price_carbon=100;   % in £/(ton CO2e)

Price_carbon= Price_carbon/1e+6; % converted in £/gr         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hourly data of Global horizontal irradiance (W/m^2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ir=xlsread(filename,6,'D18:D8778')'; 
           
% Converting W to kW
Ir=Ir/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve parameters ofa 3 room building for a lumped model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UV=xlsread(filename,1,'B4');     % overall heat transfer coefficient, U (Average U-value) ( W/m^2 K) 
UV=UV/1000;  % convertin W in kW
A=xlsread(filename,1,'B3');     % Wall surface Area (m^2)
rhoAir=xlsread(filename,1,'B8');     % Air density  (kg/ m^3)
B_V=xlsread(filename,1,'B6');     % Building Volume (m^3)
C_air=xlsread(filename,1,'B9');     % Air heat capacity (kJ/ kg/ K)
C_air=C_air/3600;  % Converting kJ in kWh
n_ac=xlsread(filename,1,'B7');     % Air changes per hour (hr^(-1))
C_Build=xlsread(filename,1,'B10');     % Building Capacity (kj/ K)
C_Build=C_Build/3600;    % Converting kJ in kWh
AFl=xlsread(filename,1,'B2');     % Floor surface Area (m^2)

uB_up=100;
ADR_max=3;
DR_min=1;
% Inputs upper bounds in kW
upSS=15*(UV*A+rhoAir*B_V*C_air*n_ac)/5.5;
uup=[7, 15, p_rate, p_rate, uB_up, 0, ADR_max, 100, 200]; %u_s_grid=0, off
% [HP(heating), HP(cooling), ud, uc, u_B_grid, u_s_grid, u_ADR, s_t, u_DR]

% Thermal solar and seasonal storage parameters

Conv_S_S=5.5;             % Conversion efficency from electricty to thermal of seasonal storage 
eff_T_S=70/100;           % Thermal solar efficiency  max should be 75% on avarage is around 65%
Area_T_S=25;              % Area (m^2) occupied by Thermal storage (assumption: thermal storage has been installed on walls.  
eff_SS=40/100;            % efficiency of seasonal storage   
K_sat=100^2;              % Saturation constant for the seasoal storage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exogenous Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 15 minutes flexible tariffs for residential customers in 2018 
% (Octopus Energy bases their prices on Nordpool day ahead option price (https://www.nordpoolgroup.com/))


Pr_elec=xlsread(filename,4,'B2:B35042')';  % (£/kWh)   From 1/1/2018   t0  end 31/12/2019  
% 15 minutes external temperature Source: National Centers for Environmental Information (National Oceanic and atmospheric administration)
%https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs

% External temperature
T_ex=xlsread(filename,2,'B2:B35042')';   % (^circ C)  From 1/1/2018   t0  end 31/12/2019 
T_ex2=xlsread(filename,2,'B1200:B36240')';   % (^circ C)  From 1/1/2018   t0  end 31/12/2019 


T_mex=[0:1/4:8760];  %Vector of times associated to the temperature measurements


% Interpolation values of Carbon emission
Carbon_em=interp1(T_carb,Carbon_e,T_mex,'linear',Carbon_e(end));



% Interpolation of of radiance data
Ir=interp1([0:8759],Ir,T_mex,'linear',Ir(end));

% Window 1000-1300 & 1800-2100
gamma1=xlsread(filename,7,'C2:C35042')';
% Calls 1100-1200, 1900-2000
alpha1=xlsread(filename,8,'C2:C35042')'; %2H
% Rand Calls // 1045, 1200, 1215, 1230 // 1845,2000,2015,2045
%alpha1=xlsread(filename,9,'C2:C35042')';

%alpha1=zeros(1,35041); %no calls case
% for i=1:35041 
%     if (gamma1(i)==1)
%         alpha1(i)=rand<7/20;
%     elseif (nnz(alpha1)/nnz(gamma1))>=(1/3)
%     end
% end

%writematrix(alpha1,filename,'Sheet',10,'Range','B2:B35042');
% Notes
%
% To Convert kJ to Wh divide by 3.6
%
% Kelvin to celsius conversion   [Â°C] = [K] âˆ’ 273.15 
%
%  The COP  for an electric driven heat pump at 7 \circ celsius is assumed equal to  3  
%
%  I assume that at 22 \circ celsius COP is 4 


COP_init=3;
T_init=7;
COP_Ta=4;
Ta=22;

% coefficient of the COP for the heat pump accounting of a linear dependence with respect to 
% the external temperature
m_COP=(COP_Ta-COP_init)/(Ta-T_init);


m_cop_cool=0.7;    % COP of the cooling system

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batteries 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dis_eff=1;     %discharge efficiency
ch_eff=0.88;   %charge efficiency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The prediction output is in KWp/m^2 

t1=xlsread(filename,3,'E4');
t2=xlsread(filename,3,'E5');
t3=xlsread(filename,3,'E6');


% Proportion of total area covered by PVs
PV_s_Ap=1/2;
% Proportion of total area covered by Thermal Solar
TSo_Ap=1/2; 


save External_data Ir Carbon_em Pr_elec T_ex T_mex gamma1 alpha1 T_ex2

save case_study_sim
