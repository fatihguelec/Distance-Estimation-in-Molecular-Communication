% Comparison of FDDE, data analysis-based and machine learning methods
close all
clear
clc

N = 1e4; %Number of iterations
train_ratio = 0.8; %Training ratio
min_d = 100; max_d = 180; %min. and maximum distances

[MC_mean_error_perc_PBE, MC_mean_error_perc_PTBE, MC_mean_error_perc_CE]  = Data_Analysis_c2(N, train_ratio); % Obtain results for data analysis based methods
MC_mean_error_perc_LR  = Linear_Regression_c(N, train_ratio); % Obtain results for linear regression
MC_mean_error_perc_NNR  = Neural_Network_Regression_c(N, train_ratio); % Obtain results for neural network regression

%% FDDE
% Parameters at 20 C
format short;
ts = 50; %Total simulation time (s)
% Te = [0.25 0.5 0.75]; %Emission time (s)
Te = 0.25;
delta_t = 10^-2; %Time step (s)
d0 = 510*10^-6; %initial droplet diameter (m)
v0 = 10.8; %Initial velocity of droplets (m/s)
theta = deg2rad([38 33 28]); %half-angle of the TX beamwidth (degrees)

rho_d = 789; %Droplet density - liquid ethanol (kg/m^3)
M_v = 46.069*10^-3; % Molecular weight of ethanol (kg/mol)
M_a = 28.9647*10^-3; %Molecular weight of air (kg/mol)
D_v = 11.81*10^-6; %Diffusion coefficient of ethanol (m^2/s)
rho_a = 1.2; % Density of air (kg/m^3)
rho_lg = 1.59*rho_a; %Droplet density - vapor ethanol (kg/m^3)
P_a = 10^5; %Partial pressure of air (kg/ms^2) (10^5 Pa)
P_v = 5800; %Ethanol vapor pressure (Pa)
delta_P = 790; %Vapour pressure difference (kg/m s^2) (790 Pa)??
T = 20+273.15; %Temperature in Kelvin
nu_a = 1.516*10^-5; % Kinematic viscosity of air
S_c = nu_a/D_v; %Schmidt number
mu_a = 17.83*10^-6; %dynamic viscosity of air (N s/ m^2)
mu_l = 1.83*10^-5; %dynamic viscosity of ethanol (N s/ m^2)
Re = 0; %Reynolds number
rho_a_t = rho_a/rho_d;
D0 = d0; %nozzle diameter
alpha_d0 = 10^-3; %volume fraction of droplets in the mixture

%% Mixture Two Phase Flow of Evaporating Droplets
t = 0:delta_t:ts;
% Changing diameter model with Eq.(39) in Sazhin et.al. 2001
% diameters of droplets(d) and volume fraction of droplets(alpha_d) until ts

for k = 1:length(theta)
    d(1) = d0; %initial droplet diameter
    alpha_d(1) = 0;
    s(k,1) = 0; 
    for i = 1:length(t)-1
        if t(i) <= Te
            d(i+1) = d(i);
            alpha_d(i+1) = t(i+1)*alpha_d0/Te;  
        else
            delta_d = -2*((M_v*D_v*rho_a*delta_P)*(2+0.6*(S_c^(1/3))*(Re^(1/3)))*delta_t/(M_a*d(i)*rho_d*P_a));
            d(i+1) = d(i) + delta_d;
            if d(i+1) <= 0
                d(i+1) = 0;
                s(k,i+1:length(t)) = s(k,i);
                break;
            end
            alpha_d(i+1) = alpha_d(i) * (d(i+1)/d(i))^3;
        end
        k3 = (16*(1-alpha_d(i+1))*rho_a*(tan(theta(k)))^2) / (D0^2*rho_d);
        k2 = (16*(1-alpha_d(i+1))*rho_a*tan(theta(k))) / (D0*rho_d);
        k1 = (1+4*(1-alpha_d(i+1))*rho_a)/rho_d;
        a1 = k3;        
        b1 = k2 - 2*k3*s(k,i); 
        c1 = k1 - 2*k2*s(k,i) + k3*s(k,i)^2 - 1; 
        f1 = -2*k1*s(k,i) + k2*s(k,i)^2 + 4*v0*delta_t + 2*s(k,i);
        g1 = k1*s(k,i)^2 - 4*v0^2*delta_t^2 - 4*v0*delta_t*s(k,i) - s(k,i)^2;
        s_roots = quartic_roots(a1, b1, c1, f1, g1);
        s(k, i+1) = max(s_roots(real(s_roots) > 0 & abs(imag(s_roots)) < 0.0000001 )); %find the real and positive root
    end
end

%% Experimental data only for Te=0.25
load('T21.mat');
Data = table2array(T21,'VariableNames',{'No','t_low','C_low','t_high','C_high','Risetime','Delta_C','Gradient','t_peak1','C_peak1','Te','E','Distance','E2','P','P2'});
Te_vec = Data(:,11); %Te
t_peak = Data(:,9); %t_peak
dist = 0.01.*Data(:,13);%distance (meter)
d_vec = 1:0.1:1.8;
x = d_vec';

%% Distance comparison
i = 1;
for m = 1:length(theta)
    for i = 1:length(t_peak)
        s_est(m,i) = real(s(m,max(find(abs(t(1,:) - t_peak(i)) < 0.01))));
    end
end

s_est = s_est(:,1:45);%Choose the estimated values for 100-180 cm
% dist_p = (dist(1:45))';

for m = 1:length(theta)
    i = 1;
    for j = 1:5:length(s_est)-4
        s_est_mean(m,i) = mean(s_est(m,j:j+4));
        s_est_std(m,i) = std(s_est(m,j:j+4));
        i = i+1;
    end
end

%% Mean Absolute Percentage Error (MAPE) comparison of distance
for m = 1:length(theta)
    i = 1;
    for j = 1:5:length(s_est)-4
        mape_FDDE(m,i) = 100*mean(abs(s_est(m,j:j+4) - d_vec(i))/d_vec(i));
        i = i+1;
    end
end

h = figure;
% plot(d_vec, mape_FDDE(1,:),'k-o',d_vec, mape_FDDE(2,:),'m:+',d_vec, mape_FDDE(3,:),'r-x','LineWidth',2); 
plot(d_vec, mape_FDDE(3,:),'r-x','LineWidth',2); 
hold on;

plot(d_vec, MC_mean_error_perc_LR,'-s','MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m')
hold on;
plot(d_vec, MC_mean_error_perc_NNR,'-.+','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k')
hold on;
plot(d_vec, MC_mean_error_perc_PBE,'-o','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on;
plot(d_vec, MC_mean_error_perc_PTBE,'--*','MarkerSize',5,'MarkerEdgeColor','g','MarkerFaceColor','g')
hold on;
plot(d_vec, MC_mean_error_perc_CE,'-.d','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b')

xlim([min_d/100-0.05 max_d/100+0.05]);
ylim([0 50]);
% legend(['FDDE (\theta=',num2str(rad2deg(theta(1))),char(176),')'],['FDDE (\theta=',num2str(rad2deg(theta(2))),char(176),')'],['FDDE (\theta=',num2str(rad2deg(theta(3))),char(176),')'],'Linear Regression','Neural Network Regression','Power Based Estimation','Peak Time Based Estimation','Combined Estimation');
legend(['FDDE (\theta=',num2str(rad2deg(theta(3))),char(176),')'],'Linear Regression','Neural Network Regression','Power Based Estimation','Peak Time Based Estimation','Combined Estimation');

grid on;
xlabel('Actual Distance (m)'); ylabel('Mean Absolute Percentage Error (%)'); 

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'MAPE_comp_plot4','-dpdf','-r0')%save as pdf
