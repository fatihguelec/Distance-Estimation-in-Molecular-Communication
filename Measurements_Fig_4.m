close all
clear 
clc

T = readmatrix('Flow_Rate_Data.xlsx');
Q = T(2:11, 7); %Flow rate
v_0 = T(2:11, 10); %Initial velocity

h = figure; hold on;
yyaxis left
eb = plot(Q,'-*','LineWidth',1.25); 
eb.Color = [1 0 0];  
set(gca,'ycolor',[1 0 0]) ; %make the left y-axis red 
ylim([1.5*10^-6 4*10^-6]);
ylabel('Flow Rate (m^3/s)');

yyaxis right
plot(v_0,'--o','Color',[0 0 1],'LineWidth',1.25);
% set(gca,'ycolor',[0.4660 0.6740 0.1880]) ; %make the left y-axis green 
set(gca,'ycolor',[0 0 1]) ; %make the left y-axis green 
ylim([8 12]); xlim([0.8 10.2]);
legend('Flow Rate','Initial Velocity');
ylabel('Initial Velocity (m/s)'); xlabel('Measurement Number');
grid on;

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h,'Measurement_plot','-dpdf','-r0')%save as pdf 
