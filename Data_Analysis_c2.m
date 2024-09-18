%Data Analysis based methods
function [MC_mean_error_perc_PBE, MC_mean_error_perc_PTBE, MC_mean_error_perc_CE] = Data_Analysis_c2(N,train_ratio)

min_d = 100; max_d = 180; %min. and maximum distances
d_vec = min_d:10:max_d;
load('T21.mat');
Data = table2array(T21,'VariableNames',{'No','t_low','C_low','t_high','C_high','Risetime','Delta_C','Gradient','t_peak1','C_peak1','Te','E','Distance','E2','P','P2'});
M = 5; %number of measurement in Data for each distance and each Te.

% E_T1 = 215.7849; %transmitted energy for Te = 0.25s - normal
% E_T2 = 498.3479; %transmitted energy for Te = 0.5s - normal
% E_T3 = 782.2501; %transmitted energy for Te = 0.75s - normal
E_T1 = 102.53; %transmitted energy for Te = 0.25s - noise level substracted
E_T2 = 234.2906; %transmitted energy for Te = 0.5s - noise level substracted
E_T3 = 493.4903; %transmitted energy for Te = 0.75s - noise level substracted

P_T1 = 7.3598; %transmitted power for Te = 0.25s - noise level substracted
P_T2 = 9.3666; %transmitted power for Te = 0.5s - noise level substracted
P_T3 = 11.0108; %transmitted power for Te = 0.75s - noise level substracted
% P_T1 = 15.6810; %transmitted power for Te = 0.25s - normal
% P_T2 = 17.6720; %transmitted power for Te = 0.5s - normal
% P_T3 = 17.4654; %transmitted power for Te = 0.75s - normal
% P_T1 = 1;P_T2 = 1;P_T3 = 1;
x = (d_vec)'; %d vector
xp = ([0 d_vec])';

% Avg matrix
Avg(1:length(d_vec),1) = 0.25;
Avg(1:length(d_vec),2) = min_d:10:max_d;

test_ratio = 1 - train_ratio + 0.0001; %If 0.0001 is not added, crossvalind does not properly work.
%Monte Carlo
for i_MC = 1:N 
    %Choose Test Data Ratio and Random Selection of Data index
    for i_tt = 1:M:length(Data)-M+1
        [Train(i_tt:i_tt+M-1,1) Test(i_tt:i_tt+M-1,1)] = crossvalind('HoldOut', M, test_ratio);   
    end
    
    %Training Data
    P_R = Data(Train,16); % received power over t_peak s P2
    t_peak = Data(Train,9); %t_peak

    %Test Data
    P_R_test = Data(Test,16);
    t_peak_test = Data(Test,9); %t_peak
    temp_d = Data(Test,13);%distance
    d = temp_d(temp_d <= max_d);
    
    %% received power average
    MT = length(find(Train(1:M) == 1)); %find the number of elements for each distance to be averaged
    j = 1;
    for i = 1:MT:length(P_R)-MT+1
        P_avg(j,1) = mean(P_R(i:i+MT-1));
        beta(j,1) = P_avg(j,1)/P_T1;
        j = j + 1;
    end
    Avg(:,3) = beta(1:length(d_vec),1);
    % h = figure;

    % Curve fitting
    y = Avg(:,3); % average data to be fitted
    fitopt = fitoptions('Method','NonlinearLeastSquares','Algorithm','Levenberg-Marquardt');
    [fun, gof, output] = fit(x, y,'exp1',fitopt);
    a1 = fun.a;
    b1 = fun.b;
%     h = figure;
%     p = plot(fun, x, y,'b*');
%     p(1).MarkerSize = 10;
%     xlabel('Distance (cm)','fontsize',12); hh = ylabel('$$\overline{P_R}$$ / $$\overline{P_T}$$');
%     set(hh,'Interpreter','latex','fontsize',14);
%     xlim([0 220]);

    %% t_peak average
    t_peak_avg = zeros(length(d_vec),1);
    j = 1;
    for i = 1:MT:length(t_peak)-MT+1
        t_peak_avg(j,1) = mean(t_peak(i:i+MT-1));
        j = j + 1;
    end
    Avg(:,4) = t_peak_avg(1:length(d_vec),1);

    yt = Avg(:,4); % average data to be fitted
    ytp = [0; yt]; % insertion of initial value
    [funt, goft, outputt] = fit(xp, ytp,'exp1',fitopt);
    a2 = funt.a;
    b2 = funt.b;
%     h = figure;
%     plot(funt, xp, ytp,'b*');
%     p(1).MarkerSize = 10;
%     xlabel('Distance (cm)','fontsize',12); hh = ylabel('$$\overline{t_{peak}}$$ (s)');
%     set(hh,'Interpreter','latex','fontsize',14);
%     xlim([0 220]);


    %% Power based estimation
    for j = 1:length(d)
        d_hat1(j,1) = (1/b1)*log(P_R_test(j)/(P_T1*a1));
    end
    
    for j = 1:length(d_vec)
        index = find(d == d_vec(j));
        mean_error_perc_PBE(i_MC,j) = nanmean(100*abs(d_hat1(index) - d(index))./d(j));
%         mean_est(i_MC,j) = nanmean(d_hat1(index));
%         std_est(i_MC,j) = std(d_hat1(index),'omitnan');  
    end
    
%     h = figure; 
%     errorbar(d_vec, mean_est, std_est,'-o','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
%     xlim([min_d-5 max_d+5]);
%     xlabel('Actual Distance (cm)'); ylabel('Estimated Distance (cm)');
%     % title('Power Based Estimation');
%     hold on; plot(d,d,'--*','LineWidth',1);
%     legend('Power Based Estimation','Actual Values');

    %% Average peak time based estimation
    for j = 1:length(d)
        d_hat2(j,1) = (1/b2)*log(t_peak_test(j)/a2);
    end
    
    for j = 1:length(d_vec)
        index = find(d == d_vec(j));
        mean_error_perc_PTBE(i_MC,j) = nanmean(100*abs(d_hat2(index) - d(index))./d(j));
%         mean_est(i_MC,j) = nanmean(d_hat1(index));
%         std_est(i_MC,j) = std(d_hat1(index),'omitnan');  
    end
    
%     h = figure; 
%     errorbar(d_vec, mean_est2, std_est2,'-o','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
%     xlim([min_d-5 max_d+5]);
%     xlabel('Actual Distance (cm)'); ylabel('Estimated Distance (cm)');
%     % title('Peak Time Based Estimation');
%     hold on; plot(d,d,'--*','LineWidth',1);
%     legend('Peak Time Based Estimation','Actual Values');

    %% Combined Estimation
    for j = 1:length(d)
        d_hat3(j,1) = (1/(b1 - b2))*log((P_R_test(j)*a2)/(P_T1*t_peak_test(j)*a1));
    end
    
    for j = 1:length(d_vec)
        index = find(d == d_vec(j));
        mean_error_perc_CE(i_MC,j) = nanmean(100*abs(d_hat3(index) - d(index))./d(j));
%         mean_est(i_MC,j) = nanmean(d_hat1(index));
%         std_est(i_MC,j) = std(d_hat1(index),'omitnan');  
    end
%     h = figure; 
%     errorbar(d_vec, mean_est3, std_est3,'-o','MarkerSize',5,'MarkerEdgeColor','red','MarkerFaceColor','red');
%     xlim([min_d-5 max_d+5]);
%     xlabel('Actual Distance (cm)'); ylabel('Estimated Distance (cm)');
%     % title('Combined Estimation');
%     hold on; plot(d,d,'--*','LineWidth',1);
%     legend('Combined Estimation','Actual Values');
end    

MC_mean_error_perc_PBE = nanmean(mean_error_perc_PBE);
MC_mean_error_perc_PTBE = nanmean(mean_error_perc_PTBE);
MC_mean_error_perc_CE = nanmean(mean_error_perc_CE);
% 
% plot(d_vec, MC_mean_error_perc_PBE, d_vec, MC_mean_error_perc_PTBE, d_vec, MC_mean_error_perc_CE);
% hold on;

end