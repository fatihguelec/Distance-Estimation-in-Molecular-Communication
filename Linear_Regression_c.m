%Linear Regression for Distance Estimation in MC
function MC_mean_error_perc_LR = Linear_Regression_c(N,train_ratio)

min_d = 100; max_d = 180; %min. and maximum distances
load('T21.mat');

Data = table2array(T21,'VariableNames',{'No','t_low','C_low','t_high','C_high','Risetime','Delta_C','Gradient','t_peak1','C_peak1','Te','E','Distance','E2','P','P2'});

features = [2 3 6:11 14]; %choose the features
% features = [2 3 9 10 11]; %choose the features
% features = 2:12; %choose all the features
d = min_d:10:max_d;
test_ratio = 1 - train_ratio + 0.0001; %If 0.0001 is not added, crossvalind does not properly work.
%Monte Carlo
for i = 1:N
    %Choose Test Data Ratio and Random Selection of Data index
    [Train, Test] = crossvalind('HoldOut', length(Data), test_ratio);   
    
    %Training Data
    X_train = Data(Train, features);
    Y_train = Data(Train,13);
%     XY_train = [X_train, Y_train];

    %Test Data
    X_test = Data(Test,features);
    Y_test = Data(Test,13);
%     XY_test = [X_test, Y_test];
%     XY_all = [XY_train; XY_test];

    lmModel = fitlm(X_train, Y_train, 'linear', 'RobustOpts', 'off');

    Y_pred = predict(lmModel, X_test);
    RMSE(i) = sqrt(mean((Y_pred - Y_test).^2));%RMSE-or cost function
%     Y = [Y_test Y_pred];
    
    %Error for each distance
%     for j = 1:length(d)
%         index = find(Y_test == d(j));
% %         RMSE_d(i,j) = sqrt(mean((Y_pred(index) - Y_test(index)).^2));
%         error_perc(i,j) = nanmean(abs((Y_pred(index) - Y_test(index))))*100 / d(j); %tekrar bak
%     end   
    for j = 1:length(d)
        index = find(Y_test == d(j));
        mean_error_perc(i,j) = nanmean(100*abs(Y_pred(index) - Y_test(index))./d(j));
        mean_est(i,j) = nanmean(Y_pred(index));
        std_est(i,j) = std(Y_pred(index),'omitnan');  
    end
    
end
% load('LR_1e4.mat');
MC_mean_est_LR = nanmean(mean_est);
MC_std_est_LR = nanmean(std_est);
MC_mean_error_perc_LR = nanmean(mean_error_perc);
% mean_error_perc = 100*abs(MC_mean_est-d)./d;
RMSE_mean_LR = mean(RMSE);

end

%% Least Squares Method - How the coefficients are calculated
% X = [ones(length(X_train),1) X_train];
% Y = Y_train;
% Theta = ((X'*X)^-1)*X'*Y;



