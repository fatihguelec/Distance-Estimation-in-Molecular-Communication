%Neural Network Regression for Distance Estimation in MC
function MC_mean_error_perc_NNR = Neural_Network_Regression_c(N,train_ratio)

min_d = 100; max_d = 180; %min. and maximum distances
load('T21.mat');

Data = table2array(T21,'VariableNames',{'No','t_low','C_low','t_high','C_high','Risetime','Delta_C','Gradient','t_peak1','C_peak1','Te','E','Distance','E2','P','P2'});

features = [2 3 6:11 14]; %choose the features
% features = 2:12; %choose all the features

input = (Data(:,features))';
output = (Data(:,13))';
d = min_d:10:max_d;
hns = 1; %Hidden node size
%Monte Carlo
for i = 1:N
%     rng(0);
%     rng(4151941);
    net = fitnet(hns,'trainlm'); %Create Neural Network with xx hidden layers and xx training algorithm
    % Setup Division of Data for Training, Validation, Testing
    % For a list of all data division functions type: help nndivide
    net.divideFcn = 'dividerand';  % Divide data randomly
    net.divideMode = 'sample';  % Divide up every sample
    net.divideParam.trainRatio = train_ratio;
    net.divideParam.valRatio = (1 - train_ratio)/2;
    net.divideParam.testRatio = (1 - train_ratio)/2;
%      net.layers{1}.transferFcn = 'purelin';
%           net.layers{2}.transferFcn = 'purelin';
    
    net.trainParam.showWindow = false;
    [net,tr] = train(net,input, output); %Training
    pred = net(input(:,tr.testInd)); %Estimation
%     pred = pred_all(tr.testInd);
    output_test = output(:,tr.testInd);
    RMSE(i) = sqrt(perform(net, pred, output_test));  %RMSE
    
%     pred2 = net(input(tr.testInd)); %Estimation test
%     RMSE2(i) = sqrt(nanmean((pred-output_test).^2));  %RMSE
      
    for j = 1:length(d)
        index = find(output_test == d(j));
        mean_error_perc(i,j) = nanmean(100*abs(pred(index) - output_test(index))./d(j));
        mean_est(i,j) = nanmean(pred(index));
        std_est(i,j) = std(pred(index),'omitnan');  
    end    
end

RMSE_mean_NNR = mean(RMSE);
MC_mean_est_NNR = nanmean(mean_est);
MC_std_est_NNR = nanmean(std_est);
MC_mean_error_perc_NNR = nanmean(mean_error_perc);

end
