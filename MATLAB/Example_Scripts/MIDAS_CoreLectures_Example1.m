%% CORE LECTURE SERIES ON MIDAS 
% EXAMPLE 1: MACRO FORECASTING 
%Purpose:
%   To follow along with the macro forecasting slides for the Core Lecture
%   Series on MIDAS given by Eric Ghysels.
%Description:
%   The following script goes through a number of macro forecasting
%   examples using MIDAS regressions. MIDAS models are compared to time
%   averaged ones both in and out-of-sample.
%Information: 
%   Author - Anessa Custovic (acustovic12@gmail.com) | 21Aug2019
%   MIDAS functions called were written by Hang Qian.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%clear workspace
clear all
close all

% add scripts to filepath
addpath(genpath('MIDASLecture/Example_Scripts'));
% add MIDAS toolbox to filepath
addpath(genpath('MIDASLecture/MIDASv2.2/MIDASv2.2/privatex'));
% add data to filepath
addpath(genpath('MIDASLecture/Example_Data'));

% Load data for macro examples:
load('example1.mat')


% log growth
logpayems = log(payems(2:end,2)./payems(1:end-1,2))*100;
logrgdp = log(rgdp(2:end,2)./rgdp(1:end-1,2))*100;

% convert to datetime
datepayems= datetime(payems(2:end,1),'ConvertFrom','datenum','Format','yyyy-MM-dd');
datergdp= datetime(rgdp(2:end,1),'ConvertFrom','datenum','Format','yyyy-MM-dd');

% Specify dates
EstStart = '1985-01-01';
EstEnd = '2009-01-01';

% Specify lag structure, horizon and method
Xlag = 9;
Ylag = 1;
Horizon = 3;
Method = 'fixedWindow';

%set data- x and y variables
DataX=logpayems;
DataY=logrgdp;

%set data dates
DataYdate = datergdp;
DataXdate = datepayems;

poly='expAlmon';

[OutputForecast1,OutputEstimate1,MixedFreqData] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);
%MIDAS weights in OutputEstimate are actually high frequency slope*weights
%so they won't sum to 1. If you want weights that sum to one you can do
%weights./slope to recover original


%now compare to time averaged employment data
%form X vars
Xvars=[MixedFreqData.EstLagY, mean(MixedFreqData.EstX(:,1:3),2), mean(MixedFreqData.EstX(:,4:6),2), mean(MixedFreqData.EstX(:,7:9),2)];
est_ta= fitlm(Xvars,MixedFreqData.EstY)

%compare in sample performance
sqrt(mean(OutputEstimate1.resid.^2))

sqrt(mean(table2array(est_ta.Residuals(:,1)).^2))


% now compare out of sample
% first predict time-averaged model
Xnew=[MixedFreqData.OutLagY, mean(MixedFreqData.OutX(:,1:3),2), mean(MixedFreqData.OutX(:,4:6),2), mean(MixedFreqData.OutX(:,7:9),2)];
ta_oos = predict(est_ta,Xnew);

% compare RMSE
% grab MIDAS RMSE
midas_oos_rmse = OutputForecast1.RMSE
ta_oos_rmse = sqrt(mean((ta_oos-MixedFreqData.OutY).^2))



%plot weights and compare
%set lags
lags=1:9;
%time averaged weights
ta_weights = [repelem(table2array(est_ta.Coefficients(3,1)),3) repelem(table2array(est_ta.Coefficients(4,1)),3) repelem(table2array(est_ta.Coefficients(5,1)),3)]';  
%now plots... 
    subplot(2,1,1)
    plot(lags,OutputEstimate1.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('Normalized exponential Almon lag polynomial');

    subplot(2,1,2)
    plot(lags,ta_weights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('Time-averaged data coefficients');
    
    
    
    
%compare fixed, rolling, expanding windows for predictions
%fixed windows first..


[OutputForecast1_Fixed,OutputEstimate1_Fixed,MixedFreqData_Fixed] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','expAlmon','Method',Method,'PlotWeights',0);

%Now time averaged MIDAS- IS
Xvars_TA=[MixedFreqData_Fixed.EstLagY, mean(MixedFreqData_Fixed.EstX(:,1:3),2), mean(MixedFreqData_Fixed.EstX(:,4:6),2), mean(MixedFreqData_Fixed.EstX(:,7:9),2)];
est_ta_Fixed= fitlm(Xvars_TA,MixedFreqData_Fixed.EstY);

%OOS
Xnew_TA_Fixed=[MixedFreqData_Fixed.OutLagY, mean(MixedFreqData_Fixed.OutX(:,1:3),2), mean(MixedFreqData_Fixed.OutX(:,4:6),2), mean(MixedFreqData_Fixed.OutX(:,7:9),2)];
ta_oos_Fixed = predict(est_ta_Fixed,Xnew_TA_Fixed);


%rolling windows 
Method='RollingWindow';
[OutputForecast1_Roll,OutputEstimate1_Roll,MixedFreqData_Roll] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','expAlmon','Method',Method,'PlotWeights',0);

%lags=1:9;
plot(lags,OutputEstimate1_Roll.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('Normalized exponential Almon lag polynomial');


%rolling time averaged
%Now time averaged MIDAS- IS
nobs=length(MixedFreqData_Roll.EstY);
nforecast=length(MixedFreqData_Roll.OutY);
nroll=nforecast;

Xvars_TA=[MixedFreqData_Roll.EstLagY, mean(MixedFreqData_Roll.EstX(:,1:3),2), mean(MixedFreqData_Roll.EstX(:,4:6),2), mean(MixedFreqData_Roll.EstX(:,7:9),2)];
Xnew_TA_Roll=[MixedFreqData_Roll.OutLagY, mean(MixedFreqData_Roll.OutX(:,1:3),2), mean(MixedFreqData_Roll.OutX(:,4:6),2), mean(MixedFreqData_Roll.OutX(:,7:9),2)];
Xvars_temp = [Xvars_TA;Xnew_TA_Roll];
Yvars_temp = [MixedFreqData_Roll.EstY;MixedFreqData_Roll.OutY];  

for t=1:nroll
    Xvars_TA_Roll = Xvars_temp(t:nobs-1+t,:);
    Yvars_Roll = Yvars_temp(t:nobs-1+t,:);
    est_ta_Roll = fitlm(Xvars_TA_Roll,Yvars_Roll);
    Xvars_OOS_Roll = Xnew_TA_Roll(t,:);
    ta_oos_Roll(t,:) = predict(est_ta_Roll,Xvars_OOS_Roll);
end

%compare forecasts of rolling
% compare RMSE
% grab MIDAS RMSE
OutputForecast1_Roll.RMSE
sqrt(mean((ta_oos_Roll-MixedFreqData_Roll.OutY).^2))


%Expanding Window
Method='Recursive';
[OutputForecast1_Exp,OutputEstimate1_Exp,MixedFreqData_Exp] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','expAlmon','Method',Method,'PlotWeights',0);

%lags=1:9;
plot(lags,OutputEstimate1_Exp.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('Normalized exponential Almon lag polynomial');



%Expanding time averaged
%Now time averaged MIDAS- IS
nobs=length(MixedFreqData_Exp.EstY);
nforecast=length(MixedFreqData_Exp.OutY);
Xvars_TA_Exp=[MixedFreqData_Exp.EstLagY, mean(MixedFreqData_Exp.EstX(:,1:3),2), mean(MixedFreqData_Exp.EstX(:,4:6),2), mean(MixedFreqData_Exp.EstX(:,7:9),2)];
Xnew_TA_Exp=[MixedFreqData_Exp.OutLagY, mean(MixedFreqData_Exp.OutX(:,1:3),2), mean(MixedFreqData_Exp.OutX(:,4:6),2), mean(MixedFreqData_Exp.OutX(:,7:9),2)];
Xvars_temp_Exp = [Xvars_TA_Exp;Xnew_TA_Exp];
Yvars_temp_Exp = [MixedFreqData_Exp.EstY;MixedFreqData_Exp.OutY];  

for t=1:nroll

    Xvars_TA_Exp = Xvars_temp_Exp(1:nobs-1+t,:);
    Yvars_Exp = Yvars_temp_Exp(1:nobs-1+t,:);
    est_ta_Exp = fitlm(Xvars_TA_Exp,Yvars_Exp);
    Xvars_OOS_Exp = Xnew_TA_Exp(t,:);
    ta_oos_Exp(t,:) = predict(est_ta_Exp,Xvars_OOS_Exp);
    
end

% compare forecasts of rolling
% compare RMSE
% grab MIDAS RMSE
OutputForecast1_Exp.RMSE
ta_oos_rmseExp = sqrt(mean((ta_oos_Exp-MixedFreqData_Exp.OutY).^2))







%%second example with cfnai data
%convert cfnai dates to datetime
datecfnai= datetime(cfnai(1:end,1),'ConvertFrom','datenum','Format','yyyy-MM-dd');

% Specify lag structure and sample size 
Xlag = 12;
Ylag = 1;
Horizon = 3;
EstStart = '1987-01-01';
EstEnd = '2011-12-01';
Method = 'fixedWindow';
%Data
DataYdate = datergdp;
DataXdate = datecfnai;
DataX=cfnai(:,2);
DataY=logrgdp;

%Estimate ADL regression model using exp almon scheme
poly='expAlmon'; 
[OutputForecast1,OutputEstimate1,MixedFreqData] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);


%Estimate ADL regression model using U-MIDAS scheme
poly='UMIDAS';
[OutputForecast_U,OutputEstimate_U,MixedFreqData_U] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);


%compare in sample performance
sqrt(mean(OutputEstimate1.resid.^2)) 
sqrt(mean(OutputEstimate_U.resid.^2))


%compare OOS performance now
%grab MIDAS RMSE
midas_oos_rmse = OutputForecast1.RMSE
umidas_oos_rmse = OutputForecast_U.RMSE



%plot weights and compare
%set lags
lags=1:12;

    subplot(2,1,1)
    plot(lags,OutputEstimate1.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('Normalized exponential Almon lag polynomial');

    subplot(2,1,2)
    plot(lags,OutputEstimate_U.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('U-MIDAS lag polynomial');
    
    
    
    
    
    
    
    
    
%Third Example now with ADS data
dateads= datetime(ads(1:end,1),'ConvertFrom','datenum','Format','yyyy-MM-dd');

% Specify lag structure and sample size 
Xlag = 66;
Ylag = 1;
Horizon = 66;
EstStart = '1987-01-01';
EstEnd = '2011-12-01';
Method = 'fixedWindow';

DataYdate = datergdp;
DataXdate = dateads;
DataX=ads(:,2);
DataY=logrgdp;

%let's estimate with exp almon polynomial first
poly='expAlmon'; 
[OutputForecast1,OutputEstimate1,MixedFreqData] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);

%Estimate ADL regression model using U-MIDAS scheme
%these match
poly='UMIDAS';
[OutputForecast_U,OutputEstimate_U,MixedFreqData_U] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);



%nowcasts
Horizon = 22;
%first, expAlmon
poly='expAlmon'; 
[OutputForecast1_Now,OutputEstimate1_Now,MixedFreqData_Now] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);

%now Umidas
poly='UMIDAS';
[OutputForecast_UNow,OutputEstimate_UNow,MixedFreqData_UNow] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial',poly,'Method',Method,'PlotWeights',0);



%forecasts OOS
midas_oos_rmsefix = OutputForecast1.RMSE
umidas_oos_rmsefix = OutputForecast_U.RMSE

%nowcasts
midas_oos_rmsenow = OutputForecast1_Now.RMSE
umidas_oos_rmsenow = OutputForecast_UNow.RMSE



%plot weights and compare
%set lags
lags=1:66;

    subplot(2,1,1)
    plot(lags,OutputEstimate1.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('Normalized exponential Almon lag polynomial');

    subplot(2,1,2)
    plot(lags,OutputEstimate_U.estWeights)
    xlabel('Lags'); 
    ylabel('Coefficient');
    title('U-MIDAS lag polynomial');
    
    
    
    

