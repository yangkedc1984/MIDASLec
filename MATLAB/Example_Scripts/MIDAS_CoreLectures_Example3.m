%% CORE LECTURE SERIES ON MIDAS 
% EXAMPLE 3: GARCH-MIDAS 
%Purpose:
%   To follow along with the GARCH-MIDAS slides for the Core Lecture
%   Series on MIDAS given by Eric Ghysels.
%Description:
%   This example shows how to estimate short and long term components of the
%   conditional varaince using the GARCH-MIDAS apprach. We show examples
%   where long-term components are driven by RV and macro data, e.g. IP
%   growth rate and the unemployment rate.
%Information: 
%   Author - Anessa Custovic (acustovic12@gmail.com) | 21Aug2019
%   GarchMidas function is called from the MIDAS toolbox written by Hang Qian.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% housekeeping
clear all
close all

% add scripts to filepath
addpath(genpath('MIDASLecture/Example_Scripts'));
% add MIDAS toolbox to filepath
addpath(genpath('MIDASLecture/MIDASv2.2/MIDASv2.2/privatex'));
% add data to filepath
addpath(genpath('MIDASLecture/Example_Data'));


%% First Example is with default RV calculated in the GARCH-MIDAS function
load('NASDAQ_daily.mat');
% trim data so start at same date..
start_indxy=find(NASDAQ.DATE=='2019-06-28');
NASDAQ_trunc=NASDAQ(1:start_indxy,:);

y = NASDAQ_trunc.NASDAQCOM_PCH./100;

period = 22;
numLags = 24;


% Estimate the GARCH-MIDAS model, and extract the volatilities
% Notice we do not enter an X variable in the arguments of GarchMidas, so
% the default option is to calculate RV from the data

[estParams,EstParamCov,Variance,LongRunVar] = GarchMidas(y,'Period',period,'NumLags',numLags); 

% Plot the conditional volatility and its long-run component
figure(1)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.2,2019.7,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1973,2019])
ylim([0,1])
hold off

%% if we want to recover the weights, we can do the following:
param1 = estParams(5);
seq = numLags:-1:1; 
weights = (1-seq./numLags+10*eps).^(param1-1);    
% now re-scale weights...   
weights = weights ./ nansum(weights);


% Plot the weights of its long-run component
figure(2)
plot(seq,weights','LineWidth',1);
legend('Beta Weights','Location','NorthEast')
xlim([1,24])
title('Beta Weights Fixed Default RV')



% Estimate the rolling window version of the GARCH-MIDAS model where RV is
% the long term component
[estParams,EstParamCov,Variance,LongRunVar]...
    = GarchMidas(y,'Period',period,'NumLags',numLags,'RollWindow',1);

% Plot the conditional volatility and its long-run component
figure(3)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.2,2019.8,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1973,2019])
ylim([0,1])
hold off


% if we want to recover the weights, we can do the following:
param1 = estParams(5);
seq = numLags:-1:1; 
weights = (1-seq./numLags+10*eps).^(param1-1);    
%now re-scale weights...   
weights = weights ./ nansum(weights);

% Plot the weights of its long-run component
figure(4)
plot(seq,weights','LineWidth',1);
legend('Beta Weights','Location','NorthEast')
xlim([1,24])
title('Beta Weights Rolling Default RV')



%% Now we run GARCH-MIDAS where exogonous variable is IP growth. 
% Industrial Production Index growth rate, 1971-2015
% Data Source: FRED database 
% https://research.stlouisfed.org/fred2/series/INDPRO
load('INDPRO_monthly.mat');

% truncate the data
yDate=NASDAQ_trunc.DATE(1:end,:);
yDateDetails = datevec(yDate);
yDateMonth = yDateDetails(:,2); %grab all months..

% trim data so start at same date..
start_indx=find(IndPro.DATE=='1971-02-01');
IndPro_trunc=IndPro(start_indx:end,:);

xMonth = IndPro_trunc.INDPRO_PCH./ 100;


nobs = size(y,1);
xDay = NaN(nobs,1);
count = 1;

% now we repeat the monthly values to match the daily data
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)    
        count = count + 1;
        if count > length(xMonth)
            break
        end
    end
    xDay(t) = xMonth(count);
end


% Estimate the GARCH-MIDAS model with the exogenous regressor
% Note: here we put 'ThetaM',1 in the inputs-- this let's the GarchMidas
% function know not to take squares for the parameter theta and m in the 
% long-run volatility component. The default is false (they are squared)
[estParams,EstParamCov,Variance,LongRunVar] = GarchMidas(y,'Period',period,'NumLags',24,'X',xDay,'ThetaM',1);

% Plot the conditional volatility and its long-run component
figure(5)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.1,2019.7,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1973,2019])
ylim([0,1])
hold off

% if we want to recover the weights, we can do the following:
param1 = estParams(5);
seq = numLags:-1:1; 
weights = (1-seq./numLags+10*eps).^(param1-1);    
% now re-scale weights...   
weights = weights ./ nansum(weights);


% Plot the weights of its long-run component
figure(6)
plot(seq,weights','LineWidth',1);
legend('Beta Weights','Location','NorthEast')
xlim([1,24])
title('Beta Weights')





%% Now we run GARCH-MIDAS where exogonous variable is Consumer Sentiment. 

load('UMCSENT_quarterly.mat')

yDate=NASDAQ_trunc.DATE(1:end,:);
yDateDetails = datevec(yDate);
yDateMonth = yDateDetails(:,2); %grab all months..


% Now need to translate months into corresponding quarters....
yQuartMonth((yDateMonth==1),:)= 1;
yQuartMonth((yDateMonth==2),:)= 1;
yQuartMonth((yDateMonth==3),:)= 1;

yQuartMonth((yDateMonth==4),:)= 4;
yQuartMonth((yDateMonth==5),:)= 4;
yQuartMonth((yDateMonth==6),:)= 4;

yQuartMonth((yDateMonth==7),:)= 7;
yQuartMonth((yDateMonth==8),:)= 7;
yQuartMonth((yDateMonth==9),:)= 7;

yQuartMonth((yDateMonth==10),:)= 10;
yQuartMonth((yDateMonth==11),:)= 10;
yQuartMonth((yDateMonth==12),:)= 10;

xMonth = UMCSENT.UMCSENT./100;
nobs = size(y,1);

xDay = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yQuartMonth(t) ~= yQuartMonth(t-1)    
        count = count + 1;
        if count > length(xMonth)
            break
        end
    end
    xDay(t) = xMonth(count);
end

% Estimate the GARCH-MIDAS model with the exogenous regressor
period = 66;
numLags = 12;

% the estimated theta is negative, meaning that an increase
% in consumer sentiment is associated with a decline in long-term volatility.
[estParams,EstParamCov,Variance,LongRunVar,ShortRunVar,logL] = GarchMidas(y,'Period',period,'NumLags',numLags,'X',xDay,'ThetaM',1);

% Plot the volatilities
figure(7)
nobs = size(y,1);
seq = (period*numLags+1:nobs)';
year = linspace(1971.1,2019.7,nobs);
plot(year(seq),sqrt(252*Variance(seq)),'g--','LineWidth',1);
hold on
plot(year(seq),sqrt(252*LongRunVar(seq)),'b-','LineWidth',2);
legend('Total Volatility','Secular Volatility','Location','SouthEast')
xlim([1973,2019])
ylim([0,1])
hold off


% if we want to recover the weights, we can do the following:
param1 = estParams(5);
seq = numLags:-1:1; 
weights = (1-seq./numLags+10*eps).^(param1-1);    
% now re-scale weights...   
weights = weights ./ nansum(weights);


% Plot the weights of its long-run component
figure(8)
plot(seq,weights','LineWidth',1);
legend('Beta Weights','Location','NorthEast')
xlim([1,numLags])
title('Beta Weights')

