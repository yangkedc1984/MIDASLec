%% CORE LECTURE SERIES ON MIDAS 
% EXAMPLE 4: MIDAS QUANTILE REGRESSION
%Purpose:
%   To follow along with the quantile regression slides for the Core Lecture
%   Series on MIDAS given by Eric Ghysels.
%Description:
%   The following script goes through a number of MIDAS quantile regression
%   examples. It compared MIDAS results to unconditional quantiles from a
%   rolling window. It also calculates CAViaR quantiles to compare to MIDAS
%   results.
%Information: 
%   Author - Anessa Custovic (acustovic12@gmail.com) | 21Aug2019
%   All code written in the Conditional_Quantile_Codes folder was written
%   by Hanwei Liu.
%   MidasQuantile_edited function was modified by Anessa Custovic from the 
%   MIDAS toolbox written by Hang Qian for this example.
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
% add conditional quantile folder to filepath
addpath(genpath('MIDASLecture/Conditional_Quantile_Codes'));


%load S&P500 Data
load('example4.mat');

%computer log returns
returns= log(snp500(2:end,2)./snp500(1:end-1,2));
%set dates
dates=snp500(2:end,1);
dates = datenum(dates);


y=returns;
% Fit 25% conditional quantile
[estParams,condQuantile3,yLowFreq,xHighFreq,yDates] = MidasQuantile_edited(y,'Dates',dates,'Period',5,'NumLags',22,'Quantile',0.25);

% Plot the 25% conditional quantile
plot(yDates,yLowFreq)
hold on;
plot(yDates,condQuantile3)
legend('low freq','cond quant')
title('MIDAS Quantile q=0.25')
dateaxis;
hold off;


% Fit 5% conditional quantile
[estParams2,condQuantile32,yLowFreq2,xHighFreq2,yDates2] = MidasQuantile_edited(y,'Dates',dates,'Period',5,'NumLags',22,'Quantile',0.05);

% Plot the 5% conditional quantile
plot(yDates2,yLowFreq2)
hold on;
plot(yDates2,condQuantile32)
dateaxis;
legend('low freq','cond quant')
title('MIDAS Quantile q=0.05')
hold off;


% now rolling window 5% unconditional quantile estimate for comparison
y = yLowFreq2;
windowsize = 100;

unquant = NaN(length(y),1);

for j = windowsize:length(y)
  unquant(j,:) = quantile(y(((j-windowsize+1):j),:), 0.05);
end



% Plot the 5% conditional quantile and 5% unconditional quantile
plot(yDates2,yLowFreq2,'k-')
hold on;
plot(yDates2,condQuantile32,'r-')
plot(yDates2,unquant,'b-')
dateaxis;
legend('low freq','cond quant','unc quant')
title('Quantile q=0.05')
hold off;





%now caviar estiamtion
empiricalQuantile =unquant(101); %to initialize loop
date=yDates;

[beta,condQ] = EstimateConditionalQuantile('C', 1, 0.05, y, empiricalQuantile, [], []);


% Plot the unconditional quantile and MIDAS
% and plot caviar results with unconditional
subplot(2,1,1)
plot(yDates2,yLowFreq2,'k-')
hold on;
plot(yDates2,condQuantile32,'r-')
plot(yDates2,unquant,'b-')
dateaxis;
legend('low freq','cond quant','unc quant')
title('MIDAS Quantile q=0.05')
hold off;

subplot(2,1,2)
plot(yDates2,yLowFreq2,'k-')
hold on;
plot(yDates2,condQ,'r-')
plot(yDates2,unquant,'b-')
dateaxis;
legend('low freq','cond quant','unc quant')
title('CaViaR Quantile q=0.05')
hold off;



% Plot quantiles for MIDAS and CAViaR to compare

plot(yDates2,condQuantile32,'k-')
hold on;
plot(yDates2,condQ,'r-')
dateaxis;
legend('MIDAS','CaViaR')
title('MIDAS vs CaViaR Quantile q=0.05')
hold off;




% compute conditional skewness at 0.95 level for MIDAS and CAViaR, plot
% first compute CAViaR for 0.5 and 0.95 levels
[beta_5,condQ_5] = EstimateConditionalQuantile('C', 1, 0.5, y, quantile(y(1:100),0.5), [], []);

[beta_95,condQ_95] = EstimateConditionalQuantile('C', 1, 0.95, y, quantile(y(1:100),0.95), [], []);

% compute conditional skewness
[cskew_CaViaR] = condskewness(condQ_95,condQ,condQ_5,0.95);


% MIDAS quantiles
% this one is tricky so we use search and set initial parameters to avoid
% falling into local mins
[estParams_50,condQuantile_50,yLowFreq_50,xHighFreq_50,yDates_50] = MidasQuantile_edited(returns,'Dates',dates,'Period',5,'NumLags',22,'Quantile',0.50,'Search',1,'Params0',[0.0009;.2;4]);

[estParams_95,condQuantile_95,yLowFreq_95,xHighFreq_95,yDates_95] = MidasQuantile_edited(returns,'Dates',dates,'Period',5,'NumLags',22,'Quantile',0.95);


% skew MIDAS
[cskew_MIDAS] = condskewness(condQuantile_95,condQuantile32,condQuantile_50,0.95);



% plot conditional skewness 
% throw out first 10% of sample due to initialization
idx=(round(length(cskew_CaViaR)*0.1)+1):length(cskew_CaViaR);

plot(yDates2(idx),cskew_MIDAS(idx),'k-')
hold on;
plot(yDates2(idx),cskew_CaViaR(idx),'r-')
dateaxis;
legend('MIDAS','CaViaR')
title('95% Conditional Skewness Est.')
hold off;







