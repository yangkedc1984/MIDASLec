%% CORE LECTURE SERIES ON MIDAS 
% EXAMPLE 2: VOLATILITY FORECASTING 
%Purpose:
%   To follow along with the volatility forecasting slides for the Core Lecture
%   Series on MIDAS given by Eric Ghysels.
%Description:
%   The following script recreates part of the output from the in-sample 
%   estimation of "Volatility Forecasting Across Asset Classes: Multi-Period 
%   Forecasts" Annual Review of Financial Economics
%   E. Ghysels, A. Plazzi, R. Valkanov, A. Rubia, and A. Dossani
%Information: 
%   All code and data called below was graciously provided by the authors
%   of the paper.
%   Adapted by - Anessa Custovic (acustovic12@gmail.com) | 21Aug2019
%   for the CORE Lecture Series on MIDAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; 
close all;

% add scripts to filepath
addpath(genpath('MIDASLecture/Example_Scripts'));
% add MIDAS toolbox to filepath
addpath(genpath('MIDASLecture/MIDASv2.2/MIDASv2.2/privatex'));
% add data to filepath
addpath(genpath('MIDASLecture/Example_Data'));
% add RV example codes to filepath
addpath(genpath('MIDASLecture/RV_Example_Codes'));


%load in data
load('example2.mat')


% General options
optionsmidas.aggrX=126; % no. of lags of RV
Horizons = [5 10]; % forecasting horizons


%%% loop across assets
for countf = 1 : 5 
    
    %load in data
    load('example2.mat')

    Rd =   AllclosingRet_2018{countf,1};
    dates = Alldates_2018{countf,1};
    RV = AllRV_2018{countf,1};
    T = size(dates,1); % no of time-series obs
    
    
    counth = 1;
    %%%     loop across horizons
    for hhh = Horizons
        
        disp([countf hhh])
        optionsmidas.aggrY=hhh;
        
        [RVh,~,dates_hhh] = datageneration_midas(RV,optionsmidas,dates);
        
        %%%     MIDAS on RV with beta weights
        [theta_MIDAS_beta, ~,s2_MIDAS_beta]= midas_optimization(RV,optionsmidas,'beta');
        %%%     MIDAS on RV with exp weights
        [theta_MIDAS_exp, objfn_MIDAS_exp,s2_MIDAS_exp]= midas_optimization(RV,optionsmidas,'exp');
                
        % collect forecasts --- 
        %[first col is beta weight forecasts second is expalmon forecasts]
        Forecasts = [s2_MIDAS_beta s2_MIDAS_exp];

        % evaluate forecasts performance
        for fff = 1 : size(Forecasts,2)
            errors2(:,fff) = (RVh - Forecasts(:,fff)).^2;
            qlike(:,fff) = log(Forecasts(:,fff))+RVh./Forecasts(:,fff);
            RMSFE(counth,fff,countf) = sqrt(mean(errors2(:,fff)));
            QLIKE(counth,fff,countf) = mean(qlike(:,fff));
        end
        
        Dates{countf,counth}=dates_hhh; % where it is stored the first date of the forecasted, non-overlapping period
        Theta_MIDAS_exp(:,counth,countf) = theta_MIDAS_exp;
        Theta_MIDAS_beta(:,counth,countf) = theta_MIDAS_beta;
        MIDAS_forecasts{countf,counth}=Forecasts;

        clear s2* theta_* RVh dates_hhh hhh_* errors2* qlike* obj* dates_*
        counth = counth+1;
    end
    
    % start over with a clean sheet/
    clearvars -except optionsmidas Horizons RMSFE QLIKE countf Theta* Dates agg AggrX MIDAS_forecasts
end


% based off of QLIKE which model had the lowest value
% the lower the value of QLIKE the better the forecasts

% QLIKE column 1 is beta model column 2 is expalmon model
% QLIKE row 1 is horz 5, row 2 is horz 10
% each cell (:,:,x) x represents asset number
% asset 1: almon model has lower QLIKE for both horizons
% asset 2: almon model has lower QLIKE for both horizons
% asset 3: almon model better for horz 5, beta for 10
% asset 4: almon has lower QLIKE (more negative)
% asset 5: almon has lower QLIKE for all horizons here

% RMSFE
% asset 1: almon model has slightly higher RMSFE for horz10
% asset 2: almon model has lower RMSFE for both horizons
% asset 3: almon model has lower RMSFE for both horizons
% asset 4: looks very similar
% asset 5: almon has lower RMSFE both horizons

