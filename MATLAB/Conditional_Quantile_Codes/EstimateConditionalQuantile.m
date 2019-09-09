function [betaH, condQ, Hit, RQ, initCond] = ...
    EstimateConditionalQuantile(flag, model, theta, r, ...
    empiricalQuantile, rHF, midasOptions)

% Inputs
% flag: specify CAViaR/MIDAS/HYBRID
% model: specify model forms 1 - 4
%        SAV/AS/Indirect GARCH/Adaptive (please see documentation)
% r: low frequency returns
% empiricalQuantile: empirical quantiles from the burn-in period
%                    initial values for the autoregressive quantile term
% rHF: high frequency returns
% midasOptions: midas lags and weighting polynomial

switch lower(flag)
    case 'c'
        % CAViaR
        [betaH, RQ, initCond] = CondQuantileOptim('C', model, theta, r, ...
            empiricalQuantile, [], []);
        VaRHit = RegressionQuantileObj('C', model, betaH, theta, r, ...
            empiricalQuantile, [], [], 'VH');
    case 'm'
        % MIDAS
        [betaH, RQ, initCond] = CondQuantileOptim('M', model, theta, r, [], ...
            rHF, midasOptions);
        VaRHit = RegressionQuantileObj('M', model, betaH, theta, r, [], ...
            rHF, midasOptions, 'VH');
    case 'h'
        % HYBRID
        [betaH, RQ, initCond] = CondQuantileOptim('H', model, theta, r, ...
            empiricalQuantile, rHF, midasOptions);
        VaRHit = RegressionQuantileObj('H', model, betaH, theta, r, ...
            empiricalQuantile, rHF, midasOptions, 'VH');
end

[condQ, Hit] = deal(VaRHit(:, 1), VaRHit(:, 2));
end