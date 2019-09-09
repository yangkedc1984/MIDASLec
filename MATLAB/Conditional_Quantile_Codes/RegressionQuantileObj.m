function [output] = RegressionQuantileObj(flag, model, beta, theta, r, ...
    empiricalQuantile, rHF, midasOptions, out)

switch lower(flag)
    case 'c'
        VaR = EvaluateVaR('C', model, beta, theta, r, ...
            empiricalQuantile, [], []);
    case 'm'
        VaR = EvaluateVaR('M', model, beta, theta, r, [], rHF, midasOptions);
    case 'h'
        VaR = EvaluateVaR('H', model, beta, theta, r, ...
            empiricalQuantile, rHF, midasOptions);
end

Hit = (r < VaR) - theta;

% Compute the Regression Quantile criterion.
RQ  = -Hit' * (r - VaR);
if (RQ == Inf || RQ ~= RQ || ~isreal(RQ))
    RQ = 1e+100;
end

if strcmp(out, 'RQ')
    output = RQ;
else
    output = [VaR Hit];
end

end