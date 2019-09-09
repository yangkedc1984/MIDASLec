function [VaR] = EvaluateVaR(flag, model, beta, theta, r, empiricalQuantile, ...
    rHF, midasOptions)

lVaR = length(r);
VaR = zeros(lVaR, 1);

if isempty(rHF)
    % CAViaR
    VaR(1) = empiricalQuantile;
    switch model 
        % SAV
        case 1
            for i = 2:lVaR
                VaR(i) = beta(1) + beta(2) * VaR(i - 1) + ...
                    beta(3) * abs(r(i - 1));
            end
                
        % AS
        case 2
            for i = 2:lVaR
                VaR(i) = beta(1) + beta(2) * VaR(i - 1) + ...
                    beta(3) * (r(i - 1) > 0) - beta(4) * (r(i - 1) < 0);
            end
                
        % GARCH
        case 3
            for i = 2:lVaR
                VaR(i) = sqrt(beta(1) + beta(2) * (VaR(i - 1)^2) + ...
                    beta(3) * (r(i - 1)^2));
            end                
            
        % Adaptive 
        case 4
            for i = 2:lVaR
                VaR(i) = VaR(i - 1) + beta(1) * (1/(1 + exp(beta(2) * ...
                    (r(i - 1) - VaR(i - 1)))) - theta);
            end
    end
    
else
    [lag, type] = deal(midasOptions.lag, midasOptions.type);
    weights = {};
    rHFPos = rHF .* (rHF > 0);
    rHFNeg = - rHF .* (rHF < 0);
    
    switch lower(flag)
        case 'm'
        % MIDAS
        if model == 2
            [weights.pos, ~] = CreateWeightingVector(lag, [1; beta(4)], type);
            [weights.neg, ~] = CreateWeightingVector(lag, [1; beta(5)], type);
        else
            [weights, ~] = CreateWeightingVector(lag, [1; beta(3)], type);
        end
        
        switch model
            % SAV
            case 1
                VaR(1) = beta(1) + beta(2) * (weights(1:lag - 1)' * ...
                    abs(rHF(lag - 1:-1:1)));
                for i = 2:lVaR
                    VaR(i) = beta(1) + beta(2) * (weights' * ...
                        abs(rHF((lag * i - 1):-1:lag * (i - 1))));
                end

            % AS    
            case 2
                VaR(1) = beta(1) + beta(2) * (weights.pos(1:lag - 1)' * ...
                        rHFPos(lag - 1:-1:1)) + ...
                        beta(3) * (weights.neg(1:lag - 1)' * ...
                        rHFNeg(lag - 1:-1:1));
                for i = 2:lVaR
                    VaR(i) = beta(1) + beta(2) * (weights.pos' * ...
                        rHFPos((lag * i - 1):-1:lag * (i - 1))) + ...
                        beta(3) * (weights.neg' * ...
                        rHFNeg((lag * i - 1):-1:lag * (i - 1)));
                end
            
            % GARCH    
            case 3
                VaR(1) = sqrt(beta(1) + beta(2) * (weights(1:lag - 1)' * ...
                        (rHF(lag - 1:-1:1)).^2));
                for i = 2:lVaR
                    VaR(i) = sqrt(beta(1) + beta(2) * (weights' * ...
                        (rHF((lag * i - 1):-1:lag * (i - 1))).^2));
                end
            
            % Adaptive
            case 4
                VaR(1) = beta(1) + beta(2) * (1/(1 + exp(beta(4) * ...
                    (weights(1:lag - 1)' * rHF(lag - 1:-1:1) - ...
                    beta(1)))) - theta);
                for i = 2:lVaR
                    VaR(i) = beta(1) + beta(2) * (1/(1 + exp(beta(4) * ...
                        (weights' * rHF((lag * i - 1):-1:lag * (i - 1)) - ...
                        beta(1)))) - theta);
                end        
        end
        
        case 'h'
        % HYBRID
        if model == 2
            [weights.pos, ~] = CreateWeightingVector(lag, [1; beta(5)], type);
            [weights.neg, ~] = CreateWeightingVector(lag, [1; beta(6)], type);
        elseif model == 4
            [weights, ~] = CreateWeightingVector(lag, [1; beta(2)], type);
        else
            [weights, ~] = CreateWeightingVector(lag, [1; beta(4)], type);
        end
        
        switch model
            % SAV
            case 1
                VaR(1) = beta(1) + beta(2) * empiricalQuantile + ...
                    beta(3) * (weights(1:lag - 1)' * ...
                    abs(rHF((lag - 1):-1:1)));
                for i = 2:lVaR
                    VaR(i) = beta(1) + beta(2) * VaR(i - 1) + ...
                        beta(3) * (weights' * ...
                        abs(rHF((lag * i - 1):-1:lag * (i - 1))));
                end
            
            % AS    
            case 2
                VaR(1) = beta(1) + beta(2) * empiricalQuantile + ...
                    beta(3) * (weights.pos(1:lag - 1)' * ...
                    rHFPos((lag - 1):-1:1)) + beta(4) * ...
                    (weights.neg(1:lag - 1)' * rHFNeg((lag - 1):-1:1));
                for i = 2:lVaR
                    VaR(i) = beta(1) + beta(2) * VaR(i - 1) + ...
                        beta(3) * (weights.pos' * ...
                        rHFPos((lag * i - 1):-1:lag * (i - 1))) + ...
                        beta(4) * (weights.neg' * ...
                        rHFNeg((lag * i - 1):-1:lag * (i - 1)));
                end
                
            % GARCH    
            case 3
                VaR(1) = sqrt(beta(1) + beta(2) * (empiricalQuantile^2) + ...
                    beta(3) * (weights(1:lag - 1)' * ...
                    (rHF((lag - 1):-1:1)).^2));
                for i = 2:lVaR
                    VaR(i) = sqrt(beta(1) + beta(2) * (VaR(i - 1)^2) + ...
                        beta(3) * (weights' * ...
                        (rHF((lag * i - 1):-1:lag * (i - 1))).^2));
                end 
             
            % Adaptive    
            case 4
                VaR(1) = empiricalQuantile + beta(1) * ...
                    (1/(1 + exp(beta(3) * (weights(1:lag - 1)' * ...
                    rHF((lag - 1):-1:1) - empiricalQuantile))) - theta);
                for i = 2:lVaR
                    VaR(i) = VaR(i - 1) + beta(1) * (1/(1 + exp(beta(3) * ...
                        (weights' * rHF((lag * i - 1):-1:lag * (i - 1)) - ...
                        VaR(i - 1)))) - theta);
                end 
        end
    end
end
end