function  [coeff_MIDAS, objfn_MIDAS,cond_var_MIDAS, dates_MIDAS]= midas_optimization(rv,optionsmidas,weights)
% main function for MIDAS optimization

N = 1;  % total # of regressors

%______________________ Data construction _____________________________%
%___________________________________________________________________%
[y,X] = datageneration_midas(rv,optionsmidas);

%______________________ Initial guess _____________________________%
%______________________________________________________________%
% Select best starting points: select numInVec random initial values, calculate the obj fn at these points,
% then select the best nInitialCond ones as candidate initial values from where to start optimizing

rng(1);      % initialize the seed for replicability
Ns=1000;    % number of trial vectors
nInitialCond = 3;             % Number of initial conditions for the optimization

if strcmp(weights,'beta') == 1
    % random vectors of MIDAS coeff
    initialTargetVectors=NaN(Ns,N*3);
    for col = 1:N % for each regressor, if more than 1
        initialTargetVectors(:,(col-1)*3+1:col*3) = [randn(Ns,1) unifrnd(0, 40,[Ns 2])];
    end
elseif strcmp(weights,'betaconstr') == 1
    % random vectors of MIDAS coeff
    initialTargetVectors=NaN(Ns,N*2);
    for col = 1:N % for each regressor, if more than 1
        initialTargetVectors(:,(col-1)*2+1:col*2) = [randn(Ns,1) unifrnd(0, 40,[Ns 1])];
    end
elseif strcmp(weights,'exp') == 1
    % random vectors of MIDAS coeff
    initialTargetVectors=NaN(Ns,N*3);
    for col = 1:N % for each regressor, if more than 1
        initialTargetVectors(:,(col-1)*3+1:col*3) = [randn(Ns,1) unifrnd(-0.1, 0,[Ns 1]) unifrnd(-0.1, 0,[Ns 1])];
    end
end
initialTargetVectors=[rand(Ns,1)/10 initialTargetVectors]; % add the constant term

% value of obj fn at these points
V = NaN(Ns, 1);
for k = 1:Ns
    V(k,1) = midas_objfn(initialTargetVectors(k,:)',y,X,1,weights);
end
Results          = [V, initialTargetVectors];
SortedResults    = sortrows(Results,1); % sort the results
BestInitialCond  = SortedResults(1:nInitialCond,2:end); % keep the top nInitialCond

% add point from genetic algorithm as initial condition
optionsga = optimset('Display', 'none');

if strcmp(weights,'beta') == 1
    lb=[-10e+4; -10e+2; 1+eps;1+eps];
    ub=[10e+4; 10e+2; 10e+2; 10e+2];
    a_ga=ga(@(a) midas_objfn(a,y,X,1,weights),1+3*N,[],[],[],[],lb,ub,[],optionsga);
elseif strcmp(weights,'betaconstr') == 1
    lb=[-10e+4; -10e+2; 1+eps];
    ub=[10e+4; 10e+2; 10e+2];
    a_ga=ga(@(a) midas_objfn(a,y,X,1,weights),1+2*N,[],[],[],[],lb,ub,[],optionsga);
elseif strcmp(weights,'exp') == 1
    bth = .75;
    lb = [ -10e+4; -10e+2; -bth; -bth];
    ub = [ 10e+4; 10e+2; bth; bth];
    a_ga=ga(@(a) midas_objfn(a,y,X,1,weights),1+3*N,[],[],[],[],lb,ub,[],optionsga);
end
BestInitialCond=[BestInitialCond;a_ga];

%______________________ Actual estimation _____________________________%
%___________________________________________________________________%
MaxFunEvals = 1000;
MaxIter     = 1000;
options = optimset('LargeScale', 'off', 'HessUpdate', 'dfp','MaxFunEvals', MaxFunEvals, ...
    'Display', 'none', 'MaxIter', MaxIter, 'TolFun', 1e-10, 'TolX', 1e-10);

a = NaN(size(BestInitialCond));
fval = NaN(size(BestInitialCond,1),1);

for i = 1:size(BestInitialCond,1)
    a(i,:) = fmincon(@(v) midas_objfn(v,y,X,1,weights),BestInitialCond(i,:),[],[],[],[],lb,ub,[],options);
    [a(i,:), fval(i,1)] = patternsearch(@(v) midas_objfn(v,y,X,1,weights),a(i,:),[],[],[],[],lb,ub,[],options);
end

SortedFval  = sortrows([fval, a, BestInitialCond], 1);
coeff_MIDAS    = SortedFval(1, 2:size(a,2)+1)'; % the max of the max
objfn_MIDAS    = SortedFval(1);
cond_var_MIDAS = midas_objfn(coeff_MIDAS,y,X,2,weights);

if SortedFval(1,2)<0 % negative intercept, may lead to negative Forecast
    lb(1,1)=eps;
    [a, fval] = patternsearch(@(v) midas_objfn(v,y,X,1,weights),SortedFval(1, 2:size(a,2)+1),[],[],[],[],lb,ub,[],options);
    coeff_MIDAS    = a;
    objfn_MIDAS    = fval;
    cond_var_MIDAS = midas_objfn(coeff_MIDAS,y,X,2,weights);
end

end