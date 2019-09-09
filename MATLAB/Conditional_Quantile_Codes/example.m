clear;clc;

load y
[beta,condQ] = EstimateConditionalQuantile('C', 1, 0.05, y, empiricalQuantile, [], []);
plot(condQ)
hold on 
plot(y)
