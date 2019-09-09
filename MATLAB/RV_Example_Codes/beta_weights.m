function [be_values]=beta_weights(dayLag,k1,k2)
%Beta weights function.
%Arthur Sinko, UNC, 2007
u=linspace(eps,1-eps,dayLag)';

k1=(k1*(k1>0)+1e-8*~(k1>0))*(k1<300)+300*~(k1<300);% 0<k1<300
k2=(k2*(k2>0)+1e-8*~(k2>0))*(k2<300)+300*~(k2<300);% 0<k1<300

be_values=u.^(k1-1).*(1-u).^(k2-1);%u.^(k1-1).*(1-u).^(k2-1)./beta(k1,k2);
be_values=be_values/sum(be_values);