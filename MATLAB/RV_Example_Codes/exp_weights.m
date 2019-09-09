function [poli i]=exp_weights(dayLag,k1,k2)

iii=transpose(1:dayLag);
poli=exp(k1*iii+k2*iii.^2)/sum(exp(k1*iii+k2*iii.^2));
if isnan(sum(poli))

    exppoly=k1*iii+k2*iii.^2;
    for i=dayLag:-1:1
        poli(i,1)=1/sum(exp(exppoly-k1*i-k2*i.^2));
    end
end

% subplot(2,1,1)
% plot(poli1)
% subplot(2,1,2)
% plot(poli)


% A=repmat(iii,1,dayLag);
% D=A-A';
% poli=sum(exp(k1*D-k2*D.^2),2);