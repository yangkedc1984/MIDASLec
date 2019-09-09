function fore = midas_filter(coeff,x,weights)
if strcmp(weights,'beta')==1
    smpl=size(x,2);
    
    be_values=beta_weights(smpl,coeff(3),coeff(4));
    fore=coeff(1)+coeff(2)*sum(x*be_values,2);
    
elseif strcmp(weights,'betaconstr')==1
    smpl=size(x,2);
    
    be_values=beta_weights(smpl,1,coeff(3));
    fore=coeff(1)+coeff(2)*sum(x*be_values,2);
    
elseif strcmp(weights,'exp')==1
    smpl=size(x,2);
    
    be_values=exp_weights(smpl,coeff(3),coeff(4));
    fore=coeff(1)+coeff(2)*sum(x*be_values,2);
    
elseif strcmp(weights,'Corsi')==1
    smpl=size(x,2);
    
    be_values=Corsi_weights(smpl,coeff(3),coeff(4),coeff(5));
    fore=coeff(1)+coeff(2)*sum(x*be_values,2);
end

end