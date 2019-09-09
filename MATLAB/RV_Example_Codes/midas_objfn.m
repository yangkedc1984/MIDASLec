function output = midas_objfn(theta,y,X,OUT,weights)
% objective fn for the MIDAS estimation

% MIDAS vol
MidasVar=theta(1)+theta(2)*midas_filter(theta(3:end),X,weights);

% Obj Fn.
Q  = (y - MidasVar)'*(y - MidasVar)*1000;

if Q == Inf || (Q ~= Q) || ~isreal(Q)
    Q = 1e+10000;
end

% Select the output of the program.
if OUT == 1
    output = Q;
elseif OUT ==2
    output = MidasVar;
else
    error('Wrong output selected. Choose OUT = 1 for Q, or OUT = 2 for MidasVol')
end

% function to create the MIDAS polynominal
    function y_h = midas_filter(a,x,weights)
        if strcmp(weights,'beta')==1
            smpl=size(x,2);
            
            be_values=beta_weights(smpl,a(1),a(2));
            y_h=sum(x*be_values,2);
            
        elseif strcmp(weights,'betaconstr')==1
            smpl=size(x,2);
            
            be_values=beta_weights(smpl,1,a(1));
            y_h=sum(x*be_values,2);
                        
        elseif strcmp(weights,'exp')==1
            
            w=exp_weights(size(x,2),a(1),a(2));
            y_h=sum(x*w,2);
            
        elseif strcmp(weights,'Corsi')==1
            
            w=Corsi_weights(size(x,2),a(1),a(2),a(3));
            y_h=sum(x*w,2);
            
        end
        
    end


end
