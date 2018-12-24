function [ scale ] = optimizeSymDirEnergyByGlobalScaling( fz, fbz, Area )
% calculate the scale that minimizes the dirichlet energy.
% this simple trick helps to reduce unecessary iterations from the newton solver 

    a=abs(fz).^2;
    b=abs(fbz).^2;
       
    x=Area'*((a+b)./(a-b).^2);
    y=Area'*(a+b);
    s=x/y;
    
    scale=s^(1/4);

%     fprintf('scale=%f\n',scale);
        
end

