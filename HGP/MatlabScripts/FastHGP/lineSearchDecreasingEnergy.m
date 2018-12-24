function [ t ] = lineSearchDecreasingEnergy( fz,fbz,d_fz,d_fbz,Area,E,G,d,tMax )
    t=tMax;
    alpha=1e-4;
    beta=0.5;
    alpha_G_d=alpha*G'*d;
    
    for ii=1:30
                
        fz_t=fz + t*d_fz;
        fbz_t=fbz + t*d_fbz;
        
        E_t=symDirEnergyByfzfbz(fz_t,fbz_t,Area);

        if(E_t < E+t*alpha_G_d)
            break;
        end
        
        t=beta*t;
    end
    
end

