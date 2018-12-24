function [E] = symDirEnergyByfzfbz(fz,fbz,Area)

    % ----------------------------
    % calcultae Energy
    % ----------------------------

    fz2=abs(fz).^2;
    fbz2=abs(fbz).^2;
    
    g1=fz2+fbz2;
    g2=fz2-fbz2;

    E=sum(Area.*(g1+g1./g2.^2));

end