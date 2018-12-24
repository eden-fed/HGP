function [ fz_local , fbz_local ] = localStep( frames , sigma2Eps , k , fz , fbz )
    y=real(fz.*frames);
    x=abs(fbz);
    assert(all(x>=0));
    
    %project jacobians to convex space
%     I1 = (y>=x./k) & (y>=x+sigma2Eps);
    I2 = (y<x./k) & (y>-k*x+sigma2Eps*((1+k^2)/(1-k)));
    I3 = (y<=-k*x+sigma2Eps*((1+k^2)/(1-k))) & (y>=-x+sigma2Eps*(1+k)/(1-k));
    I4 = (y<-x+sigma2Eps*(1+k)/(1-k)) & (y>-x+sigma2Eps) & (y<sigma2Eps+x);
    I5 = (y<=-x+sigma2Eps);
    
%     assert(all((~I2 & ~I3 & ~I4 & ~I5)==I1))
    
    y1 = y; x1 = x;
    
    y1(I2) = (y(I2)+k*x(I2))./(1+k^2);
    x1(I2) = (k*y(I2)+(k)^2*x(I2))./(1+k^2);
    
    y1(I3) = sigma2Eps/(1-k);
    x1(I3) = k*sigma2Eps/(1-k);
    
    y1(I4) = (y(I4)+x(I4)+sigma2Eps)./2;
    x1(I4) = (y(I4)+x(I4)-sigma2Eps)./2;

    y1(I5) = sigma2Eps;
    x1(I5) = 0;
    
    
    %return value
    fbz_local=x1.*fbz./abs(fbz);
    fbz_local(abs(fbz)==0)=0;
    fz_local=complex(y1,imag(fz.*frames))./frames;
    
    tol=1e-4;
    assert(all(real(fz_local.*frames) >= abs(fbz_local)+sigma2Eps-tol))
    assert(all(k*real(fz_local.*frames) >= abs(fbz_local)-tol))
end

