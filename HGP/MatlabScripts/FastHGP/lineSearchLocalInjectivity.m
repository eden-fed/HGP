function [ t ] = lineSearchLocalInjectivity( fz,fbz,d_fz,d_fbz )

    tol=1e-10;
    A=abs(d_fz).^2-abs(d_fbz).^2;
    B=2*(real(d_fz).*real(fz) + imag(d_fz).*imag(fz) - real(d_fbz).*real(fbz) - imag(d_fbz).*imag(fbz));
    C=abs(fz).^2-abs(fbz).^2;
    
    Delta = B.*B - 4.0*A.*C;
	
    Q=zeros(size(B));
    
    ind=B > 0;
    Q(ind)=-0.5 * (B(ind) + sqrt(Delta(ind)));
    ind=B<=0;
    Q(ind)=-0.5 * (B(ind) - sqrt(Delta(ind)));
    
	t1 = Q ./ A;
	t2 = C ./ Q;
    
    % linear equation
    ind = (abs(A)<=tol);
    t1(ind) = -C(ind)./B(ind);
    t2(ind) = t1(ind);
    
    % ignore complex
    ind = (Delta<=-tol);
    t1(ind)=Inf;
    t2(ind)=Inf;
     
    t1 = real(t1);
    t2 = real(t2);
          
    % get rid of negatives/imaginaries
    t1(t1<0) = inf;
    t2(t2<0) = inf;

    t = min(min(t1),min(t2));% pick the smallest one between t1,t2 and then the smallest one among all triangles
    t=min(1,t*0.9); % 0.9*t and not t because t is the limit, after thet thr triange is flipped

end

