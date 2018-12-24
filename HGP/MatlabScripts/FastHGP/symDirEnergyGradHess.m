function [ fz,fbz,E,G,H ] = symDirEnergyGradHess(x,J_fz,J_fbz,Area,FixedIndices,FixedValues,SPDHessian)

if nargin>4   
    x_small=x;
    if isa(x, 'gpuArray')
        x=zeros(length(x)+length(FixedIndices),1,'gpuArray');
    else
        x=zeros(length(x)+length(FixedIndices),1);
    end
      
    x(setdiff(1:end,FixedIndices))=x_small; % set free values 
    x(FixedIndices)=FixedValues; % set all fixed values to 0
end
if nargin<6
    SPDHessian=false;
end
    n=length(x)/2;
    
    x_complex=complex(x(1:n),x(n+1:end));
    fz=J_fz*x_complex;
    fbz=J_fbz*x_complex;

    abs_fz=abs(fz);
    abs_fbz=abs(fbz);

    SIGMA=abs_fz+abs_fbz;
    sigma=abs_fz-abs_fbz;
        
    % ----------------------------
    % calcultae Energy
    % ----------------------------
    E = calcEnergy;
    function E = calcEnergy
        E=0.5*Area'*(SIGMA.^2+SIGMA.^(-2)+sigma.^2+sigma.^(-2));
    end
  
    if nargout > 3 % gradient required
    % ----------------------------
    % calculate gradient of Energy
    % ----------------------------
        [G_complex, grad_h_SIGMA, grad_h_sigma, grad_SIGMA_x, grad_sigma_x, grad_abs_fz, grad_abs_fbz] = calcGradient;
        
        if nargin>4
            G=[real(G_complex(1:end)) ; imag(G_complex(1:end))]; 
            G=G(setdiff(1:end,FixedIndices)); % save only free values 
        else
            G=[real(G_complex) ; imag(G_complex)];
        end        
    end
    
    function [G_complex, grad_h_SIGMA, grad_h_sigma, grad_SIGMA_x, grad_sigma_x, grad_abs_fz, grad_abs_fbz] = calcGradient
        grad_h_SIGMA = Area.*((SIGMA.^4 - 1) ./ SIGMA.^3); % dim:|F| X 1
        grad_h_sigma = Area.*((sigma.^4 - 1) ./ sigma.^3); % dim:|F| X 1

        grad_abs_fz=bsxfun(@times,conj(J_fz),fz./abs_fz); % A'*A*x/norm(A*x) (where A*x=fz) for each triangle
        grad_abs_fbz=bsxfun(@times,conj(J_fbz),fbz./abs_fbz);% B'*B*x/norm(B*x) (where B*x=fbz) for each triangle

        grad_SIGMA_x=grad_abs_fz+grad_abs_fbz; % dim:|F| X n complex
        grad_sigma_x=grad_abs_fz-grad_abs_fbz; % dim:|F| X n complex

        G_complex=grad_SIGMA_x.'*grad_h_SIGMA+grad_sigma_x.'*grad_h_sigma;
    end
   
    if nargout > 4 % Hessian required
    % ----------------------------
    % calculate Hessian of Energy
    % ----------------------------        
        H=calcHessian;
        if nargin>4
            H=H(setdiff(1:end,FixedIndices),setdiff(1:end,FixedIndices)); % save only free values 
        end    
    end 
    function H = calcHessian
        grad_SIGMA_x_real=[real(grad_SIGMA_x) imag(grad_SIGMA_x)];
        grad_sigma_x_real=[real(grad_sigma_x) imag(grad_sigma_x)];
        
        grad_abs_fz_real=[real(grad_abs_fz) imag(grad_abs_fz)];
        grad_abs_fbz_real=[real(grad_abs_fbz) imag(grad_abs_fbz)];
        
        %%% hess energy %%%
        %hi is convex with respect to SIGMA,sigma > 0, so no need to fix this part
        %upper right and lower left elements of the hessian of h in each triangle are zero
        hess_h11 = Area.*(1+3*SIGMA.^(-4));
        hess_h22 = Area.*(1+3*sigma.^(-4));

        D_g1g1=grad_SIGMA_x_real'*bsxfun(@times,grad_SIGMA_x_real,hess_h11);
        D_g2g2=grad_sigma_x_real'*bsxfun(@times,grad_sigma_x_real,hess_h22);

        H1=D_g1g1+D_g2g2;
        
        %%%%%%
        
        Ar1=[real(J_fz) -imag(J_fz)]; % row 1 of A matrix in each triangle
        Ar2=[imag(J_fz) real(J_fz)]; % row 2 of A matrix in each triangle
        Br1=[real(J_fbz) -imag(J_fbz)]; % row 1 of B matrix in each triangle
        Br2=[imag(J_fbz) real(J_fbz)]; % row 2 of B matrix in each triangle
        
        %%%% term gathering %%%%
        alpha=grad_h_SIGMA+grad_h_sigma;
        beta=grad_h_SIGMA-grad_h_sigma;
        if SPDHessian
            alpha=max(alpha,0);
        end
        
        % hess first term
        scalarForMult=alpha./abs_fz;
        HA=Ar1'*bsxfun(@times,Ar1,scalarForMult) + Ar2'*bsxfun(@times,Ar2,scalarForMult); 
        HA=HA-grad_abs_fz_real'*bsxfun(@times,grad_abs_fz_real,scalarForMult);
        
        % hess second term
        scalarForMult=beta./abs_fbz;
        HB=Br1'*bsxfun(@times,Br1,scalarForMult) + Br2'*bsxfun(@times,Br2,scalarForMult); 
        HB=HB-grad_abs_fbz_real'*bsxfun(@times,grad_abs_fbz_real,scalarForMult);
               
        H2=HA+HB;
                
        H=H1+H2;
    end
end