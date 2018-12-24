function [ FastHGP ] = Newton( FastHGP )

    tic
    x=FastHGP.UVonCones;
    
    % constants
    n=(length(x)+length(FastHGP.FixedIndices))/2;
    SPDHessian=true;
    maxIter=500;
    if(~FastHGP.ATPSucsess)%if ATP failed - don't run Newton
        maxIter=0;
    end
    stepTol=1e-7;
    gradTol=1e-7;
    Etol=1e-8;
    EitersCounter=0;
    numConvIters=5;

    % check energy initial values
    [fz, fbz ,E]=symDirEnergyGradHess(x,FastHGP.J_fz,FastHGP.J_fbz,FastHGP.Area,FastHGP.FixedIndices,FastHGP.FixedValues,SPDHessian);
    FastHGP.initEnergy=E;
    
    % optimize scale
    scale = optimizeSymDirEnergyByGlobalScaling( fz, fbz, FastHGP.Area );
    x=scale*x;
    FastHGP.FixedValues=scale*FastHGP.FixedValues;

    % move arrays to gpu
    J_fz_gpu=gpuArray(FastHGP.J_fz);
    J_fbz_gpu=gpuArray(FastHGP.J_fbz);
    Area_gpu=gpuArray(FastHGP.Area);
    FixedValues_gpu=gpuArray(FastHGP.FixedValues);
    
    % run Newton
    for iter=1:maxIter
        Eprev=E;

        % ****find energy,gradient, hessian****
        x_gpu=gpuArray(x);
        [ fz_gpu,fbz_gpu,E_gpu,G_gpu,H_gpu ] = symDirEnergyGradHess(x_gpu,J_fz_gpu,J_fbz_gpu,Area_gpu,FastHGP.FixedIndices,FixedValues_gpu,SPDHessian);

        % move to cpu
        E=gather(E_gpu);
        G=gather(G_gpu);
        H=gather(H_gpu);
        fz=gather(fz_gpu);
        fbz=gather(fbz_gpu);

        % ****solve****
        [R,p]=chol(H);
        assert(p==0);
        d=R\(R'\(-G));


        % ****line search****
        d_full=zeros(length(d)+length(FastHGP.FixedIndices),1);
        d_full(setdiff(1:end,FastHGP.FixedIndices))=d; % set free values 
        d_full(FastHGP.FixedIndices)=0; % set all fixed values to 0

        d_complex=complex(d_full(1:n),d_full(n+1:end));
        d_fz=FastHGP.J_fz*d_complex;
        d_fbz=FastHGP.J_fbz*d_complex;

        % locally injective line search
        tMax = lineSearchLocalInjectivity( fz,fbz,d_fz,d_fbz);

        % decreasing energy line search
        t = lineSearchDecreasingEnergy( fz,fbz,d_fz,d_fbz,FastHGP.Area,E,G,d,tMax );

        % ****do step****
        x=x+t*d;

        % stop conditions
        if norm(G)<gradTol
            break;
        end
        if (t*norm(d)<stepTol)
            break;
        end

        if (abs(E - Eprev) < Etol*(E + 1))
            if (EitersCounter >= numConvIters)
                break;
            else
                EitersCounter = EitersCounter + 1;
            end
        else
            EitersCounter = 0;
        end

    end

    FastHGP.newtonTime=toc;
    FastHGP.newtonIters=iter;
    FastHGP.NewtonFinalEnergy=E;

    % get full solution (with fixed)
    solution=zeros(length(x)+length(FastHGP.FixedIndices),1); 
    solution(setdiff(1:end,FastHGP.FixedIndices))=x; % set free values 
    solution(FastHGP.FixedIndices)=FastHGP.FixedValues; % set fixed values
    solution=[solution(1:n),solution(n+1:end)];   
    FastHGP.UVonCones=solution;

    FastHGP.NewtonSucsess=~ isempty(iter) && ~(iter==maxIter);
    FastHGP.NewtonSucsess=double(FastHGP.NewtonSucsess);

end