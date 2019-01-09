function [ FastHGP ] = ATPForInitialValue( FastHGP )
    tic
    numCones=size(FastHGP.J_fz,2);
    numTriangles=size(FastHGP.J_fz,1);

    % ***** preprocess *****
    J=[FastHGP.J_fz;FastHGP.J_fbz];
    JWithOnes=[J;ones(1,numCones)];% add one constraint to J matrix to fix translation : the avarage of all the cones iz zero

    %Use weighted least squares instead of least squares to account for the areas of triangles  
    W_small=[FastHGP.Area;FastHGP.Area]';
    W=[W_small, mean(FastHGP.Area)];
    MInv=inv((JWithOnes'.*W)*JWithOnes); %inverse of weighted J'*J
    weightedPinvJ=(MInv*JWithOnes').*W; %weighted least squares
    
    % ***** statring point *****
    start_fz=1./FastHGP.frames;
    start_fbz=zeros(size(FastHGP.frames));
    a_start=[start_fz;start_fbz;0];
    frames=FastHGP.frames;

    % ***** first do one MAP global step *****
    [ c , a_l ] = globalStep_MAP( a_start , weightedPinvJ , JWithOnes ); 
    fz_global=a_l(1:numTriangles);
    fbz_global=a_l(numTriangles+1:end);

    %variables are: 
    %c(1..numCones) complex representation of cones and meta vetices
    %J(1..2*|F|) f_z and f_zbar in tringles near cones and near meta

    % ***** parameters *****
    maxIter=500;
    sigma2Eps=0.01;
    k=0.9;
    tol=1e-4;

    % ***** ATP iterations *****
    for iter=1:maxIter

        %%%local step%%%
        [ fz_local , fbz_local ] = localStep( frames , sigma2Eps , k , fz_global , fbz_global );
        b_l=[fz_local ; fbz_local ];

        %%%check convegrence%%%
        ni=a_l-b_l;
        norm_ni_sqr=(ni'.*W_small)*ni;
        if(sqrt(norm_ni_sqr)<tol)%stop condition 
            break;
        end

        %%%global step%%%
        [ c , a_l ] = globalStep_ATP( weightedPinvJ , JWithOnes , c , [ni;0] , norm_ni_sqr , W );
        fz_global=a_l(1:numTriangles);
        fbz_global=a_l(numTriangles+1:end);

    end

    FastHGP.ATPtime=toc;
    FastHGP.ATPiters=iter;
    FastHGP.foldsInATP=find(abs(fz_global)<abs(fbz_global));

    numFolds=length(FastHGP.foldsInATP);
    FastHGP.ATPSucsess=(numFolds==0) && all(isnan(c)==0);
    FastHGP.ATPSucsess=double(FastHGP.ATPSucsess);
    
    FastHGP.UVonCones=[real(c);imag(c)];

end