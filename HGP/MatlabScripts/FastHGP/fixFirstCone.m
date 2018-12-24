function [ x, FixedIndices, FixedValues ] = fixFirstCone( initialValue )
%fix the first cone - translation degree of freedom
    n=length(initialValue)/2; 
    FixedIndices = [1;1+n]; 
    FixedValues = [initialValue(1);initialValue(1+n)]; 
    x=initialValue(setdiff(1:end,FixedIndices));
end

