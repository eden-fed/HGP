function [x, relative_error] = general_solve(M, rhs)

    x = M \ rhs; %solve the linear system
    
    %check the accuracy of the result
    residual_matrix = M*x-rhs;
	err = norm(residual_matrix, 'fro');
    avg =  mean(mean(abs(x)));
    relative_error = err/max(avg, eps);
end