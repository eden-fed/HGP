function [ c_i_p_1 , a_l ] = globalStep_ATP( MInvJtrans , J , c_i , ni, norm_ni_sqr , W )

    M_J_ni = MInvJtrans*ni;
    
    scalar = norm_ni_sqr/(ni'*W*J*M_J_ni);
    
    c_i_p_1 = c_i - scalar*M_J_ni;
    
    a_l=J*c_i_p_1;
    a_l=a_l(1:end-1);

end


