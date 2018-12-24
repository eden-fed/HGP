function [ c , a_l ] = globalStep_MAP( a_l , MInvJtrans , J )

    c=MInvJtrans*a_l;
    a_l=J*c;
    a_l=a_l(1:end-1);

end

