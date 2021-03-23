function [angle1,angle2,s1,s2]=DecomposeAffine(A)
    [U S V]=svd(A);
    angle1=atan2(U(2,1),U(1,1));
    angle2=atan2(V(2,1),V(1,1));
    s1=S(1,1);
    s2=S(2,2);
    
    
    if (det(U)<0.0)
        s2=-s2;
    end
    
    if (det(V)<0.0)
        s2=-s2;
    end
        
end
