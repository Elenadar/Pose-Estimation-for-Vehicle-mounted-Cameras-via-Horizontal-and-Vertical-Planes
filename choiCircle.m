function [R,t] = choiCircle(P,Q)
    [A,B] = buildCoeffMtxChoiCircle(P,Q);

    C = pinv(B)*A;
    [U, S, V] = svd(transpose(C)*C);

    s1 = S(1,1);
    s2 = S(2,2);

    if s1 < 1
        y1 = [1;0];
        y2 = [1;0];
        y3 = [-1;0];
        y4 = [-1;0];
    elseif s2 > 1
         y1 = [0;1];
        y2 = [0;1];
        y3 = [0;-1];
        y4 = [0;-1];
    else
        y1 = [sqrt((1-s2)/(s1-s2));sqrt((s1-1)/(s1-s2))];
        y2 = [-sqrt((1-s2)/(s1-s2));sqrt((s1-1)/(s1-s2))];
        y3 = [sqrt((1-s2)/(s1-s2));-sqrt((s1-1)/(s1-s2))];
        y4 = [-sqrt((1-s2)/(s1-s2));-sqrt((s1-1)/(s1-s2))];
    end

    a1 = U*y1;
    a2 = U*y2;
    a3 = U*y3;
    a4 = U*y4;

    b1 = C*a1;
    b2 = C*a2;
    b3 = C*a3;
    b4 = C*a4;
    
    phi_choi = [atan2(a1(2),a1(1)), atan2(a2(2),a2(1)),atan2(a3(2),a3(1)),atan2(a4(2),a4(1))];
    theta_choi = [atan2(b1(2),b1(1))+phi_choi(1), atan2(b2(2),b2(1))+phi_choi(2),atan2(b3(2),b3(1))+phi_choi(3),atan2(b4(2),b4(1))+phi_choi(4)];
    
    R1_choi = [[cos(theta_choi(1)) 0 sin(theta_choi(1))]; [0 1 0]; [-sin(theta_choi(1)) 0 cos(theta_choi(1))]]';
    R2_choi = [[cos(theta_choi(2)) 0 sin(theta_choi(2))]; [0 1 0]; [-sin(theta_choi(2)) 0 cos(theta_choi(2))]]';
    R3_choi = [[cos(theta_choi(3)) 0 sin(theta_choi(3))]; [0 1 0]; [-sin(theta_choi(3)) 0 cos(theta_choi(3))]]';
    R4_choi = [[cos(theta_choi(4)) 0 sin(theta_choi(4))]; [0 1 0]; [-sin(theta_choi(4)) 0 cos(theta_choi(4))]]';

    
    t1_choi = [a1(2); 0 ; a1(1)];
    t2_choi = [a2(2); 0 ; a2(1)];
    t3_choi = [a3(2); 0 ; a3(1)];
    t4_choi = [a4(2); 0 ; a4(1)];

    R = [R1_choi,R2_choi,R3_choi,R4_choi];
    t = [t1_choi,t2_choi,t3_choi,t4_choi];
end