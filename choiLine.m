function [R,t,deg] = choiLine(x1,x2)
    [A_r,B_r] = buildCoeffMatricesChoiLine(x1,x2);

    C = pinv(B_r)*A_r;
    
    CTC = transpose(C)*C;
    alpha = CTC(1,1);
    beta = CTC(1,2);
    gamma = CTC(2,2);


    alpha_ = alpha-gamma;
    gamma_ = alpha+gamma-2;
    beta_ = 2*beta;
    delta = sqrt(alpha_*alpha_ + beta_*beta_ - gamma_*gamma_);

    cos2phi_1 = 1/(alpha_*alpha_+beta_*beta_)*(-alpha_*gamma_+beta_*delta);
    cos2phi_2 = 1/(alpha_*alpha_+beta_*beta_)*(-alpha_*gamma_-beta_*delta);


    sin2phi_1 = 1/(alpha_*alpha_+beta_*beta_)*(-beta_*gamma_-alpha_*delta);
    sin2phi_2 = 1/(alpha_*alpha_+beta_*beta_)*(-beta_*gamma_+alpha_*delta);


    a1 = 1/sqrt(2)*[sqrt(1+cos2phi_1);sign(sin2phi_1)*sqrt(1-cos2phi_1)];
    a2 = -1/sqrt(2)*[sqrt(1+cos2phi_1);sign(sin2phi_1)*sqrt(1-cos2phi_1)];
    a3 = 1/sqrt(2)*[sqrt(1+cos2phi_2);sign(sin2phi_2)*sqrt(1-cos2phi_2)];
    a4 = -1/sqrt(2)*[sqrt(1+cos2phi_2);sign(sin2phi_2)*sqrt(1-cos2phi_2)];

    a = [a1 a2 a3 a4];
    deg = false;
    if isreal(a) == false
        deg = true;
        cos2phi_1 = sign(gamma_)/sqrt(alpha_*alpha_+beta_*beta_)*(-alpha_);

        sin2phi_1 = sign(gamma_)/sqrt(alpha_*alpha_+beta_*beta_)*(-beta_);


        a1 = 1/sqrt(2)*[sqrt(1+cos2phi_1);sign(sin2phi_1)*sqrt(1-cos2phi_1)];
        a2 = -1/sqrt(2)*[sqrt(1+cos2phi_1);sign(sin2phi_1)*sqrt(1-cos2phi_1)];
        a3 = a1;
        a4 = a2;

        a = [a1 a2 a3 a4];
    end
    
    b = C*a;

    b1 = b(:,1);
    b2 = b(:,2);
    b3 = b(:,3);
    b4 = b(:,4);

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