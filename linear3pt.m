function [R,t] = linear3pt(P,Q)
    b = calcSolutionLinear3pt(P,Q);

    phi_Linear3pt = atan2(b(2),b(1));
    theta_Linear3p = atan2(b(4),b(3))+phi_Linear3pt;

    R = [[cos(theta_Linear3p) 0 sin(theta_Linear3p)]; [0 1 0]; [-sin(theta_Linear3p) 0 cos(theta_Linear3p)]];

    t1 = [sin(phi_Linear3pt); 0 ; cos(phi_Linear3pt)];
    t2 = -t1;
    t = [t1, t2];
end