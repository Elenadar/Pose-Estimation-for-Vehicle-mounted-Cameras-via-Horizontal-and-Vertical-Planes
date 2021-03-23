function [gamma,alpha] = calcParamsFromIntersections(x1,y1,a,b,R1,R2_T)
  eta = atan2(y1,x1);
  phi = atan2(a*y1,b*x1);
  R2_v1 = [cos(phi);sin(phi)];
  R1_T_v2 = [cos(eta);sin(eta)];
  v1 = R2_T*R2_v1;
  v2 = R1*R1_T_v2;
  gamma = atan2(v1(2),v1(1));
  alpha = atan2(v2(1),v2(2))+gamma;
end