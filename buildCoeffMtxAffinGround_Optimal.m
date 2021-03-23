function [C,b] = buildCoeffMtxAffinGround_Optimal(A,P,Q)
  sizeOfP = size(P);
  n = sizeOfP(1)*6;
  C = zeros(n,4);
  b = zeros(n,1);
  
  C_rapid = buildCoeffMtxAffinGround(A,P,Q);
  
  b = C_rapid(:,5);
  C(:,1) = C_rapid(:,1);
  C(:,3) = C_rapid(:,2);
  C(:,2) = C_rapid(:,3);
  C(:,4) = C_rapid(:,4);
end