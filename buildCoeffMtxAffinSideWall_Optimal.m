function [C,b] = buildCoeffMtxAffinSideWall_Optimal(A,P,Q)
  C_rapid = buildCoeffMtxAffinWall(A,P,Q);
  C(:,1) = C_rapid(:,1) + C_rapid(:,5);
  C(:,2) = -C_rapid(:,2) + C_rapid(:,4);

  C(:,3) = -C_rapid(:,1);
  C(:,4) = -C_rapid(:,4);

  b = C_rapid(:,3);
end