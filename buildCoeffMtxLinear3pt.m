function C = buildCoeffMtxLinear3pt(P,Q)
  sizeOfP = size(P);
  numOfPoints = sizeOfP(1);
  for i = 1:numOfPoints
    C(i,1) = P(i,1)*Q(i,2);
    C(i,2) = -P(i,3)*Q(i,2);
    C(i,3) = -Q(i,1)*P(i,2);
    C(i,4) = -Q(i,3)*P(i,2);
    
    %C(i,1) = P(i,1)*Q(i,2);
    %C(i,2) = -P(i,3)*Q(i,2);
    %C(i,3) = -Q(i,1)*P(i,2);
    %C(i,4) = -Q(i,3)*P(i,2);
  end
end