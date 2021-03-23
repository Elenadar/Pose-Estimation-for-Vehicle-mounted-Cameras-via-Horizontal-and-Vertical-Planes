function H_affinGround = calcSolutionAffinGround(A,P,Q)
  
  C = buildCoeffMtxAffinGround(A,P,Q);
  [U,S,V] = svd(C);
  sizeOfV = size(V);
  h = V(:,sizeOfV(1));
  H_affinGround = [[h(1),h(2),-h(3)];
       [0,h(5),0];
       [h(3),h(4),h(1)]];
  H_affinGround = H_affinGround/H_affinGround(2,2);
end