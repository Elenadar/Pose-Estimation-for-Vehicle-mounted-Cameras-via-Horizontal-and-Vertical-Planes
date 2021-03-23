function H_2affin = calcSolution2Affin(A,P,Q)
  C = buildCoeffMtx2Affin(A,P,Q);
  [U,S,V] = svd(C);
  sizeOfV = size(V);
  h = V(:,sizeOfV(1));
  H_2affin = [[h(1),h(2),h(3)];[h(4),h(5),h(6)];[h(7),h(8),h(9)]];
  H_2affin = H_2affin/H_2affin(3,3); 
end