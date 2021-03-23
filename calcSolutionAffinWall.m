function H_affin = calcSolutionAffinFrontWall_Optimal(A,P,Q)
  C = buildCoeffMtxAffinWall(A,P,Q);
  [U,S,V] = svd(C);
  sizeOfV = size(V); 
  h = V(:,sizeOfV(1));
  H_affin = [[h(1),0,h(2)];
      [0,h(3),0];
      [h(4),0,h(5)]];
end