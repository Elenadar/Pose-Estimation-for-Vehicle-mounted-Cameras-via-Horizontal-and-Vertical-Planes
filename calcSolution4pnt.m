function H_4pt = calcSolution4pnt(P,Q)
  A = buildCoeffMtx4pnt(P,Q);
  [U,S,V] = svd(A);
  sizeOfV = size(V);
  h = V(:,sizeOfV(1));
  H_4pt = [[h(1),h(2),h(3)];[h(4),h(5),h(6)];[h(7),h(8),h(9)]];
  H_4pt = H_4pt/H_4pt(3,3); 
end