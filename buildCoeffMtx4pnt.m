function C = buildCoeffMtx4pnt(P,Q)

  sizeOfP = size(P); 
  n = sizeOfP(1)*2;
  C = zeros(n,9);
  C(1:2:n,1) = P(:,1);
  C(1:2:n,2) = P(:,2);
  C(1:2:n,3) = 1;
  %C(1:2:n,4) = 0;
  %C(1:2:n,5) = 0;
  %C(1:2:n,6) = 0;
  C(1:2:n,7) = -P(:,1).*Q(:,1);
  C(1:2:n,8) = -P(:,2).*Q(:,1);
  C(1:2:n,9) = -Q(:,1);
  
  %C(2:2:n,1) = 0;
  %C(2:2:n,2) = 0;
  %C(2:2:n,3) = 0;
  C(2:2:n,4) = P(:,1);
  C(2:2:n,5) = P(:,2);
  C(2:2:n,6) = 1;
  C(2:2:n,7) = -P(:,1).*Q(:,2);
  C(2:2:n,8) = -P(:,2).*Q(:,2);
  C(2:2:n,9) = -Q(:,2);
end