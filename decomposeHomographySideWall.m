function [R_decomp,t_decomp,n_decomp] = decomposeHomographySideWall(H)
  %normalizing the H(2,2) to be 1.
  scale = H(2,2);
  H = H/scale;

  h1 = H(1,1);
  h2 = H(1,3);
  h3 = H(3,1);
  h4 = H(3,3);
  
  alpha = atan2(-h2,h4);

  R = [[cos(alpha) 0 sin(alpha)];[0 1 0];[-sin(alpha) 0 cos(alpha)]];

  n = [1; 0 ; 0];

  p = -h1 + cos(alpha);
  q = -h3 + sin(alpha);


  t1 = [-p; 0 ; -q];
  t = t1/norm(t1);
  
  
  R_decomp = [R,R,R,R];
  t_decomp = [t,t,t,t];
  n_decomp = [n,n,n,n];
end