function [R,t,n] = decomposeHomographyGround(H)
  %normalizing the H(2,2) to be 1.
  scale = H(2,2);
  H = H/scale;

  h1 = H(1,1);
  h2 = H(1,2);
  h3 = -H(1,3);
  h4 = H(3,2);
  h5 = H(2,2);
  
  alpha = atan2(h3,h1);

  %R = [[H(1,1) 0 H(1,3)];[0 1 0];[H(3,1) 0 H(3,3)]];
  R = [[cos(alpha) 0 -sin(alpha)];[0 1 0];[sin(alpha) 0 cos(alpha)]];

  n = [0; 1 ; 0];

  p = h2;
  q = h4;

  beta = atan2(q,p);

  %t1 = -[H(1,2); 0 ; H(3,2)]
  %t = t1/norm(t1);
  t = -[cos(beta) ; 0 ; sin(beta)];

end