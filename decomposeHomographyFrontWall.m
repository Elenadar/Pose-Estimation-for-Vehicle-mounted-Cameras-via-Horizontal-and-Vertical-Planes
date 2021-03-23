function [R,t] = decomposeHomographyFrontWall(H)
  %normalizing the H(2,2) to be 1.
  scale = H(2,2);
  H = H/scale;

  h1 = H(1,1);
  h2 = H(1,3);
  h3 = H(3,1);
  h4 = H(3,3);
  
  alpha = atan2(h3,h1);

  R = [[cos(alpha) 0 -sin(alpha)];[0 1 0];[sin(alpha) 0 cos(alpha)]];


  p = -h2 - sin(alpha);
  q = -h4 + cos(alpha);

  
  beta = atan2(q,p);
  t = [cos(beta) ; 0 ; sin(beta)];
end