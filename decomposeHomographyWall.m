function [R_decomp,t_decomp,n_decomp] = decomposeHomographyWall(H)
  %normalizing the H(2,2) to be 1.
  scale = H(2,2);
  H = H/scale;

  h1 = H(1,1);
  h2 = H(1,3);
  h3 = H(3,1);
  h4 = H(3,3);
  
  B = [h2 -h1; h4 -h3];

  [R1,S,R2] = svd(B);

  %calculating the intersection points of the ellipse and the circle
  a = S(1,1);
  b = S(2,2);
  if a*a < 1
    i_x1 = 1;
    i_x2 = -1;
    i_y1 = 0;
    i_y2 = 0;
  elseif a*a >=1 && b*b <= 1
    i_x1 = a*sqrt((1-b*b)/(a*a-b*b));
    i_x2 = -a*sqrt((1-b*b)/(a*a-b*b));
    i_y1 = b*sqrt((a*a-1)/(a*a-b*b));
    i_y2 = -b*sqrt((a*a-1)/(a*a-b*b));
  else
    i_x1 = 0;
    i_x2 = 0;
    i_y1 = 1;
    i_y2 = -1;
  end

  %calculating the angle of rotation and the angle of the normal vector
  [gamma1,alpha1] = calcParamsFromIntersections(i_x1,i_y1,a,b,R1,transpose(R2));
  [gamma2,alpha2] = calcParamsFromIntersections(i_x2,i_y1,a,b,R1,transpose(R2));
  [gamma3,alpha3] = calcParamsFromIntersections(i_x1,i_y2,a,b,R1,transpose(R2));
  [gamma4,alpha4] = calcParamsFromIntersections(i_x2,i_y2,a,b,R1,transpose(R2));


  R1 = [[cos(alpha1) 0 sin(alpha1)];[0 1 0];[-sin(alpha1) 0 cos(alpha1)]];
  R2 = [[cos(alpha2) 0 sin(alpha2)];[0 1 0];[-sin(alpha2) 0 cos(alpha2)]];
  R3 = [[cos(alpha3) 0 sin(alpha3)];[0 1 0];[-sin(alpha3) 0 cos(alpha3)]];
  R4 = [[cos(alpha4) 0 sin(alpha4)];[0 1 0];[-sin(alpha4) 0 cos(alpha4)]];

  n1 = [cos(gamma1); 0 ; sin(gamma1)];
  n2 = [cos(gamma2); 0 ; sin(gamma2)];
  n3 = [cos(gamma3); 0 ; sin(gamma3)];
  n4 = [cos(gamma4); 0 ; sin(gamma4)];

  %Using the alternative solution if needed.
  if abs(cos(gamma1)) > abs(sin(gamma1))
    x1 = (-h1+cos(alpha1))/cos(gamma1);
    y1 = (-h3-sin(alpha1))/cos(gamma1);
  else
    x1 = (-h2+sin(alpha1))/sin(gamma1);
    y1 = (-h4+cos(alpha1))/sin(gamma1);
  end

  if abs(cos(gamma2)) > abs(sin(gamma2))
    x2 = (-h1+cos(alpha2))/cos(gamma2);
    y2 = (-h3-sin(alpha2))/cos(gamma2);
  else
    x2 = (-h2+sin(alpha2))/sin(gamma2);
    y2 = (-h4+cos(alpha2))/sin(gamma2);
  end

  if abs(cos(gamma3)) > abs(sin(gamma3))
    x3 = (-h1+cos(alpha3))/cos(gamma3);
    y3 = (-h3-sin(alpha3))/cos(gamma3);
  else
    x3 = (-h2+sin(alpha3))/sin(gamma3);
    y3 = (-h4+cos(alpha3))/sin(gamma3);
  end

  if abs(cos(gamma4)) > abs(sin(gamma4))
    x4 = (-h1+cos(alpha4))/cos(gamma4);
    y4 = (-h3-sin(alpha4))/cos(gamma4);
  else
    x4 = (-h2+sin(alpha4))/sin(gamma4);
    y4 = (-h4+cos(alpha4))/sin(gamma4);
  end


  t1 = [x1; 0 ; y1];
  t2 = [x2; 0 ; y2];
  t3 = [x3; 0 ; y3];
  t4 = [x4; 0 ; y4];
  
  
  R_decomp = [R1,R2,R3,R4];
  t_decomp = [t1,t2,t3,t4];
  n_decomp = [n1,n2,n3,n4];
end