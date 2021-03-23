function angleError = calcAngleErrorVectors(gt,est)
  gt = gt/norm(gt);
  est = est/norm(est);
  angleError = abs(180/pi*acos(clamp(dot(gt,est),0,1)));
end