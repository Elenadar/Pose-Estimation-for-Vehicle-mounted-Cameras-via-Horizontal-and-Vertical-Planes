function d = distanceFromPlane(planeNormal,planePoint,point)
  D = dot(planeNormal,planePoint);
  d = abs(planeNormal(1)*point(1)+planeNormal(2)*point(2)+planeNormal(3)*point(3)+D)/norm(planeNormal);
end
