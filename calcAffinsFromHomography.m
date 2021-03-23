function A = calcAffinsFromHomography(H,P,Q)
  s=H(3,:)*P';
  
  a11=(H(1,1)-H(3,1).*Q(:,1))./s(:);
  a12=(H(1,2)-H(3,2).*Q(:,1))./s(:);
  a21=(H(2,1)-H(3,1).*Q(:,2))./s(:);
  a22=(H(2,2)-H(3,2).*Q(:,2))./s(:);
  A=[a11,a12,a21,a22];
end
