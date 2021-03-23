function [best_R_Error,best_t_Error] = calcAlphaBetaErrors(R_est,t_est,R_gt,t_gt)
  best_R_Error = 1000;
  best_t_Error = 1000;
  
  sizeofEstimationsT = size(t_est);
  numOfT = sizeofEstimationsT(2);
  

  
  sizeofEstimationsR = size(R_est);
  if sizeofEstimationsR(1) > sizeofEstimationsR(2)
      R_est = R_est';
  end
  sizeofEstimationsR = size(R_est);
  numOfR = sizeofEstimationsR(2);
  t_best = [];
  R_best = [];
  for i = 1:numOfT
    t = t_est(:,i);
    if numOfR == 3
        R = R_est;
    else
        R = R_est(1:3,(i-1)*3+1:i*3);
    end
    
    t_Error = calcAngleErrorVectors(t_gt,t);

    if t_Error < best_t_Error
      t_best = t;
      best_t_Error = t_Error;
    end
    if trace(transpose(R)*R_gt) < -1
      R = -R;
    end
    %clampet használjak
    
    R_Error = 180/pi*acos(clamp((trace(transpose(R_gt)*R)-1)/2,0,1));
    
    if R_Error < best_R_Error
      if isreal(R_Error)
          R_best = R;
          best_R_Error = R_Error;
      end
    end
  end

end