function H = calcSolutionAffinFrontWall_Optimal(A,P,Q)
  [C,b] = buildCoeffMtxAffinFrontWall_Optimal(A,P,Q);
  x = Solver_1Const(C,b);
  %C*x-b
  if size(x) == 0
      H = calcSolutionAffinWall(A,P,Q);
  elseif max(C(:)) > 100000
      H = calcSolutionAffinWall(A,P,Q);
  else
      H = [-x(1) 0 x(3)+x(2); 
           0 1 0;
           -x(2) 0 x(4)-x(1)];
  end
end