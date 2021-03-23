function H = calcSolutionAffinGround_Optimal(A,P,Q)
  [C,b] = buildCoeffMtxAffinGround_Optimal(A,P,Q); 
  x = Solver_1Const(C,b);
  x = -x;
  %C*x-b
  if size(x) == 0
      H = calcSolutionAffinGround(A,P,Q);
  else
      H = [x(1) x(3) -x(2); 
           0 1 0;
           x(2) x(4) x(1)];
  end
end