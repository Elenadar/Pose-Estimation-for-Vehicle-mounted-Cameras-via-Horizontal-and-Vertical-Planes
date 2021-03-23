function [A,B] = buildCoeffMatricesChoiLine(P,Q)
  sizeOfP = size(P);
  n = sizeOfP(1);
  A = zeros(n,2);
  B = zeros(n,2);
  for i=1:n
      [A(i,:), B(i,:)] = buildChoiRow(P(i,:),Q(i,:));
  end
 end