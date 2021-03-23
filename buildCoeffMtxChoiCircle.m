function [A,B] = buildCoeffMtxChoiCircle(P,Q)
    sizeOfP = size(P);
    numOfPoints = sizeOfP(1);
    A = zeros(numOfPoints,2);
    B = zeros(numOfPoints,2);
    for i = 1: numOfPoints
        A(i,1) = P(i,1)*Q(i,2);
        A(i,2) = -P(i,3)*Q(i,2);
        B(i,1) = Q(i,1)*P(i,2);
        B(i,2) = Q(i,3)*P(i,2);
    end
end