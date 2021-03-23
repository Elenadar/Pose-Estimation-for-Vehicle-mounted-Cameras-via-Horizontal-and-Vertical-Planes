function b = calcSolutionLinear3pt(P,Q)
    A = buildCoeffMtxLinear3pt(P,Q);
    [U, S, V] = svd(A);
    sizeOfV = size(V);
    b = V(:,sizeOfV(1));
end