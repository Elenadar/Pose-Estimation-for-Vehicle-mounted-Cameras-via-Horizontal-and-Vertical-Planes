function C = buildCoeffMtxAffinWall(A,P,Q)
  C = [];
  sizeOfP = size(P);
  for k = 1:sizeOfP(1)
    C = [ C;
        [1,      0,       0,       -Q(k,1)-A(k,1)*P(k,1),   -A(k,1)];
        [0,      0,       0,       -A(k,2)*P(k,1),          -A(k,2)];
        [0,      0,       0,       -Q(k,2)-A(k,3)*P(k,1),   -A(k,3)];
        [0,      0,       1,       -A(k,4)*P(k,1),          -A(k,4)];
        [P(k,1), 1,       0,       -P(k,1)*Q(k,1),          -Q(k,1)];
        [0,      0,       P(k,2),  -P(k,1)*Q(k,2),          -Q(k,2)];
  ];
  end
end