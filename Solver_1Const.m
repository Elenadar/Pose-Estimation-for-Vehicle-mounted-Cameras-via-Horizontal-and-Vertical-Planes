%Solves the problem argmin(Ax-b) w.r.t. x, subject to x(1)^2+x(2)^2=1
function vec=Solver_1Const(A,b)
    sizeOfA = size(A);
    N = sizeOfA(2);
    M = sizeOfA(1);
    
    specMtx=zeros(N,2);
    specMtx(1,1)=1.0;
    specMtx(2,2)=1.0;
    mtx=[A'*A,specMtx,-A'*b];
    
    nullVcs=null(mtx);
    n1=nullVcs(:,1);
    n2=nullVcs(:,2);
    n3=nullVcs(:,3);
    
    


    sc_alpha=-n1(N+3)/n3(N+3);
    sc_beta=-n2(N+3)/n3(N+3);
    sc_null=1/n3(N+3);
    
    
    x1poly=[n1(1)+n3(1)*sc_alpha,n2(1)+n3(1)*sc_beta,n3(1)*sc_null];
    x2poly=[n1(2)+n3(2)*sc_alpha,n2(2)+n3(2)*sc_beta,n3(2)*sc_null];
    xNp1poly=[n1(N+1)+n3(N+1)*sc_alpha,n2(N+1)+n3(N+1)*sc_beta,n3(N+1)*sc_null];
    xNp2poly=[n1(N+2)+n3(N+2)*sc_alpha,n2(N+2)+n3(N+2)*sc_beta,n3(N+2)*sc_null];
    
    
    
    A1left=xNp1poly(1)*x2poly(1);
    B1left=xNp1poly(1)*x2poly(2)+xNp1poly(2)*x2poly(1);
    C1left=xNp1poly(2)*x2poly(2);
    D1left=xNp1poly(1)*x2poly(3)+xNp1poly(3)*x2poly(1);
    E1left=xNp1poly(2)*x2poly(3)+xNp1poly(3)*x2poly(2);
    F1left=xNp1poly(3)*x2poly(3);

    A1right=xNp2poly(1)*x1poly(1);
    B1right=xNp2poly(1)*x1poly(2)+xNp2poly(2)*x1poly(1);
    C1right=xNp2poly(2)*x1poly(2);
    D1right=xNp2poly(1)*x1poly(3)+xNp2poly(3)*x1poly(1);
    E1right=xNp2poly(2)*x1poly(3)+xNp2poly(3)*x1poly(2);
    F1right=xNp2poly(3)*x1poly(3);
    
    
    A1=A1left-A1right;
    B1=B1left-B1right;
    C1=C1left-C1right;
    D1=D1left-D1right;
    E1=E1left-E1right;
    F1=F1left-F1right;

    

    A2first=x1poly(1)*x1poly(1);
    B2first=2*x1poly(1)*x1poly(2);
    C2first=x1poly(2)*x1poly(2);
    D2first=2*x1poly(1)*x1poly(3);
    E2first=2*x1poly(2)*x1poly(3);
    F2first=x1poly(3)*x1poly(3);

    A2second=x2poly(1)*x2poly(1);
    B2second=2*x2poly(1)*x2poly(2);
    C2second=x2poly(2)*x2poly(2);
    D2second=2*x2poly(1)*x2poly(3);
    E2second=2*x2poly(2)*x2poly(3);
    F2second=x2poly(3)*x2poly(3);
    
    A2=A2first+A2second;
    B2=B2first+B2second;
    C2=C2first+C2second;
    D2=D2first+D2second;
    E2=E2first+E2second;
    F2=F2first+F2second-1.0;
    
    
    first=[A1,B1,C1,D1,E1,F1];
    second=[A2,B2,C2,D2,E2,F2];
    
    params=IntersectionTwoQudratics(first, second);
    
    
    errMin=inf;
    vec=[];
    sizeOfParams = size(params);
    for idx=1:sizeOfParams(1)
        alpha=params(idx,1);
        beta=params(idx,2);
        gamma=1/n3(N+3)-n1(N+3)*alpha/n3(N+3)-n2(N+3)*beta/n3(N+3);
        
        cand=alpha*n1+beta*n2+gamma*n3;
        cand=cand(1:N);
        err=norm(A*cand-b);

        if (err<errMin)
            vec=cand;
            errMin=err;
        end   
    end
end
