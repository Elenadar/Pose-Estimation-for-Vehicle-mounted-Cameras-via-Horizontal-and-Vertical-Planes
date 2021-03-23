function pts=IntersectionTwoQudratics(first, second)

%l1*%alpha^2+l2*%alpha*%*%beta+l3*beta^2+l4*alpha+l5*beta+l6


k1=first(1);
k2=first(2);
k3=first(3);
k4=first(4);
k5=first(5);
k6=first(6);

l1=second(1);
l2=second(2);
l3=second(3);
l4=second(4);
l5=second(5);
l6=second(6);



n1=4*k3*k3*l1-2*k3*l2*k2+2*l3*k2*k2-4*l3*k3*k1;



n2=-2*k3*l5*k2-2*k3*l2*k5+4*l3*k2*k5-4*l3*k3*k4+4*k3*k3*l4;

n3=-2*k3*l5*k5+2*l3*k5*k5-4*l3*k3*k6+4*k3*k3*l6;



%Left side:

left4=n1*n1;
left3=2*n1*n2;
left2=n2*n2+2*n1*n3;
left1=2*n2*n3;
left0=n3*n3;
left=[left4,left3,left2,left1,left0];


%Right side 
poly1=[k2*k2-4*k3*k1,2*k2*k5-4*k3*k4,k5*k5-4*k3*k6];

poly2=[(-l3*k2+k3*l2)*(-l3*k2+k3*l2),2*(-l3*k2+k3*l2)*(-l3*k5+k3*l5),(-l3*k5+k3*l5)*(-l3*k5+k3*l5)];

right=4*conv(poly1,poly2);

poly=left-right;

rss=roots(poly);


cnt=1;
pts=[];

for idx=1:length(rss)
    if (isreal(rss(idx)))
        a=rss(idx);
        b1=(-1*(k2*a+k5)+sqrt((k2*a+k5)*(k2*a+k5)-4*k3*(k1*a*a+k4*a+k6)))/(2*k3);
        b2=(-1*(k2*a+k5)-sqrt((k2*a+k5)*(k2*a+k5)-4*k3*(k1*a*a+k4*a+k6)))/(2*k3);
        
        err1=k1*a*a+k2*a*b1+k3*b1*b1+k4*a+k5*b1+k6;
        err2=l1*a*a+l2*a*b1+l3*b1*b1+l4*a+l5*b1+l6;
        
        if ((abs(err1)<1e-7)&&abs(err2)<1e-7)
            pts (cnt,1:2)=[a,b1];
            cnt=cnt+1;
        end
        
                
        err1=k1*a*a+k2*a*b2+k3*b2*b2+k4*a+k5*b2+k6;
        err2=l1*a*a+l2*a*b2+l3*b2*b2+l4*a+l5*b2+l6;

        if ((abs(err1)<1e-7)&&abs(err2)<1e-7)
            pts(cnt,1:2)=[a,b2];
            cnt=cnt+1;
        end
        
    end
end






end
