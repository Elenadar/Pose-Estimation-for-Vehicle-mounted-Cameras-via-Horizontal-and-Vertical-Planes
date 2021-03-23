%NoiseScale: %
%noiseAngle: radian
function A_noisy=AddNoiseToAffine(A,noiseScale,noiseAngle)

    [angle1,angle2,s1,s2]=DecomposeAffine(A);
    
    angle1=angle1+noiseAngle*randn;
    angle2=angle2+noiseAngle*randn;
    
    s1=s1*(1.0+randn*noiseScale/100.0);
    s2=s2*(1.0+randn*noiseScale/100.0);
    
    U=[cos(angle1),-sin(angle1);sin(angle1),cos(angle1)];
    V=[cos(angle2),-sin(angle2);sin(angle2),cos(angle2)];
    S=[s1,0;0,s2];

    
    A_noisy=U*S*V';

end
