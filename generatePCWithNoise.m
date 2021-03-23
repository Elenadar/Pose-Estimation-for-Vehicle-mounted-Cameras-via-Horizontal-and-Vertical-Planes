function [pt2D_1,pt2D_2,pt2D_2_noise] = generatePCWithNoise(plane,P1,P2,C,noiseLength)
  sizeOfPlane = size(plane); 
  NUM = sizeOfPlane(1);
  
  pt2D_1 = [];
  pt2D_2 = [];
  for i = 1:NUM
    pc1 = P1*[transpose(plane(i,:));1];
    pc2 = P2*[transpose(plane(i,:));1];
    pc1 = pc1/pc1(3);
    pc2 = pc2/pc2(3);
    pt2D_1 = [pt2D_1;transpose(pc1)];
    pt2D_2 = [pt2D_2;transpose(pc2)];
  end 
    
  noiseAngles = 2*pi*rand(NUM,1);
  noiseVector = zeros(NUM,3);

  noiseVector(:,1) = cos(noiseAngles(:));
  noiseVector(:,2) = sin(noiseAngles(:));
   
  pt2D_2_noise = pt2D_2 + noiseLength * noiseVector;
end