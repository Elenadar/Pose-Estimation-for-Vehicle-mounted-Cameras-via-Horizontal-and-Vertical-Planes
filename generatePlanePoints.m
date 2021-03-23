function points = generatePlanePoints(NUM,r,t,n)
  n_T = transpose(n);
  v = null(n_T);
  v1 = v(:,1);
  v2 = v(:,2);
  points = [];
  for idx = 1:NUM
    alpha = rand();
    point = rand()*v1+rand()*v2-(v1+v2)*0.5;
    points = [points, point];
  end
  
  points = r*points;
  points = [points',ones(NUM,1)];
  points = ([eye(3),t]*points')';
end
