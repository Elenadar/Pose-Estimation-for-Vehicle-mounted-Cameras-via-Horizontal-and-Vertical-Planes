%number of points used for estimating the homography matrix
numOfPoints = 4;

%testCase 0 : Ground plane
%testCase 1 : Wall in the fron
%testCase 2 : Wall on the side
%testCase 3 : Purely forward movement, random wall in the front
%testCase 4 : Purely forward movement, ground plane
%testCase 5 : Ground plane, gravity
%testCase 6 : Random wall in the front, gravity 

testCase = 0;


destFolder = 'output\\';
if ~exist(destFolder, 'dir')
       mkdir(destFolder)
end

%1: use only purely forward motion
%2: random planar motion
purelyForwardMovement = 0;

%number of iteration
iterNum = 1000;
numOfMethods = 7;

%Data structures to store errors
errors_stat = zeros(1,numOfMethods*4);
errors = zeros(iterNum,numOfMethods*2);

%Maximal errors
maxNoiseInPixels = 1;
maxNoiseScale = 1;
maxNoiseAngle = 1/180*pi;
maxGravityRotAngle = 5;

%Don't save results to file, only plot
plotOnly  = 0;

%Don't have to set these.
numOfSteps = 11;
noiseArray = zeros(numOfSteps,1);


%These are controlled by the testCase variable
proposedMethod = 0;
planeVariation = 0;
use3pt = 0;
useOptimal = 0;
gravity = 0;

switch testCase
    case 0 %ground
        proposedMethod = 0;
        planeVariation = 0;
        purelyForwardMovement = 0;
        use3pt = 1;
        useOptimal = 1;
        propColor = 'red';
        propMarker = 'o';
        propColorO = [0.58 0 0.83];
        propLineStyle = ':';
        propMarkerO = 'v';
        propMethodRapid = '1AC Ground';
        propMethodOptimal = '1AC Ground Optimal';
    case 1 %front wall
        proposedMethod = 1;
        planeVariation = 1;
        purelyForwardMovement = 0;
        use3pt = 0;
        useOptimal = 1;
        propColor = 'magenta';
        propColorO = [148,0,211]/255;
        propMarker = '^';
        propLineStyle = '-';
        propMarkerO = 'v';
        propMethodRapid = '1AC FV';
        propMethodOptimal = '1AC FV Optimal';
    case 2 %side wall
        proposedMethod = 2;
        planeVariation = 2;
        purelyForwardMovement = 0;
        use3pt = 0;
        useOptimal = 1;
        propColor = 'magenta';
        propColorO = [255,69,0]/255;
        propMarker = '^';
        propLineStyle = '-';
        propMarkerO = 'v';
        propMethodRapid = '1AC SV';
        propMethodOptimal = '1AC SV Optimal';
    case 3 %front wall random normal
        proposedMethod = 3;
        planeVariation = 3;
        purelyForwardMovement = 1;
        use3pt = 0;
        useOptimal = 0;
        propColor = 'magenta';
        propColorO = [0.58 0 0.83];
        propMarker = '^';
        propLineStyle = '-';
        propMarkerO = 'v';
        propMethodRapid = '1AC Vertical';
     case 4 %ground
        proposedMethod = 0;
        planeVariation = 0;
        purelyForwardMovement = 1;
        use3pt = 1;
        useOptimal = 1;
        propColor = 'red';
        propColorO = [0.58 0 0.83];
        propMarker = 'o';
        propLineStyle = ':';
        propMarkerO = 'v';
        propMethodRapid = '1AC Ground';
        propMethodOptimal = '1AC Ground Optimal';
    case 5
        proposedMethod = 0;
        planeVariation = 0;
        use3pt = 1;
        gravity = 1;
        useOptimal = 1;
        propColor = 'red';
        propColorO = [0.58 0 0.83];
        propMarker = 'o';
        propLineStyle = ':';
        propMarkerO = 'v';
        propMethodRapid = '1AC Ground';
        propMethodOptimal = '1AC Ground Optimal';
    case 6
        proposedMethod = 3;
        planeVariation = 3;
        use3pt = 0;
        gravity = 1;
        useOptimal = 0;
        propColor = 'magenta';
        propColorO = [0.58 0 0.83];
        propMarker = '^';
        propMarkerO = 'v';
        propLineStyle = '-';
        propMethodRapid = '1AC Vertical';
end


    for k = 1:numOfSteps

        (k-1)/(numOfSteps-1)

        if gravity == 0
            noiseInPixels = maxNoiseInPixels*(k-1)/(numOfSteps-1);
            noiseArray(k) = noiseInPixels;
            noiseScale = maxNoiseScale*(k-1)/(numOfSteps-1);
            noiseAngle = maxNoiseAngle*(k-1)/(numOfSteps-1);
        else
            gravityRotAngle = maxGravityRotAngle*(k-1)/(numOfSteps-1)
            noiseInPixels = maxNoiseInPixels;
            noiseArray(k) = gravityRotAngle;
            noiseScale = maxNoiseScale;
            noiseAngle = maxNoiseAngle;
        end
        
        if plotOnly == 0

        for j = 1:iterNum
            phi = 2*pi*rand();
            theta = 2*pi*rand();
            if purelyForwardMovement == 1 %purely forward movement
                phi = 0;
                theta = 0;
            end
            t1 = [0;0;0];
            R1 = eye(3);
            t2 = [sin(phi);0;cos(phi)];
            R2 = [cos(theta) 0 sin(theta); 0 1 0;-sin(theta) 0 cos(theta)];
            if gravity == 1
                R2_orig = R2;
                gravityRotAngleRad = randn*gravityRotAngle/180*pi;
                R_z = [[cos(gravityRotAngleRad) -sin(gravityRotAngleRad) 0];[sin(gravityRotAngleRad) cos(gravityRotAngleRad) 0];[0 0 1]];
                R_x = [[1 0 0];[0 cos(gravityRotAngleRad) -sin(gravityRotAngleRad)];[0 sin(gravityRotAngleRad) cos(gravityRotAngleRad)]];
                R2 = [[cos(theta) 0 sin(theta)];[0 1 0];[-sin(theta) 0 cos(theta)]];
                R2 = R_x*R_z*R2;
            else
                R2 = [[cos(theta) 0 sin(theta)];[0 1 0];[-sin(theta) 0 cos(theta)]];
            end

            d = 1;

            switch planeVariation
                case 0 %ground
                    X1 = [randn();-d;randn()-10;1];
                    X2 = [randn();-d;randn()-10;1];
                    X3 = [randn();-d;randn()-10;1];
                    X4 = [randn();-d;randn()-10;1];
                    n = [0;1;0];
                case 1 %front wall
                    d = 10;
                    X1 = [randn();randn();-d;1];
                    X2 = [randn();randn();-d;1];
                    X3 = [randn();randn();-d;1];
                    X4 = [randn();randn();-d;1];
                    n = [0;0;1];
                case 2 %side wall
                    d = 10;
                    X1 = [-d;randn();randn();1];
                    X2 = [-d;randn();randn();1];
                    X3 = [-d;randn();randn();1];
                    X4 = [-d;randn();randn();1];
                    n = [1;0;0];
                case 3 %front wall random normal
                    l = -10;
                    gamma = 80/180*pi*(rand-0.5)*2;
                    d = cos(gamma)*abs(l);

                    R_y = [cos(gamma) 0 -sin(gamma) 0;
                        0 1 0 0;
                        sin(gamma) 0 cos(gamma) 0;
                        0 0 0 1];
                    l = -10;
                    L = l*[0;0;1;0];

                    X1 = [randn();randn();0;1];
                    X2 = [randn();randn();0;1];
                    X3 = [randn();randn();0;1];
                    X4 = [randn();randn();0;1];
                    X1 = R_y*X1+L;
                    X2 = R_y*X2+L;
                    X3 = R_y*X3+L;
                    X4 = R_y*X4+L;

                    if gamma > 0
                        n = [-sin(gamma);0;cos(gamma)];
                    else
                        n = [-sin(gamma);0;cos(gamma)];
                    end

            end

            X = [X1,X2,X3,X4];
            %focal lenght in pixels
            f = 1000;
            %Principal point
            u0 = 960;%640;%960;
            v0 = 540;%360;%540;

            %Camera calibration matrix
            K = [[f 0 u0]; [0 f v0]; [0 0 1]];

            P1 = K*[eye(3) ,zeros(3,1)]*[[R1 , t1]; [zeros(1,3), 1]];
            P2 = K*[eye(3) ,zeros(3,1)]*[[R2 , t2]; [zeros(1,3), 1]];

            p1 = P1*X1;
            p2 = P1*X2;
            p3 = P1*X3;
            p4 = P1*X4;
            p1 = p1/p1(3);
            p2 = p2/p2(3);
            p3 = p3/p3(3);
            p4 = p4/p4(3);

            pt = [p1,p2,p3,p4]';

            q1 = P2*X1 + noiseInPixels*[randn();randn();0];
            q2 = P2*X2 + noiseInPixels*[randn();randn();0];
            q3 = P2*X3 + noiseInPixels*[randn();randn();0];
            q4 = P2*X4 + noiseInPixels*[randn();randn();0];
            q1 = q1/q1(3);
            q2 = q2/q2(3);
            q3 = q3/q3(3);
            q4 = q4/q4(3);

            qt = [q1,q2,q3,q4]';

            p1_IK = inv(K)*p1;
            p2_IK = inv(K)*p2;
            p3_IK = inv(K)*p3;
            p4_IK = inv(K)*p4;

            q1_IK = inv(K)*q1;
            q2_IK = inv(K)*q2;
            q3_IK = inv(K)*q3;
            q4_IK = inv(K)*q4;

            p1_IK = p1_IK/p1_IK(3);
            p2_IK = p2_IK/p2_IK(3);
            p3_IK = p3_IK/p3_IK(3);
            p4_IK = p4_IK/p4_IK(3);

            q1_IK = q1_IK/q1_IK(3);
            q2_IK = q2_IK/q2_IK(3);
            q3_IK = q3_IK/q3_IK(3);
            q4_IK = q4_IK/q4_IK(3);

            pt_IK = [p1_IK,p2_IK,p3_IK,p4_IK]';
            qt_IK = [q1_IK,q2_IK,q3_IK,q4_IK]';

            H_exact_IK = R2-1/d*t2*n';
            A = calcAffinsFromHomography(H_exact_IK,pt_IK,qt_IK);

            AN = A;
            for i = 1:numOfPoints
                A_tmp = [A(i,1) A(i,2); A(i,3) A(i,4)];
                A_noisy = AddNoiseToAffine(A_tmp,noiseScale,noiseAngle);
                AN(i,:) = [A_noisy(1,1) A_noisy(1,2) A_noisy(2,1) A_noisy(2,2)];
            end
            A = AN;


            %4 point method
            H_4pt_IK = calcSolution4pnt(pt_IK,qt_IK);
            H_4pt_IK = H_4pt_IK/H_4pt_IK(2,2);
            [R_est_4pt,t_est_4pt,n_est_4pt,alpha_est_4pt] = decomposeHomography(H_4pt_IK);
            [error_R_4pt,error_t_4pt] = calcAlphaBetaErrors(R_est_4pt,t_est_4pt,R2,t2);
            errors(j,1) = error_R_4pt;
            errors(j,2) = error_t_4pt;

            %2Affin method
            H_2A_IK = calcSolution2Affin(A,pt_IK,qt_IK);
            H_2A_IK = H_2A_IK/H_2A_IK(2,2);


            [R_est_2A,t_est_2A,n_est_2A,alpha_est_2A]=decomposeHomography(H_2A_IK);
            [error_R_2A,error_t_2A] = calcAlphaBetaErrors(R_est_2A,t_est_2A,R2,t2);
            errors(j,3) = error_R_2A;
            errors(j,4) = error_t_2A;

            %Proposed methods
            switch proposedMethod
                case 0
                    H_Prop_IK_rapid = calcSolutionAffinGround(A,pt_IK,qt_IK);
                    H_Prop_IK_rapid = H_Prop_IK_rapid/H_Prop_IK_rapid(2,2);

                    [R_est_Prop_rapid, t_est_Prop_rapid, n_est_Prop_rapid] =  decomposeHomographyGround(H_Prop_IK_rapid);

                    H_Prop_IK_optimal = calcSolutionAffinGround_Optimal(A,pt_IK,qt_IK);
                    H_Prop_IK_optimal = H_Prop_IK_optimal/H_Prop_IK_optimal(2,2);

                    [R_est_Prop_optimal, t_est_Prop_optimal] = decomposeHomographyGround(H_Prop_IK_optimal);


                case 1
                    H_Prop_IK_rapid = calcSolutionAffinWall(A,pt_IK,qt_IK);
                    H_Prop_IK_rapid = H_Prop_IK_rapid/H_Prop_IK_rapid(2,2);
                    [R_est_Prop_rapid, t_est_Prop_rapid, n_est_Prop_rapid] = decomposeHomography(H_Prop_IK_rapid);

                    H_Prop_IK_optimal = calcSolutionAffinFrontWall_Optimal(A,pt_IK,qt_IK);
                    H_Prop_IK_optimal = H_Prop_IK_optimal/H_Prop_IK_optimal(2,2);
                    %[R_est_Prop_optimal, t_est_Prop_optimal] = decomposeHomographyFrontWall(H_Prop_IK_optimal);
                    [R_est_Prop_optimal, t_est_Prop_optimal] = decomposeHomography(H_Prop_IK_optimal);

                case 2
                    H_Prop_IK_rapid = calcSolutionAffinWall(A,pt_IK,qt_IK);
                    H_Prop_IK_rapid = H_Prop_IK_rapid/H_Prop_IK_rapid(2,2);
                    [R_est_Prop_rapid, t_est_Prop_rapid, n_est_Prop_rapid] = decomposeHomography(H_Prop_IK_rapid);

                    H_Prop_IK_optimal = calcSolutionAffinSideWall_Optimal(A,pt_IK,qt_IK);
                    H_Prop_IK_optimal = H_Prop_IK_optimal/H_Prop_IK_optimal(2,2);
                    [R_est_Prop_optimal, t_est_Prop_optimal, n_est_Prop_optimal] = decomposeHomography(H_Prop_IK_optimal);

                case 3
                    H_Prop_IK_rapid = calcSolutionAffinWall(A,pt_IK,qt_IK);
                    H_Prop_IK_rapid = H_Prop_IK_rapid/H_Prop_IK_rapid(2,2);
                    [R_est_Prop_rapid, t_est_Prop_rapid, n_est_Prop_rapid] = decomposeHomography(H_Prop_IK_rapid);

            end

            [error_R_Prop_rapid,error_t_Prop_rapid] = calcAlphaBetaErrors(R_est_Prop_rapid,t_est_Prop_rapid,R2,t2);
            errors(j,5) = error_R_Prop_rapid;
            errors(j,6) = error_t_Prop_rapid;

            if useOptimal == 1
                [error_R_Prop_optimal,error_t_Prop_optimal] = calcAlphaBetaErrors(R_est_Prop_optimal,t_est_Prop_optimal,R2,t2);
                errors(j,13) = error_R_Prop_optimal;
                errors(j,14) = error_t_Prop_optimal;
            end

            if use3pt == 1
                [R_3pt,t_3pt] = linear3pt(pt_IK,qt_IK);
                R_3pt = R_3pt';
                [error_R_3pt,error_t_3pt] = calcAlphaBetaErrors(R_3pt,t_3pt,R2,R2'*t2);
                errors(j,7) = error_R_3pt;
                errors(j,8) = error_t_3pt;
            end

            %Choi Line
            [R_ChoiLine,t_ChoiLine] = choiLine(pt_IK,qt_IK);
            [error_R_ChoiLine,error_t_ChoiLine] = calcAlphaBetaErrors(R_ChoiLine,t_ChoiLine,R2,R2'*t2);
            errors(j,9) = error_R_ChoiLine;
            errors(j,10) = error_t_ChoiLine;

            %Choi Circle
            [R_ChoiCircle,t_ChoiCircle] = choiCircle(pt_IK,qt_IK);
            [error_R_ChoiCircle,error_t_ChoiCircle] = calcAlphaBetaErrors(R_ChoiCircle,t_ChoiCircle,R2,R2'*t2);
            errors(j,11) = error_R_ChoiCircle;
            errors(j,12) = error_t_ChoiCircle;

        end

        for i = 1:numOfMethods
            errors_stat(k,(i-1)*4+1) = mean(errors(:,(i-1)*2+1));
            errors_stat(k,(i-1)*4+2) = var(errors(:,(i-1)*2+1));
            errors_stat(k,(i-1)*4+3) = mean(errors(:,(i-1)*2+2));
            errors_stat(k,(i-1)*4+4) = var(errors(:,(i-1)*2+2));
        end
    end

    
end
   
if plotOnly == 1
    fileName = strcat(destFolder,'testCase',num2str(testCase),'.mat');
    data = matfile(fileName);
    errors_stat = data.errors_stat
else
    save(strcat(destFolder,'testCase',num2str(testCase),'.mat'),'errors_stat');
end
errors_stat(1,17) = 0;
errors_stat(1,19) = 0;


lw = 3;
figure
set(gca,'FontSize',14)
hold on
plot(noiseArray,errors_stat(:,9),'LineWidth',lw,'Color',propColor, 'LineStyle', propLineStyle, 'Marker', propMarker, 'MarkerSize', 12)
hold on
if useOptimal == 1
    plot(noiseArray,errors_stat(:,25),'LineWidth',lw,'Color',propColorO, 'LineStyle',  propLineStyle, 'Marker', propMarkerO, 'MarkerSize', 12)
    hold on
end
if use3pt == 1
    plot(noiseArray,errors_stat(:,13),'LineWidth',lw,'Color',[218/255 165/255 32/255], 'LineStyle', '-.', 'Marker', 's', 'MarkerSize', 12)
    hold on
end
plot(noiseArray,errors_stat(:,17),'LineWidth',lw,'Color',[33 64 154] / 255, 'LineStyle', '-', 'Marker', '>', 'MarkerSize', 12)
hold on
plot(noiseArray,errors_stat(:,21),'LineWidth',lw,'Color',[127.5 0 0] / 255, 'LineStyle', ':', 'Marker', '<', 'MarkerSize', 12)
hold on
plot(noiseArray,errors_stat(:,1),'LineWidth',lw,'Color','b', 'LineStyle', '-.', 'Marker', 'd', 'MarkerSize', 12)
hold on
plot(noiseArray,errors_stat(:,5),'LineWidth',lw,'Color','g', 'LineStyle', '--', 'Marker', '>', 'MarkerSize', 12)
grid()
if (use3pt == 1) && (useOptimal == 1)
    lgnd = legend({propMethodRapid,propMethodOptimal,'3PC Linear','2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 0) && (useOptimal == 1)
    lgnd = legend({propMethodRapid,propMethodOptimal,'2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 1) && (useOptimal == 0)
    lgnd = legend({propMethodRapid,'3PC Linear','2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
else
    lgnd = legend({propMethodRapid,'2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
end
set(lgnd,'color','none');
ylabel('Rotation Error (degrees)','FontSize', 14)
xlabel("Noise",'FontSize', 14) 
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Rot.png'))
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Rot.fig'))
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Rot.svg'))

hold off

lw = 3;
figure
set(gca,'FontSize',14)
hold on
plot(noiseArray,errors_stat(:,11),'LineWidth',lw,'Color',propColor, 'LineStyle',  propLineStyle, 'Marker', propMarker, 'MarkerSize', 12)
hold on
if useOptimal == 1
    plot(noiseArray,errors_stat(:,27),'LineWidth',lw,'Color',propColorO, 'LineStyle',  propLineStyle, 'Marker', propMarkerO, 'MarkerSize', 12)
    hold on
end
if use3pt == 1
    plot(noiseArray,errors_stat(:,15),'LineWidth',lw,'Color',[218/255 165/255 32/255], 'LineStyle', '-.', 'Marker', 's', 'MarkerSize', 12)
    hold on
end
plot(noiseArray,errors_stat(:,19),'LineWidth',lw,'Color',[33 64 154] / 255, 'LineStyle', '-', 'Marker', '>', 'MarkerSize', 12)
hold on
plot(noiseArray,errors_stat(:,23),'LineWidth',lw,'Color',[127.5 0 0] / 255, 'LineStyle', ':', 'Marker', '<', 'MarkerSize', 12)
hold on
plot(noiseArray,errors_stat(:,3),'LineWidth',lw,'Color','b', 'LineStyle', '-.', 'Marker', 'd', 'MarkerSize', 12)
hold on
plot(noiseArray,errors_stat(:,7),'LineWidth',lw,'Color','g', 'LineStyle', '--', 'Marker', '>', 'MarkerSize', 12)
grid()
if (use3pt == 1) && (useOptimal == 1)
   lgnd = legend({propMethodRapid,propMethodOptimal,'3PC Linear','2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 0) && (useOptimal == 1)
    lgnd = legend({propMethodRapid,propMethodOptimal,'2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 1) && (useOptimal == 0)
   lgnd =  legend({propMethodRapid,'3PC Linear','2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
else
    lgnd = legend({propMethodRapid,'2PC Line','2PC Circle','4 PC','2AC'},'FontSize', 12,'Location','northwest');
end
set(lgnd,'color','none');
xlabel("Noise",'FontSize', 14) 
ylabel('Translation Error (degrees)','FontSize', 14)
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Tr.png'))
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Tr.fig'))
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Tr.svg'))
