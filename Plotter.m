
numOfPoints = 4;

%testCase 0 : Ground plane
%testCase 1 : Wall in the fron
%testCase 2 : Wall on the side
%testCase 3 : Purely forward movement, random wall in the front - 1 fig
%testCase 4 : Purely forward movement, ground plane - 1 fig
%testCase 5 : Ground plane, gravity - 1 fig
%testCase 6 : Random wall in the front, gravity -1 fig

testCase = 6;

final = 0;
if final == 0
    destFolder = 'test\\';
else
    destFolder = 'results\\';
end
purelyForwardMovement = 0;

numOfMethods = 7;
errors = zeros(iterNum,numOfMethods*2);

numOfSteps = 11;
maxNoiseInPixels = 1;
maxNoiseScale = 1;
maxNoiseAngle = 1/180*pi;
noiseArray = zeros(numOfSteps,1);
maxGravityRotAngle = 5;

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
        
end
   
fileName = strcat(destFolder,'testCase',num2str(testCase),'.mat');
data = matfile(fileName);
errors_stat = data.errors_stat

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
    lgnd = legend({propMethodRapid,propMethodOptimal,'3PC Linear','2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 0) && (useOptimal == 1)
    lgnd = legend({propMethodRapid,propMethodOptimal,'2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 1) && (useOptimal == 0)
    lgnd = legend({propMethodRapid,'3PC Linear','2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
else
    lgnd = legend({propMethodRapid,'2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
end
set(lgnd,'color','none');
xlabel("Noise",'FontSize', 14) 
ylabel('Rotation Error (degrees)','FontSize', 14)
switch testCase
    case 0
        ylim([0,6]);
    case 1 
        ylim([0,6])
    case 2
        ylim([0 4])
    case 3
        ylim([0 1])
    case 4
        ylim([0,3.5])
    case 5 
        ylim([0,25])
    case 6
        ylim([0,10])
end
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
   lgnd = legend({propMethodRapid,propMethodOptimal,'3PC Linear','2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 0) && (useOptimal == 1)
    lgnd = legend({propMethodRapid,propMethodOptimal,'2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
elseif (use3pt == 1) && (useOptimal == 0)
   lgnd =  legend({propMethodRapid,'3PC Linear','2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
else
    lgnd = legend({propMethodRapid,'2PC Line','2PC Circle','4PC','2AC'},'FontSize', 12,'Location','northwest');
end
set(lgnd,'color','none');
switch testCase
    case 0
        ylim([0,50]);
    case 1 
        ylim([0,30]);
    case 2
        ylim([0,30]);
    case 3
        ylim([0,7]);
    case 4
        ylim([0,15]);
    case 5 
        ylim([0,80]);
    case 6
        ylim([0,80]);
end

xlabel("Noise",'FontSize', 14) 
ylabel('Translation Error (degrees)','FontSize', 14)
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Tr.png'))
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Tr.fig'))
saveas(gcf,strcat(destFolder,'testCase',num2str(testCase),'_Tr.svg'))
