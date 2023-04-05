% This is a basic version of the QD segregation method which does not
% include any dauphane twin removal and assumes that the only phase
% collected is 'Quartz-new'.  
% 
% Installation of the MTEX toolbox is required, and can be found at:
% https://mtex-toolbox.github.io/
% 
% Knee.m function from Cross et al., 2017 included below. 
% Also available at https://www.researchgate.net/publication/326341597_RexRelict_version_22_Cross_et_al_2017_GRL
% 
%
% file_name = the name of your .ctf file with single quotes, e.g. 'hello_world.ctf'

close all
clear all
file_name = '170815-2.ctf'     % <======= input your file name here
GrainCond = .30;               % <======= Choose your GBF threshold here, e.g. 30% shared boundary = .3 

%% color convention
color_Connected = [106, 166, 219]./255;                                     %selected from adobe color, triad
color_Isolated = [219, 106, 101]./255;

%Trigonal Symmetry --- use hexagonal symmetry to remove dauphane twins
CSTrig = {...
'notIndexed',...
crystalSymmetry('-3m1', [4.913 4.913 5.504], 'X||a*', 'Y||b', 'Z||c',...
'mineral', 'Quartz-new', 'color', 'blue'),...
crystalSymmetry('-3m1', [4.913 4.913 5.504], 'X||a*', 'Y||b', 'Z||c',...
'mineral', 'Connected', 'color', color_Connected),...
crystalSymmetry('-3m1', [4.913 4.913 5.504], 'X||a*', 'Y||b', 'Z||c',...
'mineral', 'Isolated', 'color', color_Isolated)};

%% Setting Plotting Convention - Specific to your EBSD data
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','OutOfPlane');

%% load EBSD data  ---- this section will need to be tailored to your system
ebsd = EBSD.load(file_name,CSTrig,'interface','ctf',...
    'convertEuler2SpatialReferenceFrame');
ebsd = flipud(ebsd);                                                        %this line may need to be removed depending on your SEM setup

%% data cleaning
ebsd(ebsd.mad >0.9)=[];                                                     %remove high MAD value pixels
[grains,ebsd.grainId] =...
    calcGrains(ebsd,'angle',10*degree, 'boundary', 'tight');
ebsd(grains(grains.grainSize< 4)) = [];                                     %remove small grains
%% re-calculate grains
[grains,ebsd.grainId] =...
    calcGrains(ebsd,'angle',10*degree, 'boundary', 'tight');



%% grain segregation - Quartz Domain method
[temp_grains,ebsd.grainId] = calcGrains(ebsd,'angle',360*degree);           %calculate grains using 360 angle
temp_grains=temp_grains('Quartz-new');
knee = tradeOff(temp_grains.grainSize);                                     % Use a trade-off curve to find the cutoff between low and high grainSize values
set(gcf,'Visible','off');                                                   % makes tradeoff curve figure invisible - comment out to see curve
ebsd(temp_grains(temp_grains.grainSize>=knee)).phase=2;                     % reassign large quartz domains to "Connected" phase
ebsd(temp_grains(temp_grains.grainSize<knee)).phase=3;                      % reassign small quartz domain to "Isolated" phase


%% Recalculate grains
[post_grains ebsd.grainId ebsd.mis2mean] =...
     calcGrains(ebsd,'angle',10*degree); 

%% plot Map
figure;
plot(post_grains('Connected'),'linecolor','black','linewidth',.2);
title('GBF Map','Interpreter','none')
hold on
plot(post_grains('Isolated'),'linecolor','black','linewidth',.2);
grains_Map_Attributes_All.micronBar.backgroundAlpha=1; 
legend('Connected', 'Isolated')


%% Calculate ODFs

odf_All = calcDensity(grains('Quartz-new').meanOrientation,...
    'halfwidth', 5*degree);
odf_Connected = calcDensity(post_grains('Connected').meanOrientation,...
    'halfwidth', 5*degree);
odf_Isolated = calcDensity(post_grains('Isolated').meanOrientation,...
    'halfwidth', 5*degree);
h_All = [Miller(0,0,0,1,odf_All.CS),Miller(1,1,-2,0,odf_All.CS)];

%% Plot Pole figures
% pole figure of all grains
fig_PDF_Connected = figure();
plotPDF(odf_All,h_All,'lower')
mtexColorbar('Title','M.U.D.')
caption=['n=', num2str(length(grains('Quartz-new')))];
text(-1.505,-1.23,0,caption,'FontSize',24,'HorizontalAlignment','center');
text(-1.505,1.60,0, ' All Grains','FontSize',24,'HorizontalAlignment',...
    'center','Interpreter','none');

% pole figure of Connected Grains
fig_PDF_Connected = figure();
plotPDF(odf_Connected,h_All,'lower')
mtexColorbar('Title','M.U.D.')
caption=['n=', num2str(length(post_grains('Connected')))];
text(-1.505,-1.23,0,caption,'FontSize',24,'HorizontalAlignment','center');
text(-1.505,1.60,0, 'Connected Grains','FontSize',24,...
    'HorizontalAlignment','center','Interpreter','none');

% pole figure of Isolated Grains
fig_PDF_Isolated = figure();
plotPDF(odf_Isolated,h_All,'lower')
mtexColorbar('Title','M.U.D.')
caption=['n=', num2str(length(post_grains('Isolated')))];
text(-1.505,-1.23,0,caption,'FontSize',24,'HorizontalAlignment','center');
text(-1.505,1.60,0, 'Isolated Grains','FontSize',24,...
    'HorizontalAlignment','center','Interpreter','none');


%% FUNCTIONS

%==========================================================================
%% tradeOff.m - Andrew J. Cross

% Cross, A. J., Kidder, S., & Prior, D. J. (2015). Using microstructures and
% TitaniQ thermobarometry of quartz sheared around garnet porphyroclasts to
% evaluate microstructural evolution and constrain an Alpine Fault Zone
% geotherm. Journal of Structural Geology, 75, 17-31.
%
% A function to find the 'knee' of a trade-off curve. A trade-off curve is used
% to visualise two populations with different magnitudes of a given
% variable. Trade-off curves are ideally shaped like an exponential function 
% - the 'knee' of the curve (the point furthest from a line connecting the 
% first and last points of the curve) provides the critical threshold value
% of the given variable, to separate the two populations. 
%
%
%         ^  
%         |                        o
%         |                        o
%         |                        o
%     y   |                       o 
%         |                      x  <== knee in the curve
%         |                   o
%         |    o    o   o
%         |     
%         +------------------------------>
%              number of measurements
%
% 
% Here, the knee in the curve is found by drawing a line between the first
% and last points in the curve; the point on the curve which is furthest
% from the line is the knee. 


function knee = tradeOff(y) 

    y = sort(y)'; % Variable used to separate two populations
    x = 1:length(y); % Number of measurements

    % Find 'the line' passing through the first and last points 
    % of the trade-off curve
    line = polyfit([x(1) x(end)],[y(1) y(end)],1);

    % Find distance between each point and the curve
    dist = abs((line(1).*x)+(-1.*y)+line(2))/...
            sqrt(line(1)^2+(-1)^2);

    % The 'knee' in the trade-off curve is at the point furthest from 'the line'
    [c id] = max(dist);

    % Find y value at this point
    knee = y(id);

    % Plot trade-off curve and the knee in the curve
    figure, plot(x,y,'linewidth',3)
    hold on
    scatter(x(id),knee,100,'or','linewidth',2)
    plot([0 x(id)],[knee knee],'--k','linewidth',2)

end



