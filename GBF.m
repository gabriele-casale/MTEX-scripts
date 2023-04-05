% This is a basic version of the GBF segregation method which does not
% include any dauphane twin removal and assumes that the only phase
% collected is 'Quartz-new'.  
% 
% Installation of the MTEX toolbox is required, and can be found at:
% https://mtex-toolbox.github.io/
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

%% Grain Segregation - GBF method
grains = grains('Quartz-new');
GBfrac = [];
for i = 1:length(grains)                                                    %loops through grains to calc grain boundary fraction - computationally expensive
GBfrac(i) = sum(grains(i).boundary('Quartz-new','Quartz-new').segLength)...
    ./grains(i).perimeter; 
end
GBF_Pass=find(GBfrac>=GrainCond);                                           %grains that pass condition
GBF_Fail=find(GBfrac<GrainCond);                                            %grains that fail condition
ebsd(grains(GBF_Pass)).phase=2;                                             %reassign connected grains phase
ebsd(grains(GBF_Fail)).phase=3;                                             %reassign isolate grains phase

%% Recalculate grains
post_grains = calcGrains(ebsd,'angle',10*degree); 

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


