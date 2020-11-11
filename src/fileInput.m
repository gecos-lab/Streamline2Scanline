function [RidgesShp, OutcropsShp, StreamlnsShp, IntersctsShp, nRidges, nOutcrops, nStreamlns, nInterscts] = fileInput
% streamline2scanline: Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia
%
% @ 2020 by Andrea Bistacchi
% distributed under the GNU AGPL v3.0 license.
%
% last updated 30/6/2020

%% Import files.

disp(' ')
disp('1- Importing input files')

% Get file names
[filename, pathname] = uigetfile('.shp','Select SHP with RIDGES');
inFileRidges = fullfile(pathname,filename);
[filename, pathname] = uigetfile('.shp','Select SHP with OUTCROPS');
inFileOutcrops = fullfile(pathname,filename);
clear pathname filename

% Import shapefiles with arcuate ridges and outcrops as structures
try
    RidgesShp = shaperead(inFileRidges);
catch
    disp('Unexpected file format')
end

try
    OutcropsShp = shaperead(inFileOutcrops);
catch
    disp('Unexpected file format')
end

% rename outcrop id fields (if needed)
try
    RidgesShp = rnfield(RidgesShp,'OUTCROPID','OutcropId');
end

try
    OutcropsShp = rnfield(OutcropsShp,'OUTCROPID','OutcropId');
end

% number of ridge objects (polylines) in shapefile
nRidges = length(RidgesShp);
disp(['N. of ridge polylines = ' num2str(nRidges)])

% number of outcrops
nOutcrops = length(OutcropsShp);
disp(['N. of outcrops = ' num2str(nOutcrops)])

%% Plot overview map.

figure(1)  % 'Overview map'
clf(1)
hold on
title({'Thumbprint terranes of Arcadia Planitia' ; 'mapped in continuous outcrops'})
axis tight
axis equal
%xlabel('X [m]');
%ylabel('Y [m]');
set(gca,'XTickLabel',[],'YTickLabel',[])
grid on
box on
for i=1:nOutcrops
    plot(OutcropsShp(i).X,OutcropsShp(i).Y,'LineWidth',1,'Color',[0,0.7,1])
    text(mean(OutcropsShp(i).X(1:end-1)),mean(OutcropsShp(i).Y(1:end-1)),num2str(OutcropsShp(i).OutcropId),'HorizontalAlignment','center','FontWeight','bold','Color','b','BackgroundColor','w')
end
for i=1:nRidges
    plot(RidgesShp(i).X,RidgesShp(i).Y,'LineWidth',1,'Color',[0.7,0,0]);
end
scalebar('Unit','m','Location','southeast') % set this here to preserve correct scale
hold off

disp('File import completed.')

%% Create empty structures to store project data.

StreamlnsShp = struct('Geometry',{},'BoundingBox',{},'X',{},'Y',{},'OutcropId',{},'StreamlnId',{},'IntersctsNum',{});
IntersctsShp = struct('Geometry',{},'X',{},'Y',{},'OutcropId',{},'StreamlnId',{},'IntersctId',{},'ScanDist',{});
nStreamlns = 0;
nInterscts = 0;

disp('Empty project variables created')

end

