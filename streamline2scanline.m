function streamline2scanline
% streamline2scanline: Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia
%
% @ 2020 by Andrea Bistacchi
% distributed under the GNU AGPL v3.0 license.
%
% last updated 30/6/2020

% Clear and clean
clear all; close all; clc;

% Find path of this function and add sub-folders
addpath(genpath(fileparts(which(mfilename))));

% Initialize figures etc.
set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultAxesFontName','Times','DefaultAxesFontSize',10,'DefaultFigureColor','w');

fig1 = figure(1); fig1.Name = 'Overview map';
fig2 = figure(2); fig2.Name = 'Outcrop map';
fig3 = figure(3); fig3.Name = 'Single scanline statistics';
fig4 = figure(4); fig4.Name = 'Multi-scanline statistics';

% Focus on command window
drawnow; commandwindow

% Define default value of median filter of derivative
MedianFilterValue = 10;

% Welcome message
disp(' ');
disp('Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia - Andrea Bistacchi 30/6/2020');

% Root menu

while 1
    
    disp(' ');
    disp('__________');
    disp('ROOT MENU:');
    disp('1- Clear current data and import input files');
    disp('2- Select outcrop to be processed');
    disp('3- Digitize scanlines from streamlines');
    disp('4- Review single-scanline stats (and delete bad scanlines)');
    disp('5- Multi-scanline statistics');
    disp('6- Save project');
    disp('7- Load previously saved project');
    disp('8- Quit');
    disp(' ');
    
    action = 0;
    while (action < 1) || (action > 9)
        action = input(' > ');
        if isempty(action), action = 0; end
        if not(isnumeric(action)), action = 0; end
        action = round(action);
    end
    
    if     action == 1,	[RidgesShp, OutcropsShp, StreamlnsShp, IntersctsShp, nRidges, nOutcrops, nStreamlns, nInterscts] = fileInput;
        
    elseif action == 2,	[thisOutcropRidgesShp, thisOutcropBoundaryShp, thisOutcropStreamlnsShp, thisOutcropIntersctsShp, thisOutcropID, nThisOutcropStreamlns, nThisOutcropInterscts, thisOutcropPolygon, thisOutcropGridArray, thisOutcropXgrid, thisOutcropYgrid, thisOutcropLevelGrid, thisOutcropGradLevelX, thisOutcropGradLevelY] = selectOutcrop(RidgesShp, OutcropsShp, StreamlnsShp, IntersctsShp);
        
    elseif action == 3,	[thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, nThisOutcropInterscts] = digitizeScanlines(thisOutcropRidgesShp, thisOutcropBoundaryShp, thisOutcropStreamlnsShp, thisOutcropIntersctsShp, thisOutcropID, nThisOutcropStreamlns, nThisOutcropInterscts, thisOutcropPolygon, thisOutcropXgrid, thisOutcropYgrid, thisOutcropLevelGrid, thisOutcropGradLevelX, thisOutcropGradLevelY);
        
    elseif action == 4,	[thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, nThisOutcropInterscts] = singleScanStats(thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, nThisOutcropInterscts, thisOutcropID, MedianFilterValue);
        
    elseif action == 5,	multiScanStats(thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, thisOutcropID);
        
    elseif action == 6
        % save project as mat file
        disp(' ');
        disp('6- Save project to .mat file');
        [file, path] = uiputfile('*.mat');
        save([path file]);
        disp(' ');
        disp([' -> file ' file ' successfully saved.']);
        
    elseif action == 7
        % load project from mat file
        disp(' ');
        disp('7- Load previously saved project from .mat file');
        [file, path] = uigetfile('*.mat');
        load([path file]);
        disp(' ');
        disp([' -> file ' file ' successfully loaded.']);
        
    elseif action == 8
        % quit
        break
    end
end

end