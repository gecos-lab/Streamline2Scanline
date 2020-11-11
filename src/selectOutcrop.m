function [thisRidgesShp, thisOutcropShp, thisStreamlnsShp, thisIntersctsShp, thisOutcrop, nThisStreamlns, nThisInterscts, thisOutcropPolygon, gridArray, Xgrid, Ygrid, LevelGrid, gradLevelX, gradLevelY] = selectOutcrop(RidgesShp, OutcropsShp, StreamlnsShp, IntersctsShp)
% streamline2scanline: Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia
%
% @ 2020 by Andrea Bistacchi
% distributed under the GNU AGPL v3.0 license.
%
% last updated 30/6/2020

%% Select outcrop to be processed.
disp(' ')
disp('2- Selecting outcrop to be processed')

OutcropsList = [OutcropsShp.OutcropId];
disp(['Valid outcrop IDs are: ' num2str(OutcropsList)])

thisOutcrop = NaN;

while not(thisOutcrop==OutcropsList)
    thisOutcrop = input('Type outcrop ID: ');
end

disp(['Selected outcrop: ' num2str(thisOutcrop)])

% Extract ridges of this outcrop
thisRidgesShp = RidgesShp([RidgesShp.OutcropId] == thisOutcrop);
disp(['N. of ridges in outcrop ' num2str(length(thisRidgesShp))])

% Extract boundary of this outcrop
thisOutcropShp = OutcropsShp([OutcropsShp.OutcropId] == thisOutcrop);
if length(thisOutcropShp)>1
    disp(['ERROR: this file includes more than one polyline with index ' num2str(thisOutcrop)]);
    disp('   -> CHOSE ANOTHER OUTCROP OR CHECK INPUT FILES')
end

%% Create level field for scanline interpolation

% initialize
X = [];
Y = [];
Level = [];

% loop over ridges in outcrop
for j = 1:length(thisRidgesShp)
    % extract X Y data points and assign a LEVEL value to each point
    newX = thisRidgesShp(j).X;
    newY = thisRidgesShp(j).Y;
    newLevel = thisRidgesShp(j).LEVEL*ones(size(thisRidgesShp(j).X));
    
    % check for NaNs
    okDataId = find(~isnan(newX).* ~isnan(newY).* ~isnan(newLevel));
    newX = newX(okDataId);
    newY = newY(okDataId);
    newLevel = newLevel(okDataId);
    
    % add data to arrays
    X = [X newX];
    Y = [Y newY];
    Level = [Level newLevel];
end

% define interpolation area from bounday and ridges
maxX = max(max(X),max(thisOutcropShp.X));
minX = min(min(X),min(thisOutcropShp.X));
maxY = max(max(Y),max(thisOutcropShp.Y));
minY = min(min(Y),min(thisOutcropShp.Y));

% create interpolant object
levelFieldInterpolant = scatteredInterpolant(X',Y',Level','natural');

% create regular grid with twice cells in the "long" direction than input polylines
longAxis = 2 - (maxX - minX) <= (maxY - minY);  % 1 -> X is long axis, 2 -> Y is long axis
switch longAxis
    case 1
        gridSpac = (maxX-minX) / length(thisRidgesShp);
    case 2
        gridSpac = (maxY-minY) / length(thisRidgesShp);
end
gridSpac = min(gridSpac,1400); % added for smaller outcrops

% correct spacing along X and Y in order to have almost the same spacing between max and min
Xspac = (maxX-minX) / round((maxX-minX)/gridSpac);
Yspac = (maxY-minY) / round((maxY-minY)/gridSpac);

% % check
% disp(['longAxis ' num2str(longAxis)])
% disp(['maxX-minX ' num2str(maxX-minX)])
% disp(['maxY-minY ' num2str(maxY-minY)])
% disp(['gridSpac ' num2str(gridSpac)])
% disp(['Xspac ' num2str(Xspac)])
% disp(['Yspac ' num2str(Yspac)])

% grid and interpolate LEVEL
[Xgrid,Ygrid] = meshgrid(minX:Xspac:maxX, minY:Yspac:maxY);
LevelGrid = levelFieldInterpolant(Xgrid,Ygrid);

% interpolate gradient
[gradLevelX,gradLevelY] = gradient(LevelGrid,Xspac,Yspac);

% polygon from outcrop boundary
thisOutcropPolygon = polyshape(thisOutcropShp.X,thisOutcropShp.Y);

% find grid points strictly inside polygon
gridArray = [reshape(Xgrid,[],1) reshape(Ygrid,[],1)];
[isInside,isOnBnd] = isinterior(thisOutcropPolygon,gridArray);
isInside = find(isInside - isOnBnd); % standard isInside includes boundary LEAVE THIS???
gridArray = gridArray(isInside,:);

%% Extract ridges and intersections of this outcrop.

thisStreamlnsShp = StreamlnsShp([StreamlnsShp.OutcropId] == thisOutcrop);
nThisStreamlns = length(thisStreamlnsShp);
disp(['N. of streamlines in outcrop ' num2str(nThisStreamlns)])

thisIntersctsShp = IntersctsShp([IntersctsShp.OutcropId] == thisOutcrop);
nThisInterscts = length(thisIntersctsShp);
disp(['N. of intersections in outcrop ' num2str(nThisInterscts)])


%% Plot outcrop boundary and ridges
fig2 = figure(2);
clf(fig2, 'reset')
fig2.Name = 'Outcrop map';
hold on
title(['Arcadia Planitia Outcrop ' num2str(thisOutcrop)])
axis tight
axis equal
%xlabel('X [m]');
%ylabel('Y [m]');
set(gca,'XTickLabel',[],'YTickLabel',[])
grid on
box on
plot(thisOutcropShp.X,thisOutcropShp.Y,'LineWidth',1.5,'Color',[0,0.7,1]);
for j=1:length(thisRidgesShp)
    plot(thisRidgesShp(j).X,thisRidgesShp(j).Y,'LineWidth',1.5,'Color',[0.7,0,0]);
end
scalebar('Unit','m','Location','southeast') % set this here to preserve correct scale
drawnow
% plot field
%plot(gridArray(:,1),gridArray(:,2),'.')
contour(Xgrid,Ygrid,LevelGrid,'--','Color',[0,0.7,1])
quiver(Xgrid,Ygrid,gradLevelX,gradLevelY,'Color',[0,0.7,1])
drawnow
hold off

% SEE WHAT TO DO HERE - USE gridArray THAT IS TRIMMED TO OUTCROP
% BOUNDARY??? OR ALL OTHER GRIDS THAT ARE NOT TRIMMED?

end

