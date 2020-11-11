function [thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, nThisOutcropInterscts] = digitizeScanlines(thisOutcropRidgesShp, thisOutcropBoundaryShp, thisOutcropStreamlnsShp, thisOutcropIntersctsShp, thisOutcropID, nThisOutcropStreamlns, nThisOutcropInterscts, thisOutcropPolygon, thisOutcropXgrid, thisOutcropYgrid, thisOutcropLevelGrid, thisOutcropGradLevelX, thisOutcropGradLevelY)
% streamline2scanline: Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia
%
% @ 2020 by Andrea Bistacchi
% distributed under the GNU AGPL v3.0 license.
%
% last updated 30/6/2020

%% Create streamlines starting from digitized seed points (X,Y)

disp(' ')
disp('3- Digitizing scanlines from streamlines')

%% (Re-)Plot outcrop boundary and ridges
fig2 = figure(2);
clf(fig2, 'reset')
fig2.Name = 'Outcrop map';
hold on
title(['Arcadia Planitia Outcrop ' num2str(thisOutcropID)])
axis tight
axis equal
%xlabel('X [m]');
%ylabel('Y [m]');
set(gca,'XTickLabel',[],'YTickLabel',[])
grid on
box on
plot(thisOutcropBoundaryShp.X,thisOutcropBoundaryShp.Y,'LineWidth',1.5,'Color',[0,0.7,1]);
for j=1:length(thisOutcropRidgesShp)
    plot(thisOutcropRidgesShp(j).X,thisOutcropRidgesShp(j).Y,'LineWidth',1.5,'Color',[0.7,0,0]);
end
scalebar('Unit','m','Location','southeast') % set this here to preserve correct scale
drawnow
% plot field
%plot(gridArray(:,1),gridArray(:,2),'.')  % this is trimmed to the outcrop boundary, the others not
contour(thisOutcropXgrid,thisOutcropYgrid,thisOutcropLevelGrid,'--','Color',[0,0.7,1])
quiver(thisOutcropXgrid,thisOutcropYgrid,thisOutcropGradLevelX,thisOutcropGradLevelY,'Color',[0,0.7,1])
drawnow

% plot previous scnlines and intersections
for i = 1:nThisOutcropStreamlns
    plot([thisOutcropStreamlnsShp(i).X],[thisOutcropStreamlnsShp(i).Y],'Color',[0,0,0],'LineStyle','--');
end
plot([thisOutcropIntersctsShp.X],[thisOutcropIntersctsShp.Y],'ok')

%% Digitize scanlines
while 1
    digitAnswer = questdlg('Digitize a streamline?', 'Streamline', 'YES', 'NO', 'YES');
    if strncmp(digitAnswer, 'NO',2)
        % terminate digitizing
        break
    else
        % digitize seed point
        x = NaN;
        y = NaN;
        button = 0;
        while not(isinterior(thisOutcropPolygon,x,y))
            while  button~=1
                [x,y,button] = ginput(1);
                if isempty(button), button = 0; end
            end
        end
        
        % streamline positive side
        XYstreamPos = stream2(thisOutcropXgrid,thisOutcropYgrid,thisOutcropGradLevelX,thisOutcropGradLevelY,x,y);
        XYstreamPosMat = [XYstreamPos{:}];
        
        % streamline negative side
        XYstreamNeg = stream2(thisOutcropXgrid,thisOutcropYgrid,-thisOutcropGradLevelX,-thisOutcropGradLevelY,x,y);
        XYstreamNegMat = flipud([XYstreamNeg{:}]);
        
        % combine negative and positive sides
        if size(XYstreamPosMat,1)>=2 && size(XYstreamNegMat,1)>0
            XYstreamNegMat = XYstreamNegMat(1:end-1,:);
        end
        XYstream = [XYstreamNegMat ; XYstreamPosMat];        
        
        % delete streamline out of outcrop
        [XYstream,~] = intersect(thisOutcropPolygon,XYstream);
        
        % find node closer to seed point
        [~,closeToSeed] = min(sqrt((XYstream(:,1)-x).^2 + (XYstream(:,2)-y).^2));
        
        % select first and last node of segment close to seed point (needed in case of multiparts)
        nanPoints = find(isnan(XYstream(1:end-1,1)));  % end-1 to excude the last closing NaN
        if nanPoints
            firstPoints = [1 ; nanPoints+1];
            [~,closestIndex] = min(abs(closeToSeed - firstPoints));
            firstPoint = firstPoints(closestIndex);
            
            lastPoints = [nanPoints-1 ; size(XYstream,1)-1];
            [~,closestIndex] = min(abs(closeToSeed - lastPoints));
            lastPoint = lastPoints(closestIndex)+1;  % +1 to add the last closing NaN
        else
            firstPoint = 1;
            lastPoint = length(XYstream(:,1));
        end
        
        % fix XYstream
        XYstream = XYstream(firstPoint:lastPoint,:);
        
        % plot if streamline is not empty
        if not(isempty(XYstream))            
            % plot
            lastLine = plot(XYstream(:,1),XYstream(:,2),'Color',[0,0,0],'LineStyle','--');
            
            keepAnswer = questdlg('Keep this streamline for analysis?', 'Streamline', 'YES', 'NO', 'YES');
            if strncmp(keepAnswer, 'YES',3)
                % update index
                if nThisOutcropStreamlns == 0
                    streamId = 1;
                else
                    streamId = max([thisOutcropStreamlnsShp.StreamlnId])+1;
                end
                nThisOutcropStreamlns = nThisOutcropStreamlns + 1;
                
                % copy streamline in structure
                thisOutcropStreamlnsShp(end+1) = struct(...
                    'Geometry','Line',...
                    'BoundingBox',[min(XYstream(:,1)) min(XYstream(:,2)); max(XYstream(:,1)) max(XYstream(:,2))],...
                    'X',XYstream(:,1)',...
                    'Y',XYstream(:,2)',...
                    'OutcropId',thisOutcropID,...
                    'StreamlnId',streamId,...
                    'IntersctsNum',0);
                
                % extract coordinates from SHP
                ThisStream = thisOutcropStreamlnsShp([thisOutcropStreamlnsShp.StreamlnId] == streamId);
                
                % find intersections - loop over ridges in outcrop and find intersection with this streamline
                intThisStrm = [];
                for k = 1:length(thisOutcropRidgesShp)
                    [tempOutX,tempOutY,tempOutID] = polyxpoly(ThisStream.X,ThisStream.Y,thisOutcropRidgesShp(k).X,thisOutcropRidgesShp(k).Y);
                    intThisStrm = [intThisStrm ; [tempOutX tempOutY tempOutID(:,1)]]; % we are interested in indexes of segments along the streamline only
                end
                
                if ~isempty(intThisStrm)
                    % clean and sort intersections along the streamline
                    intThisStrm = sortrows(intThisStrm,3);
                    
                    % calculate spacing
                    spacThisStrm = sqrt((intThisStrm(2:end,1)-intThisStrm(1:end-1,1)).^2 + (intThisStrm(2:end,2)-intThisStrm(1:end-1,2)).^2);
                    spacThisStrm = [spacThisStrm ; NaN];
                    
                    % calculate incrmental distance by adding spacings
                    distThisStrm(1) = 0;
                    for l = 2:length(spacThisStrm)
                        distThisStrm(l) = distThisStrm(l-1) + spacThisStrm(l-1);
                    end
                    
                    % save streamline intersections in structure
                    for k = 1:size(intThisStrm,1)
                        thisOutcropIntersctsShp(end+1) = struct('Geometry','Point',...
                            'X',intThisStrm(k,1),...
                            'Y',intThisStrm(k,2),...
                            'OutcropId',thisOutcropID,...
                            'StreamlnId',streamId,...
                            'IntersctId',intThisStrm(k,3),...
                            'ScanDist',distThisStrm(k));
                        nThisOutcropInterscts = nThisOutcropInterscts+1;
                    end
                    
                    % update number of intersections in scanline
                    thisOutcropStreamlnsShp([thisOutcropStreamlnsShp.StreamlnId] == streamId).IntersctsNum = size(intThisStrm,1);
                    
                    % plot
                    plot(intThisStrm(:,1),intThisStrm(:,2),'ok')
                end
            else
                delete(lastLine)
            end
        else
            disp('streamline is empty')
        end
    end
end

% summary of digitizing session
disp(['N. of streamlines: ' num2str(nThisOutcropStreamlns)])
disp(['N. of intersections: ' num2str(nThisOutcropInterscts)])

% hold off for outcrop figure
hold off


end

