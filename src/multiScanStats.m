function multiScanStats(thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, thisOutcropID)
% streamline2scanline: Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia
%
% @ 2020 by Andrea Bistacchi
% distributed under the GNU AGPL v3.0 license.
%
% last updated 30/6/2020

%% 4. Merged-scanline statistics
% In this section we calculate cumulative spacing statistics of all the scanlines of one outcrop

disp(' ')
disp('5- Multi-scanline statistics')

fig4 = figure(4);
clf(fig4, 'reset')
fig4.Name = 'Multi-scanline statistics';
hold on
sgtitle(['Outcrop: ' num2str(thisOutcropID)],'FontName','Times','FontSize',12,'FontWeight','bold')

%% Extract all distances
% retrieve streamline IDs that after editing can be no more seqeuntial
currStreamlnIds = [thisOutcropStreamlnsShp.StreamlnId];
S = [];
% loop over scanlines
for j=1:nThisOutcropStreamlns
    % extract intersections of this scanline and populate Dist vector of distances
    % note that X was used in 1D scanline App instead of Dist, but here it
    % might create confusion with 2D map coordinates
    currStreamlnId = currStreamlnIds(j);
    thisScanIntersctsShp = thisOutcropIntersctsShp([thisOutcropIntersctsShp.StreamlnId] == currStreamlnId);
    Dist = [thisScanIntersctsShp.ScanDist]; % square brackets to convert to matrix
    % calculate spacing as length and 'baricenter' coordinate of 'bricks' between two fractures
    distS = (Dist(2:end)+Dist(1:end-1))/2; % baricenter coordinate of bricks
    thisS = (Dist(2:end)-Dist(1:end-1)); % length of bricks
    S = [S thisS];
end

%% Statistics

%calculate stats and P10
meanS = mean(S);
medianS = median(S);
modeS = mode(S);
stdS = std(S);
skewS = skewness(S);
kurtS = kurtosis(S);
P10 = 1/meanS;

%_____________________________________________________________________
% empirical cumulative distribution of spacing for all following tests
[observSfreq,Secdf] = ecdf(S);

%____________________________________________  %%%%%%%%%%%%%%%%%%%%%%
% spacing histogram
subplot(2,2,1)
hold on
histogram(S,'Normalization','pdf','FaceColor',[.6 .6 .6]);
box on
grid on
xlabel 'fracture spacing S [m]';
ylabel 'frequency'
title('PDF fit');
plot([meanS meanS],ylim,'k-','LineWidth',1)
plot([medianS medianS],ylim,'k--','LineWidth',1)
plot([modeS modeS],ylim,'k:','LineWidth',1)

%____________________________________________  %%%%%%%%%%%%%%%%%%%%%%
% spacing cumulative
subplot(2,2,3)
hold on
plot(Secdf,observSfreq,'k-','LineWidth',1)
title('CDF fit')
box on
grid on

%________________________________________
% K-S test for EXPONENTIAL spacing distribution

% fit EXPONENTIAL distribution
ExpSdist = fitdist(S','Exponential');

% create expected EXPONENTIAL CDF and PDF
expectExpScdf = cdf(ExpSdist,Secdf);
expectExpSpdf = pdf(ExpSdist,Secdf);

% Lilliefors test
[ExpLHo,ExpLPval,~,~] = lillietest(S,'Distr','exp');

if ExpLHo == 0
    ExpLOutcome = 'Exponential spacing dist. detected at 5% sign.';
else
    ExpLOutcome = 'NO Exponential spacing dist. detected at 5% sign.';
end

% K-S test empirical vs. expected EXPONENTIAL
[ExpKSHo,ExpKSPval,~,~] = kstest(S,'CDF',[Secdf expectExpScdf]);

if ExpKSHo == 0
    ExpKSOutcome = 'Exponential spacing dist. detected at 5% sign.';
else
    ExpKSOutcome = 'NO Exponential spacing dist. detected at 5% sign.';
end

% plotting PDF
subplot(2,2,1)
plot(Secdf,expectExpSpdf,'r')

% plotting CDF
subplot(2,2,3)
plot(Secdf,expectExpScdf,'r')

%____________________________________________
% K-S test for LOGNORMAL spacing distribution

% fit LOGNORMAL distribution
LognormSdist = fitdist(S','Lognormal');

% create expected LOGNORMAL CDF and PDF
expectLognormScdf = cdf(LognormSdist,Secdf);
expectLognormSpdf = pdf(LognormSdist,Secdf);

% K-S test empirical vs. expected LOGNORMAL
[SLognormKSHo,SLognormKSPval,~,~] = kstest(S,'CDF',[Secdf expectLognormScdf]);

if SLognormKSHo == 0
    SLognormKSOutcome = 'Lognormal spacing dist. detected at 5% sign.';
else
    SLognormKSOutcome = 'NO Lognormal spacing dist. detected at 5% sign.';
end

% plotting PDF
subplot(2,2,1)
plot(Secdf,expectLognormSpdf,'b')

% plotting CDF
subplot(2,2,3)
plot(Secdf,expectLognormScdf,'b')

%_________________________________________
% K-S test for NORMAL spacing distribution

% fit NORMAL distribution
NormSdist = fitdist(S','Normal');

% create expected NORMAL CDF and PDF
expectNormScdf = cdf(NormSdist,Secdf);
expectNormSpdf = pdf(NormSdist,Secdf);

% Lilliefors test
[NormLHo,NormLPval,~,~] = lillietest(S,'Distr','norm');

if NormLHo == 0
    NormLOutcome = 'Normal spacing dist. detected at 5% sign.';
else
    NormLOutcome = 'NO Normal spacing dist. detected at 5% sign.';
end

% K-S test empirical vs. expected NORMAL
[NormKSHo,NormKSPval,~,~] = kstest(S,'CDF',[Secdf expectNormScdf]);

if NormKSHo == 0
    NormKSOutcome = 'Normal spacing dist. detected at 5% sign.';
else
    NormKSOutcome = 'NO Normal spacing dist. detected at 5% sign.';
end

% plotting PDF
subplot(2,2,1)
plot(Secdf,expectNormSpdf,'Color',[0 0.5 0])

% plotting CDF
subplot(2,2,3)
plot(Secdf,expectNormScdf,'Color',[0 0.5 0])

%____________________________________________  %%%%%%%%%%%%%%%%%%%%%%
% text

%     disp([' Scanline: ' num2str(currStreamlnId)])
%
%     disp(['Spacing statistics:',...
%         ' - mean = ' num2str(meanS),...
%         ' - median = ' num2str(medianS),...
%         ' - mode = ' num2str(modeS),...
%         ' - std.dev. = ' num2str(stdS),...
%         ' - skewness = ' num2str(skewS),...
%         ' - kurtosis = ' num2str(kurtS),...
%         ' - P10 = ' num2str(P10),...
%         'Tests:',...
%         'Lilliefors Exponential spacing dist. p-value = ' num2str(ExpLPval),...
%         'Lilliefors Exponential spacing dist. outcome = ' ExpLOutcome,...
%         'K-S Exponential spacing dist. p-value = ' num2str(ExpKSPval),...
%         'K-S Exponential spacing dist. outcome = ' ExpKSOutcome,...
%         'K-S Lognormal spacing dist. p-value = ' num2str(SLognormKSPval),...
%         'K-S Lognormal spacing dist. outcome = ' SLognormKSOutcome,...
%         'Lilliefors Normal spacing dist. p-value = ' num2str(NormLPval),...
%         'Lilliefors Normal spacing dist. p-value = ' NormLOutcome,...
%         'K-S Normal spacing dist. p-value = ' num2str(NormKSPval),...
%         'Lilliefors Normal spacing dist. p-value = ' NormKSOutcome]);

% create table
Statistics = {...
    'mean';...
    'median';...
    'mode';...
    'std.dev.';...
    'skewness';...
    'kurtosis';...
    'P10';...
    'Lilliefors Exp. S';...
    'K-S Exp. S';...
    'K-S Lognorm. S';...
    'Lilliefors Norm. S';...
    'K-S Normal S'};

Values = [...
    meanS;
    medianS;
    modeS;
    stdS;
    skewS;
    kurtS;
    P10;
    ExpLPval;
    ExpKSPval;
    SLognormKSPval;
    NormLPval;
    NormKSPval];

Outcomes = {...
    ' ';...
    ' ';...
    ' ';...
    ' ';...
    ' ';...
    ' ';...
    ' ';...
    ExpLOutcome;...
    ExpKSOutcome;...
    SLognormKSOutcome;...
    NormLOutcome;...
    NormKSOutcome};

Table = table(Values, Outcomes, 'RowNames',Statistics);

% Get the table in string form
TString = evalc('disp(Table)');

% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');

% Get a fixed-width font
FixedWidth = get(0,'FixedWidthFontName');

% show table in subplot
subplot(2,2,[2;4])
axis off
text(0,1,TString,'Interpreter','Tex','FontName',FixedWidth,'Units','Normalized','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',10);


end