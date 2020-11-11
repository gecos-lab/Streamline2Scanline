function [thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, nThisOutcropInterscts] = singleScanStats(thisOutcropStreamlnsShp, thisOutcropIntersctsShp, nThisOutcropStreamlns, nThisOutcropInterscts, thisOutcropID, MedianFilterValue)
% streamline2scanline: Spacing of arcuate ridges from thumbprint terranes - Arcadia Planitia
%
% @ 2020 by Andrea Bistacchi
% distributed under the GNU AGPL v3.0 license.
%
% last updated 30/6/2020

%% Single-scanline statistics
% Functions from DomStudioFracStat1D version 2019/08/02

disp(' ')
disp('4- Reviewing single-scanline stats (and deleting bad scanlines)')

fig3 = figure(3);
clf(fig3, 'reset')
fig3.Name = 'Single scanline statistics';
hold on

% use scrollsubplot as follows, for a grid with nThisOutcropStreamlns rows (one for each scanline) and 8 columns:
% column 1: barcode and scan name as title      > scrollsubplot(nThisOutcropStreamlns, 8, j*8-7)
% column 2: cumulative spacing function (CSF)	> scrollsubplot(nThisOutcropStreamlns, 8, j*8-6)
% column 3: cumulative spacing deruvative (CSD) > scrollsubplot(nThisOutcropStreamlns, 8, j*8-5)
% column 4: Poisson distribution test           > scrollsubplot(nThisOutcropStreamlns, 8, j*8-4)
% column 5: histogram and PDF fitting           > scrollsubplot(nThisOutcropStreamlns, 8, j*8-3)
% column 6: ECDF and CDF fitting                > scrollsubplot(nThisOutcropStreamlns, 8, j*8-2)
% column 7-8: text with summary about tests       > scrollsubplot(nThisOutcropStreamlns, 8, j*8-1)

% retrieve streamline IDs that after editing can be no more seqeuntial
currStreamlnIds = [thisOutcropStreamlnsShp.StreamlnId];

% loop over scanlines
for j=1:nThisOutcropStreamlns
    %____________________________________________
    % calculate stats
    
    % extract intersections of this scanline and populate Dist vector of distances
    % note that X was used in 1D scanline App instead of Dist, but here it
    % might create confusion with 2D map coordinates
    currStreamlnId = currStreamlnIds(j);
    thisScanIntersctsShp = thisOutcropIntersctsShp([thisOutcropIntersctsShp.StreamlnId] == currStreamlnId);
    Dist = [thisScanIntersctsShp.ScanDist]; % square brackets to convert to matrix
    
    % scaline length
    L = Dist(end);
    
    % calculate spacing as length and 'baricenter' coordinate of 'bricks' between two fractures
    distS = (Dist(2:end)+Dist(1:end-1))/2; % baricenter coordinate of bricks
    S = (Dist(2:end)-Dist(1:end-1)); % length of bricks
    
    % rank order of S NOT NEEDED SINCE corr(... 'Spearman') IS USED BELOW
    % [~,ii]=sort(S,'Ascend');
    % [~,rankS]=sort(ii); clear ii
    
    %calculate stats and P10
    meanS = mean(S);
    medianS = median(S);
    modeS = mode(S);
    stdS = std(S);
    skewS = skewness(S);
    kurtS = kurtosis(S);
    P10 = 1/meanS;
    
    % Spearman rank correlation test for TREND
    % Ho = no correlation or no TREND
    % Ha = positive or negative correlation or TREND
    % calculate correlation coefficient - xS is baricenter coordinate of bricks
    [~,TrendPval] = corr(distS',S','Type','Spearman'); % data are converted to column vectors as required by corr function
    if TrendPval < 0.05
        TrendOutcome = 'TREND detected at 5% sign.';
    else
        TrendOutcome = 'NO TREND detected at 5% sign.';
    end
    
    % Spearman rank correlation test for PATTERN
    % Ho = no correlation or no PATTERN
    % Ha = positive or negative correlation or PATTERN
    % calculate correlation coefficient
    [~,PatternPval] = corr(S(1:end-1)',S(2:end)','Type','Spearman'); % data are converted to column vectors as required by corr function
    if PatternPval < 0.05
        PatternOutcome = 'PATTERN detected at 5% sign.';
    else
        PatternOutcome = 'NO PATTERN detected at 5% sign.';
    end
    
    %____________________________________________
    % 'barcode' plot
    scrollsubplot(3, 8, j*8-7)
    
    % plot X as a row vector
    barPlotX = Dist;
    if size(barPlotX,1)>1, barPlotX = barPlotX'; end
    
    % aspect ratio of plot
    AR = L/10;
    
    % plot
    plot([barPlotX; barPlotX],[ones(size(barPlotX)); -ones(size(barPlotX))]*AR,'-k','LineWidth',1)
    xlim([0 L]);
    xlabel 'X [m]'
    axis tight
    axis equal
    title(['Outcrop: ' num2str(thisOutcropID) ', Scanline: ' num2str(currStreamlnId)])
    set(gca,'YTick',[])
    box on
    
    %___________________________________________________
    % K-S test for UNIFORM fracture spatial distribution with
    % Cumulative Spacing Function (CSF) and Cumulative Spacing Derivative (CSD)
    
    % create empirical cumulative distribution
    [observDistFreq,DistEcdf] = ecdf(Dist);
    
    % make UNIFORM distribution
    UnifDistDistrb = makedist('Uniform','lower',Dist(1),'upper',Dist(end));
    
    % create expected NORMAL CDF and PDF
    expectDistCdf = cdf(UnifDistDistrb,DistEcdf);
    
    % K-S test empirical vs. expected exponentital
    [UnifDistKsHo,UnifDistKsPval,~,~] = kstest(Dist,'CDF',[DistEcdf expectDistCdf]);
    
    if UnifDistKsHo == 0
        UnifDistKsOutcome = 'UNIFORM spatial distr. detected at 5% sign.';
    else
        UnifDistKsOutcome = 'NO UNIFORM spatial distr. detected at 5% sign.';
    end
    
    % derivatives of observXfreq and expectXcdf
    DobservDistFreq = diff(observDistFreq)./diff(DistEcdf);
    DexpectDistCdf = diff(expectDistCdf)./diff(DistEcdf);
    
    % median filter of derivative
    MedianDobservDistFreq = medfilt1(DobservDistFreq, MedianFilterValue);
    
    % cumulative spacing function (CSF)
    scrollsubplot(3, 8, j*8-6)
    hold on
    plot(DistEcdf, observDistFreq,'k-')
    plot(DistEcdf, expectDistCdf,'r-')
    %legend('Empirical CDF', 'UNIFORM CDF', 'Location','SE')
    xlabel 'Position [m]'
    ylabel 'Proportion of position'
    title('CSF');
    box on
    grid on
    
    % cumulative spacing derivative (CSD)
    scrollsubplot(3, 8, j*8-5)
    hold on
    plot(DistEcdf(2:end), DobservDistFreq,'k-')
    plot(DistEcdf(2:end), DexpectDistCdf,'r-')
    plot(DistEcdf(2:end), MedianDobservDistFreq,'Color',[.3 .3 .3],'LineWidth',1)
    %legend('Empirical CDF derivative', 'UNIFORM CDF derivative', 'Median filter of derivative', 'Location','NE')
    xlabel 'Position [m]';
    ylabel 'Proportion of position'
    title('CSD');
    box on
    grid on
    
    %___________________________________________________
    % Chi-squared test for POISSON fracture distribution
    
    % Poisson parameters
    PoissonT = floor(length(Dist)/5); % number of scanline segments for POISSON distribution counts
    PoissonSegmBnd = linspace(0,L,PoissonT+1);
    
    % assign fractures to segments
    for k=1:PoissonT
        PoissonJ(k) = length(find(and(Dist>=PoissonSegmBnd(k),Dist<PoissonSegmBnd(k+1))));
    end
    
    % create expected distribution from:
    % from http://faculty.business.utsa.edu/manderso//MATLAB/SalmonDigging/SalmonDigging.html
    % not using standard matlab
    
    % descriptive stats (observed counts in segments)
    PoissonObsFreqTable = tabulate(PoissonJ);
    PoissonLambda = mean(PoissonJ);
    PoissonN = length(PoissonJ);
    PoissonBins = PoissonObsFreqTable(:,1);
    
    % calculate expected counts according on POISSON distribution
    PoissonExpectPDF = pdf('Poisson', PoissonObsFreqTable(:,1), PoissonLambda);
    PoissonExpFreq = PoissonN*PoissonExpectPDF;
    
    % GOF test 1
    [PoissHo,PoissPval,~] = chi2gof(PoissonBins, 'Ctrs', PoissonBins, 'Frequency', PoissonObsFreqTable(:,2), 'Expected', PoissonExpFreq, 'NParams', 1, 'EMin',2);
    
    if PoissHo == 0
        PoissOutcome = 'Poisson spatial dist. detected at 5% sign.';
    else
        PoissOutcome = 'NO Poisson spatial dist. detected at 5% sign.';
    end
    
    % plotting standard binning
    scrollsubplot(3, 8, j*8-4)
    bbar = bar(PoissonBins, [PoissonObsFreqTable(:,2), PoissonExpFreq], 'grouped');
    bbar(1).FaceColor = [.6 .6 .6];
    bbar(2).FaceColor = [1 .2 0];
    %legend('Observed', 'Expected', 'Location', 'NorthWest')
    title('Poisson')
    grid on
    
    %     % plotting pooled binning
    %     for j=1:length(PoissStats.O)
    %         PoissonPoolBins{j} = [num2str(PoissStats.edges(j)) '-' num2str(PoissStats.edges(j+1))];
    %     end
    %     PoissonPoolBins = categorical({PoissonPoolBins{1:end}});
    %     bar(PoissonPoolBins, [PoissStats.O', PoissStats.E'], 'grouped')
    %     %legend('Observed', 'Expected', 'Location', 'NorthWest')
    %     title('Poisson')
    
    %_____________________________________________________________________
    % empirical cumulative distribution of spacing for all following tests
    [observSfreq,Secdf] = ecdf(S);
    
    %____________________________________________  %%%%%%%%%%%%%%%%%%%%%%
    % spacing histogram
    scrollsubplot(3, 8, j*8-3)
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
    scrollsubplot(3, 8, j*8-2)
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
    scrollsubplot(3, 8, j*8-3)
    plot(Secdf,expectExpSpdf,'r')
    
    % plotting CDF
    scrollsubplot(3, 8, j*8-2)
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
    scrollsubplot(3, 8, j*8-3)
    plot(Secdf,expectLognormSpdf,'b')
    
    % plotting CDF
    scrollsubplot(3, 8, j*8-2)
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
    scrollsubplot(3, 8, j*8-3)
    plot(Secdf,expectNormSpdf,'Color',[0 0.5 0])
    
    % plotting CDF
    scrollsubplot(3, 8, j*8-2)
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
%         'Spearman Trend p-value = ' num2str(TrendPval),...
%         'Spearman Trend outcome = ' TrendOutcome,...
%         'Spearman Pattern p-value = ' num2str(PatternPval),...
%         'Spearman Pattern outcome = ' PatternOutcome,...
%         'K-S Uniform spatial dist. p-value = ' num2str(UnifDistKsPval),...
%         'K-S Uniform spatial dist.n outcome = ' UnifDistKsOutcome,...
%         'Chi-squared Poisson spatial dist. p-value = ' num2str(PoissPval),...
%         'Chi-squared Poisson spatial dist. outcome = ' PoissOutcome,...
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
        'Spearman Trend';...
        'Spearman Pattern';...
        'K-S Uniform';...
        'Chi-squared Poisson';...
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
        TrendPval;
        PatternPval;
        UnifDistKsPval;
        PoissPval;
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
        TrendOutcome;...
        PatternOutcome;...
        UnifDistKsOutcome;...
        PoissOutcome;...
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
    scrollsubplot(3, 8, j*8-1)
    axis off
    text(0,1,TString,'Interpreter','Tex','FontName',FixedWidth,'Units','Normalized','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',6);
        
end

%% Delete scanlines

while 1
    
    disp(' ');
    disp('__________');
    disp('Delete scanline (type its number, 0 to exit):');
    
    scanToDel = -1;
    while 1
        scanToDel = input(' > ');
        if isnumeric(scanToDel)
            scanToDel = round(scanToDel);
            break
        end
    end
    
    if scanToDel == 0
        break;  % BREAK
    elseif any([thisOutcropStreamlnsShp.StreamlnId] == scanToDel)
        nThisOutcropInterscts = nThisOutcropInterscts - thisOutcropStreamlnsShp([thisOutcropStreamlnsShp.StreamlnId] == scanToDel).IntersctsNum;
        nThisOutcropStreamlns = nThisOutcropStreamlns - 1;
        thisOutcropStreamlnsShp([thisOutcropStreamlnsShp.StreamlnId] == scanToDel) = [];
        thisOutcropIntersctsShp([thisOutcropIntersctsShp.StreamlnId] == scanToDel) = [];
    else
        disp('Scanline not found')
    end
end

%% summary of stats session
disp(['N. of streamlines: ' num2str(nThisOutcropStreamlns)])
disp(['N. of intersections: ' num2str(nThisOutcropInterscts)])

% hold off for outcrop figure
hold off

end

