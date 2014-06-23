% TestRansac runs the ransac solver for different constelaltions of outliers and n´missing data. The
% result is thereafter presented in histograms of relative errors for 1. raw outdata, 2. after Gauss-Newton oprimization and 3. for a comparative L1-optimizer with randomized initialization.
%
%This code can be used free of charge for academic research and private use, provided
% that you cite "Simon Burgess, Yubin Kuang, Kalle Åström, TOA sensor network self-calibration for
% receiver and transmitter spaces with difference in dimension, Signal Processing (2014),
% 10.1016/j.sigpro.2014.05.034"  
%
%For questions, mail simonb(äth)maths.lth.se .
clearvars, close all
addpath('..\..\..\sharedRoutines')
addpath('..\..\..\\bundleadjusters')
addpath('..\')

%%


%% Constants
d=2; %[ 5];% %different dimensions that is to be stested. NOT IMPLEMENTED FOR LOOP FOR THIS YET
nbrRuns=10; %numer of constellations (i.e. testruns) for each distance and numer receiver-transmittor setup
freqOutliersVec=[ 0.04 0.08 0.12];
freqMissingData=0.01;
%numer of receivers and transmittors each run
nbrR=10;
nbrS=15;
tol=0.005; %tolerance for the RANSAC method
opts.v=0.0001;
opts.dim=d+1;
nbrRansacRunsOuter=300;
nbrRansaRunsInner=50;
colors=['b' ;'g' ; 'r'; 'y'; 'k'];
pointStyle=['d'; 's'; '<'];
%FLags for different execution options
doBundling=1;
doGuessL1Comparison=1;
totTime=0;

maxIterOuter=10;
maxIterInner=5;
%%
%falseMat__falsePosOutlier_falseNehOutliers=[sum(sum(abs((~indInliers-indOutliers))))/nbrR/nbrS sum(sum(abs((~indInliers-indOutliers)==1)))/nbrR/nbrS sum(sum(abs((~indInliers-indOutliers)==-1)))/nbrR/nbrS]
relErr          =zeros(length(freqOutliersVec),nbrRuns);
relErrBundled   =zeros(length(freqOutliersVec),nbrRuns);
relErrGuess     =zeros(length(freqOutliersVec),nbrRuns);
%relErrBefore=zeros(length(freqOutliersVec),nbrRuns); %for debug
isFailures=zeros(length(freqOutliersVec),nbrRuns);
falseOutlierMat=zeros(length(freqOutliersVec),3,nbrRuns);
falseOutlierMatReally =zeros(length(freqOutliersVec),3,nbrRuns);
isErrorSender=zeros(length(freqOutliersVec),nbrRuns);

for indFreqOutl =1:length(freqOutliersVec)
    
    for ii=1:nbrRuns
        indFreqOutl
        ii
        ri=[[zeros(1,d);rand(nbrR-1,d)-0.5] zeros(nbrR,1)]';
        opts.centers=ri;
        
        [D, gt]=simulateSynteticTOA(nbrR,nbrS,[0.2*sqrt(d) sqrt(d)/2],opts);
        gt.coords(end,:)=gt.coords(end,:).* sign(gt.coords(end,:));
        
        %adding outliers
        indOutliers=rand(nbrR,nbrS)<freqOutliersVec(indFreqOutl);
        nbrOutliers=sum(sum(indOutliers));
        outliersMeas= (rand(nbrOutliers,1))*sqrt(d); %All meaus rents, outlers or not, should be >0
        D(indOutliers)=outliersMeas;
        indOutliersReally=(D-gt.distances)>tol; %Yeah, really!
        %Adding missing data
        indMissingData=rand(nbrR,nbrS)<freqMissingData;
        indMissingData(indOutliers)=0; %DOn't mark data as missing if it is an outlier
        D(indMissingData)=NaN;
        
        indInliersGt=(indOutliers+ indMissingData) ==0;  %Missing data is
        indInliersGtReally=(indOutliersReally+ indMissingData)==0;
        %clear indOutliers, indMissingData
        assert(sum(sum(D<0))==0)
        tic;
        [bestRec,bestSend,nbrValidCases,indInliers]=diffdimRansac(D,d,nbrRansacRunsOuter, nbrRansaRunsInner ,tol);
        thisRansacTime=toc;
        totTime=totTime+thisRansacTime;
        
        
        
        if nbrValidCases >0
            %[~, T,bestModelTransformed]=rigidTransform([bestRec, bestSend],[gt.centers gt.coords]); %fordebug
            resBefore=resErr(D,bestRec,bestSend,indInliers); %for debug. Seems to be working thoug. coincides with rec inside bundler
            %relErrBefore(indFreqOutl,ii) = norm(bestModelTransformed-[gt.centers gt.coords],'fro')/norm([gt.centers gt.coords],'fro'); %fordebug
            if doBundling %Bundle over the onmes classified as inliers
                %Obs! Not sure if bundler works as it should
                Dvec=D(indInliers(:));
                [indRec, indSend]=find(indInliers);
                [xopt,yopt]=bundleToaMics2DSounds3D(Dvec,indRec,indSend,bestRec(1:(end-1),:),bestSend,0);
                resAfter=resErr(D ,[xopt; zeros(1,nbrR)],yopt,indInliers); %For debug
                bestRecBund=[xopt; zeros(1,nbrR)]; %add the zeros again
                bestSendBund=yopt;
                [~, T,bestModelTransformedBundeled]=rigidTransform([bestRecBund, bestSendBund],[gt.centers gt.coords]);
                relErrBundled(indFreqOutl,ii) = norm(bestModelTransformedBundeled-[gt.centers gt.coords],'fro')/norm([gt.centers gt.coords],'fro');
           
                   if resBefore < resAfter   %debug
                       keyboard
                   end
            end
            
         
            [~, T,bestModelTransformed]=rigidTransform([bestRec, bestSend],[gt.centers gt.coords]);
            relErr(indFreqOutl,ii) = norm(bestModelTransformed-[gt.centers gt.coords],'fro')/norm([gt.centers gt.coords],'fro');
            
            falseOutlierMat(indFreqOutl,:,ii)=[sum(sum(abs((~indInliers-~indInliersGt))))/nbrR/nbrS    sum(sum(abs((~indInliers-~indInliersGt)==1)))/nbrR/nbrS    sum(sum(abs((~indInliers-~indInliersGt)==-1)))/nbrR/nbrS];
            falseOutlierMatReally(indFreqOutl,:,ii)=[sum(sum(abs((~indInliers-~indInliersGtReally))))/nbrR/nbrS    sum(sum(abs((~indInliers-~indInliersGtReally)==1)))/nbrR/nbrS    sum(sum(abs((~indInliers-~indInliersGtReally)==-1)))/nbrR/nbrS];
            assert(sum(indInliers(indMissingData))==0) %for debug 
            
            if doGuessL1Comparison
               %Ok, so we want to time i  
               
               totTimeGuess=0;
               nbrGuesses=0;
               relErrGuess(indFreqOutl,ii)=Inf;
               while totTimeGuess < thisRansacTime
                   tic
                   riGuess=[[zeros(1,d);rand(nbrR-1,d)-0.5] zeros(nbrR,1)]';
                   
                   [~, gtGuess]=simulateSynteticTOA(nbrR,nbrS,[0.2*sqrt(d) sqrt(d)/2],opts);
                   sjGuess=gtGuess.coords;
                   [I,J] = find(isfinite(D));
                   Dvec = D(isfinite(D));
                   [xoptGuess,yoptGuess,resVec]= bundleL1IteratedToaMics2DSounds3D(Dvec, I,J, riGuess(1:(end-1),:), sjGuess, 1, maxIterOuter, maxIterInner);
                   
                   
                    [~, ~,ModTransBundGuess]=rigidTransform([[xoptGuess; zeros(1,nbrR)] yoptGuess],[gt.centers gt.coords]);
                   nbrGuesses=nbrGuesses+1;
                   currRelErr=norm(ModTransBundGuess-[gt.centers gt.coords],'fro')/norm([gt.centers gt.coords],'fro');
                   if currRelErr< relErrGuess(indFreqOutl,ii)
                       relErrGuess(indFreqOutl,ii)=currRelErr;
                       
                   end
                   totTimeGuess=totTimeGuess + toc;

               end
                
            end
        else
            isFailures(indFreqOutl,ii)=1;
        end
    end
end


%%Histograms
%OBS! FIXA FALIURE CASES SÅ DI PRESNTERSA MEN EJ PLOTTAS!
maxbin=max(log10(relErr(isFailures(:)==0)));
minbin=min(log10(relErr(isFailures(:)==0)));
figure
hold on
for ii=1:size(freqOutliersVec,2)
       nbins=12
    linewidth = 2;
    fontsize = 14;
    relErrsToPlot=relErr(ii,~isFailures(ii,:));
    [aa,bb] = hist(log10(relErrsToPlot),linspace(minbin,maxbin,nbins));
    plot(bb,aa,strcat(colors(ii),'-',pointStyle(ii)),'linewidth',linewidth);
    set(gca,'fontsize',fontsize);
      set(gca,'XTick',[-16:0.5:0])
    xlabel('log_{10}(rel errors)','fontsize',19);
    ylabel('Frequency','fontsize',19); 
end
axis([-4 0,0,35])
legend([ num2str(100*freqOutliersVec(1)) '% outliers'],[ num2str(100*freqOutliersVec(2)) '% outliers'],[ num2str(100*freqOutliersVec(3)) '% outliers'])
hold off


if doBundling  %Make a hist of the bundled rsluts as well
    maxbinBund=max(log10(relErrBundled(isFailures(:)==0)));
    minbinBund=min(log10(relErrBundled(isFailures(:)==0)));
    figure
    hold on
    for ii=1:size(freqOutliersVec,2)
        nbins=12
        linewidth = 2;
        fontsize = 14;
        relErrsToPlot=relErrBundled(ii,~isFailures(ii,:));
        [aa,bb] = hist(log10(relErrsToPlot),linspace(minbinBund,maxbinBund,nbins));
        plot(bb,aa,strcat(colors(ii),'-',pointStyle(ii)),'linewidth',linewidth);
        set(gca,'fontsize',fontsize);
        set(gca,'XTick',[-16:0.5:0])
        xlabel('log_{10}(rel errors)','fontsize',19);
        ylabel('Frequency','fontsize',19);
    end
    axis([-4 0,0,35])
    legend([ num2str(100*freqOutliersVec(1)) '% outliers'],[ num2str(100*freqOutliersVec(2)) '% outliers'],[ num2str(100*freqOutliersVec(3)) '% outliers'])
    hold off
    
end

if doGuessL1Comparison  %Make a hist of the bundled rsluts as well
    maxbinBund=max(log10(relErrGuess(isFailures(:)==0)));
    minbinBund=min(log10(relErrGuess(isFailures(:)==0)));
    figure
    hold on
    for ii=1:size(freqOutliersVec,2)
        nbins=12
        linewidth = 2;
        fontsize = 14;
        relErrsToPlot=relErrGuess(ii,~isFailures(ii,:));
        [aa,bb] = hist(log10(relErrsToPlot),linspace(minbinBund,maxbinBund,nbins));
        plot(bb,aa,strcat(colors(ii),'-',pointStyle(ii)),'linewidth',linewidth);
        set(gca,'fontsize',fontsize);
        %set(gca,'XTick',[-16:0.5:0])
        xlabel('log_{10}(rel errors)','fontsize',19);
        ylabel('Frequency','fontsize',19);
    end
    %axis([-4 0,0,35])
    legend([ num2str(100*freqOutliersVec(1)) '% outliers'],[ num2str(100*freqOutliersVec(2)) '% outliers'],[ num2str(100*freqOutliersVec(3)) '% outliers'])
    hold off
    
end

   failurePercentage = sum(isFailures,2)./nbrRuns
falseOutlierAverage=mean(falseOutlierMat,3)
falseOutlierAverageReally=mean(falseOutlierMatReally,3)
nbrRelErrsBig=sum(sum(abs(isErrorSender-(relErr>10^-2))))

