% A script to see the histograms of the minimal solvers relative error.
%
%This code can be used free of charge for academic research and private use, provided
% that you cite "Simon Burgess, Yubin Kuang, Kalle Åström, TOA sensor network self-calibration for
% receiver and transmitter spaces with difference in dimension, Signal Processing (2014),
% 10.1016/j.sigpro.2014.05.034"
%
%For questions, mail simonb(äth)maths.lth.se .
clearvars

addpath('..\..\..\sharedRoutines\.')
addpath('..\')

%% Constants
dvec=[ 2 3 4]; %[ 5];% %Dimension of underlying subspace
nbrConst=100; %numer of constellations (i.e. testruns) for each distance and numer receiver-transmittor setup
stdErr=0;%0.001; %WHite gaussion noise added to the TDOA measurements



%%
errsDiffdim=zeros(length(dvec),nbrConst);
errsSameDim=zeros(length(dvec),nbrConst);
realSols=zeros(length(dvec),nbrConst);
colors=['b' ;'g' ; 'r'; 'y'];
pointStyle=['d'; 's'; '<'];
totTime=0;
hold on

for ii=1:length(dvec)
    
    nbrR= 1+dvec(ii) + dvec(ii)*(dvec(ii)+1)/2 ;
    nbrS=dvec(ii)+1 ;
    
    for jj=1:nbrConst
        
        
        ri=[[zeros(1, dvec(ii));rand(nbrR-1,dvec(ii))-0.5] zeros(nbrR,1)]';
        
        opts.v=stdErr; %0.001;
        opts.centers=ri;
        opts.dim=dvec(ii)+1;
        [D, gt]=simulateSynteticTOA(nbrR,nbrS,[sqrt(dvec(ii))*0.1 sqrt(dvec(ii))],opts);
        gt.coords(end,:)=gt.coords(end,:).* ((gt.coords(end,:) >=0) -0.5)*2;
        %% Run tests
        
        
        
        
        tic
        [ri,sj]=solveRecLowerDimTOAv2(D,dvec(ii));%toa_diffdim(D,dvec(ii));%
        %[riWorking,sjWorking]=solveRecLowerDimTOAv2(D,dvec(ii))
        t1=toc;
        totTime=totTime+t1;
                %keyboard
                %% solution and gt are only matched up till rotation
                if isreal([ri,sj]) 
                    realSols(ii,jj)=1;
                    [~, T]=rigidTransform([ri sj],[gt.centers gt.coords]);
                    %errsDiffdim(ii,jj)=sqrt(9)*errsDiffdim(ii,jj)/norm([gt.centers gt.coords]); %9* to make it the norm
                    
                    %TEST so it really is the relative norm
                    solAligned=T*[ri sj; ones(1,size([ri sj],2))];
                    solAligned=solAligned(1:end-1,:);
                    errsDiffdim(ii,jj)=norm(solAligned-[gt.centers gt.coords])/norm([gt.centers gt.coords]);
                end


        % calculate average errors and failure rate
        
    end
    

    
end

maxbin=max(log10(errsDiffdim(:)));
minbin=min(log10(errsDiffdim(:)));
figure; hold on
for ii=1:size(dvec,2)
       nbins=20
    linewidth = 2;
    fontsize = 14;
    [aa,bb] = hist(log10(errsDiffdim(ii,:)),linspace(minbin,maxbin,nbins));
    plot(bb,aa,strcat(colors(ii),'-',pointStyle(ii)),'linewidth',linewidth);
    set(gca,'fontsize',fontsize);
    set(gca,'XTick',[-16:1:-6])
    xlabel('log_{10}(rel errors)','fontsize',19);
    ylabel('Frequency','fontsize',19); 
end
 axis([-16 -9,0,160])
legend(strcat('p=',num2str(dvec(1))),strcat('p=',num2str(dvec(2))),strcat('p=',num2str(dvec(3))))
hold off