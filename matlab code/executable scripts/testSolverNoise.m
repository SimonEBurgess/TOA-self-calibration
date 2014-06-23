% A script for comparing soluting for different number or receivers and transmitters, the reliatve
% error under different additive noise and the rate of complex solutions.
%
%This code can be used free of charge for academic research and private use, provided
% that you cite "Simon Burgess, Yubin Kuang, Kalle Åström, TOA sensor network self-calibration for
% receiver and transmitter spaces with difference in dimension, Signal Processing (2014),
% 10.1016/j.sigpro.2014.05.034"  
%
%For questions, mail simonb(äth)maths.lth.se .
clearvars, %close all
addpath('..\..\..\sharedRoutines')
addpath('.\..')
%% Constants
dvec=2; %[ 5];% %different dimensions that is to be stested. NOT IMPLEMENTED FOR LOOP FOR THIS YET
nbrConst=100; %numer of constellations (i.e. testruns) for each distance and numer receiver-transmittor setup
stdErr=10.^[-9  -7 -5 -4 ];% %WHite gaussion noise added to the TDOA measurements
nbrRS=[  6 3 ;7 4 ; 20 10]; %number or ris ans sj respectively
 %numer of receivers and transmittors each run


%%
errsDiffdim=zeros(length(stdErr),nbrConst,size(nbrRS,1));
errsSameDim=zeros(length(stdErr),nbrConst,size(nbrRS,1));
fail=zeros(length(stdErr),nbrConst,size(nbrRS,1));
colors=['b' ;'g' ; 'r'];
pStyle=['*-' ; '<-' ; '>-'];
totTime=0;




for ii=1:length(stdErr)
    
    for kk=1:size(nbrRS,1)
        
        nbrR=nbrRS(kk,1);
        nbrS=nbrRS(kk,2);
        
        for jj=1:nbrConst
            
            
            ri=[[zeros(1,dvec);rand(nbrR-1,dvec)-0.5] zeros(nbrR,1)]';
            
            opts.v=stdErr(ii); %0.001;
            opts.centers=ri;
            opts.dim=dvec+1;
            [D, gt]=simulateSynteticTOA(nbrR,nbrS,[0.1*sqrt(dvec) sqrt(dvec)],opts); %
            gt.coords(end,:)=gt.coords(end,:).* ((gt.coords(end,:) >=0) -0.5)*2;
            %% Run tests
            
            
            
            
            tic
            [ri,sj]=toa_diffdim(D,dvec); %solveRecLowerDimTOAv2(D,dvec);%toa_diffdim(D,dvec);% solveRecLowerDimTOAv2(D,dvec);%toa_diffdim(D,dvec);%solveRecLowerDimTOAv2(D,dvec);
            
            if( ~(isreal(ri) & isreal(sj) ))
                if(~isreal(ri))
                    
                    %keyboard
                end
                fail(ii,jj,kk)=1;
            end
            t1=toc;
            totTime=totTime+t1;
            %keyboard
            %% solution and gt are only matched up till rotation
            if ~fail(ii,jj,kk)
                [~, T]=rigidTransform([ri sj],[gt.centers gt.coords]);
                %errsDiffdim(ii,jj)=errsDiffdim(ii,jj)/norm([gt.centers gt.coords]);
                solAligned=T*[ri sj; ones(1,size([ri sj],2))];
                solAligned=solAligned(1:end-1,:);
                errsDiffdim(ii,jj,kk)=norm(solAligned-[gt.centers gt.coords])/norm([gt.centers gt.coords]);
            end
            
            
            % calculate average errors and failure rate
            
        end
        
    end
    
    
end

%to make RMSE
meanRelErr=mean(errsDiffdim,2);
meanRelErr=reshape(meanRelErr,length(stdErr),size(nbrRS,1))
flre=zeros(length(stdErr),size(nbrRS,1));
for ii=1:length(stdErr)
    for kk=1:size(nbrRS,1)
        flre(ii,kk)=sum(fail(ii,:,kk))/length(fail(ii,:,kk));
    end
end
fontsize = 14;
set(gca,'fontsize',fontsize);


for ii=1:size(meanRelErr,2)
    loglog(stdErr,meanRelErr(:,ii)',strcat(colors(ii),pStyle(ii,:)),'LineWidth',2);
    hold on;
end
h1=gca;
set(h1,'YAxisLocation','right')
set(h1,'XLim',[stdErr(1)*0.5 stdErr(end)*2])


xlabel('standard deviation of error','fontsize',19);
ylabel('LINES - mean rel. error ','fontsize',19);
%title('Farfield Approximation performance');
legend('6 r - 3 s','7 r - 4 s','20 r - 10 s','Location','NorthWest')

h2 = axes('Position',get(h1,'Position'));
bar(log(stdErr),flre(),0.4);
set(h2,'YAxisLocation','left','Color','none','XTickLabel',[]);
set(h2,'XLim',log(get(h1,'XLim')),'Layer','top');
set(h2,'YLim',[0 0.5],'Layer','top');
set(gca,'fontsize',fontsize);
ylabel('BARS - Complex solution rate','fontsize',19); 



hold off



