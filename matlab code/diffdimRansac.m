function [bestRec,bestSend,nbrValidCases,isInlier]=diffdimRansac(D,d,nbrRunsOuter, nbrRunsInner, tol,nbrOutliersPerColTol,nbrRunReestSend)
%DIFFDIMRANSAC runs ransac on the measurement matrix D, and returns the
%best found receiver and transmitter position coordinates for the case when receivers lie in a
%d-dimensional affine subspace and transmitters lie in a higher dimensional space which receivers
%are embedded in. 
%
%Transmitters, sj, recide in a d+1-dimensional subspace, and  receivers, ri, are assumed to embedded
%in a lower dimensional subspace of dimension d. D(i,j) is distance measurement between receiver i
%and transmitter j. 
%
%INPUT
%D          m x n       The m rows represent distance measurements from m receivers,
%                       and the n columns represent measurenemts from n
%                       transmitters. If no distance measurement is
%                       available between receiver i and transmiter j, set
%                       D(i,j)=NaN.
%d          scalar      Positive integer of dimensions of the affine
%                       subspace that receivers are located in.
%nbrRuns    scalar      Positive interger of how many randomizations the ransac loop should do.
%
%nbrOutliersPerColTol scalar    Postivie integer that sets how many
%                               outliers we can accept per column without trying to
%                               reestimate the sender for that colum in th end.
%                               (THis exists as each sender is only estiamted by
%                               triangulating form the minimal number of receivers.
%                               THis can be slightly unstable for high number of outliers, as for d=2, we only have six receivers, so six measurments to do ransac over each sender, whcih has three unknowns).
%                               (default=20% of m, rounded up)
%nbrRunReestSend        scalar  Positive integer, of how many ransac tries one
%                               should do to reestimate the sender with outliers
%                               above nbrOutliersPerColTol. 
%                               (default=nbrRunsInner)
%OUTPUT
%bestRec    d+1 x m     The best found constellation of coordinates for
%                       receiver positions. Will be just one 0 if no real
%                       solutions were found amongst the nbrRuns
%                       testruns. Last row should be zeros.
%sjBest     d+1 x n     The best found constellation of coordinates for
%                       transmitter positions. Will be just one 0 if no real
%                       solutions were found amongst the nbrRuns
%                       testruns.
%nbrValidCases scalar   Positive integer of how many of the nbrRuns
%                       runs gave a valid starting minimal case
%indInliers m x n       Boolean matrix that has 1 if corresponding element
%                       in Measurement matrix D is determined to be inler, or 0 for outlier or missing.

%% Initialization
[m,n]=size(D);

if nargin <6
    nbrOutliersPerColTol=ceil(m*0.2);
end

if nargin <7
   nbrRunReestSend=nbrRunsInner;
end

bestSend=zeros(d+1,n); %Will hold all the senders, one for each column of D, for the best found RANSAC fit.
%check so that we have the minimum number or receivers
bestNbrInilers=-Inf;
nbrValidCases=0;
bestRecMin=0;
bestIndRec=0;

assert(m>=1+d+d*(d+1)/2)
assert(n>=d+1)

%each loop is one randomization of a minimal case, and then see how good that fits to the unused columns corresponding to the used rows.
for ii=1:nbrRunsOuter
    indRec=randsample(m,1+d+d*(d+1)/2); %select the receivers indices for minimal case. (Naming receiver index --->
    indSend=randsample(n,d+1);  %Select transmitter indices. Naming: SendMin index ---> dCOl index
    
    
    if sum(sum(~isfinite(D(indRec,indSend))))~=0 %if thereare elements that are not finite, ex NaN
        continue
    end
    
    [recMin,sendMin]=solveRecLowerDimTOAv2(D(indRec,indSend),d);
    
    if ~isreal([recMin sendMin]) %Bad minimal case with no real solution
        continue
    end
    
    
    %OBS! RIght now, senders atree returned every time. It really just
    %needs to be returned for the best run.
    [send,nbrInliersSendExtend,allSendsIdent]=extendSenders(indRec,indSend,recMin,sendMin, D,d,tol);
    
    if allSendsIdent==0
        continue
    end
    
    nbrValidCases=nbrValidCases+1;
    
    if nbrInliersSendExtend>bestNbrInilers
        bestRecMin=recMin;
        bestIndRec=indRec;
        bestSend=send;
        bestIndSend=indSend; %For debug Purposes
        bestNbrInilers=nbrInliersSendExtend;
        %bestNbrInliersSendExtend=nbrInliersSendExtend;
    end
    
end

if nbrValidCases==0
    bestRec=0;
    bestSend=0;
    isInlier=0;
    warning('Could not find any valid cases of the nbrRunsOuter cases' )
    return
end

%Now ransac-triangulate up the receivers from the best fit
[bestRec,~,allReceiversIdent]=extendReceivers(bestIndRec,bestRecMin,bestSend,D,d,tol,nbrRunsInner);

if allReceiversIdent==0
    bestRec=0;
    bestSend=0;
    isInlier=0;
    warning('Could not find all the receivers. Possibly use more outerRuns, or too many missing data in D' )
    return
end

%saker kvar: *Skicka tillbaks inlier indices
isInlier=cast(zeros(m,n),'logical');
for ii=1:m
    for jj=1:n
        isInlier(ii,jj)= abs(D(ii,jj) -sqrt(sum((bestRec(:,ii)-bestSend(:,jj)).^2))) <tol;
    end 
end

nbrInliersPerCol=sum(isInlier);
isReestimateSends =  nbrInliersPerCol < (m-nbrOutliersPerColTol); 

isReestimated=zeros(1,n);
for jj=find(isReestimateSends) %reestimate senders for these columns
    [sj,nbrInliersInCol]=reestimateSend(bestRec,bestSend,D,d,tol,jj,nbrRunReestSend);
    
    if size(sj,1)~=1 && (nbrInliersInCol >nbrInliersPerCol(jj))
        bestSend(:,jj)=sj;
        isReestimated(jj)=1;
    end
    
end

%%AND update isInlier forthis column DOOO IT! OBS!
for jj=find(isReestimated)
    for ii=1:m
        
        isInlier(ii,jj)= abs(D(ii,jj) -sqrt(sum((bestRec(:,ii)-bestSend(:,jj)).^2))) <tol;
        
    end
end

end

function [sendReest,nbrInliers]=reestimateSend(recs,sends,D,d,tol,jj,nbrRuns)

%If no reestimation was found, zero is returned
[m,n]=size(D);

bestNbrFits=-Inf;
nbrComplexSj=0;
for kk=1:nbrRuns %Do a ransac loop to find better sender with more iniliers
    indRec = randsample(m,d+1);
    
    if sum(~isfinite(D(indRec,jj)))~=0
        continue
    end
    
    sj=triangulateSendMin(recs(:,indRec), indRec ,D,jj);
    
    if ~isreal(sj)
        nbrComplexSj=nbrComplexSj+1; %for debug
        continue
    end
    
        diff=D(:,jj)-sqrt(sum((recs-repmat(sj,1,m)).^2,1)).'; 
        nbrFits=sum(abs(diff)<tol);
        
        
        if  nbrFits>bestNbrFits
            bestSj=sj;
            bestNbrFits=nbrFits;
        end
        
end

if exist('bestSj','var')
   sendReest=bestSj; 
   nbrInliers=bestNbrFits;
else
    sendReest=0;
    nbrInliers=0;
end


end



function [rec,nbrInliers,allReceiversIdent]=extendReceivers(indRecKnown,recKnown, sendKnown,D,d,tol,nbrRansacTriangs)
%Assuming that we have all senders,  but only a subset of the receivers,
%specified which exalty in indRecKnown, ransac-triangulate up the unknown
%senders, and count inliers
%OUTPUT
%rec            d+1 x m     All receivers
%nbrInliers  scalar      Number of inliers for the rows of D that is to
%                           be found new receiers for.
[m,n]=size(D);
rec=zeros(d+1,m); %Could be done to just get the missing receivers, but this makes the main code look nicer
rec(:,indRecKnown)=recKnown;
nbrInliers=0;
allReceiversIdent=1;
for ii=setdiff(1:m,indRecKnown)
    %For each unkown receiver, do a small ransac
    bestNbrInliers=-Inf; %Number of How many of the 1+d+d*(d+1)/2 mesurements that fits best sender Hypothesis
    clear bestRi
    for kk=1:nbrRansacTriangs
        indSend=randsample(n,d+1); %holds which d+1 of the m known senders to be used for tiangulation. (In reality, only d is keeded, but then we get two solution. With d+1, we can use a linear solver.
        if sum(~isfinite(D(ii,indSend)))~=0
            continue
        end
        
        ri=triangulateRec(sendKnown(:,indSend),indSend,D,ii);
        
        diffVec=D(ii,:)-sqrt(sum((repmat(ri,1,n)-sendKnown).^2)); %difference between the measiruemts and the model
        nbrFits=sum(abs(diffVec)<tol);
        
        if nbrFits > bestNbrInliers
            bestNbrInliers=nbrFits;
            bestRi=ri;
            
        end
    end
    
    
    if exist('bestRi','var')
        rec(:,ii) = bestRi;
        nbrInliers = nbrInliers + bestNbrInliers;
    else
        %If we did not find ANY receiver to triangulate up (slim, but
        %theoretically possible)
        allReceiversIdent=0;
        return
    end
    
    
    
end

end


function ri=triangulateRec(send, indSends ,D,indRec)
%THe minimal case requires d senders. THis case, using a linear solver,
%requires >=d+1 senders. SO not really minimal, but close to. THis gies
%only one solution for th receiver triangulated , instead of two.
%SOlve for r,1 sj,2 ... sj,d
d=size(send,1)-1;

A=-2*(send(1:(end-1),2:end) - repmat(send(1:(end-1),1),1,d)).'; 
y=(   D(indRec,indSends(2:end)).^2 - D(indRec,indSends(1)).^2   ).';
y=y + sum(send(:,1).^2) - sum(send(:,2:end).^2,1).';
%This linear equations has d unknowns and >=d equations

riD=A\y;





%riD=A\y; %TECKENFEL! Med riktiga lösningar, har y precis omvänt tecken!

ri=[riD; 0];
end




function [send,nbrInliers,allSendersIdent]=extendSenders(indRecKnown,indSendKnown,recKnown,sendKnown,D,d,tol)
%for each new sender, we have 1+d+d*(d+1)/2 equations/measurements and d+1 unknowns, as
%some of the measurements might be outliers, we do a small, exhaustive, ransac loop on
%this too
%OBS! DOES NOT WORK WELL WITH MISSING DATA RIGHT NOW!
%INPUT
% indRecTot (1+d+d*(d+1)/2 ) x 1    The indices of D for the receivers
%                                   that we have coordinatec in rec for
%indSend    1+d x 1                 The indices of senders that we allready know
%                                   coordinates for. (the ones not to be
%                                   tiangulated from the known receivers)
%rec        1+d x (1+d+d*(d+1)/2 )  The known receivers, where rows corresponds to coordinates and columns to receivers 
%D          m x n                   The distance measuremnt matrix
%d          scalar                  Dimension of the supspace receivers
%                                   inhabit
%tol        scalar                  Distance tolerance.
%                                   ||rec_i-send_j|| -D(i,j) < tol, it is
%                                   considerend an inlier.
%
%OUTPUT

%send       1+d x (n)               The extended senders, plus the original sender found by the
%                                   exahustive ransac loop. Rows are
%                                   coordinates and columns are different
%                                   senders, corresponding to columns in D.
%nbrInliers scalar                  Looking at D(:,setdiff(1:size(D,2),indSendKnown)) (ie the part in D that is used to
%                                   extend senders), count how many of the
%                                   measurements are within tol of
%                                   ||rec_i-send_j||. (THe minimum number
%                                   will be
%                                   (d+1)*length(setdiff(1:size(D,2),indSendKnown)),
%                                   as for each column, d+1 measurements
%                                   will always fit as this is the minimal
%                                   triangulation case).
%notAllSendersIdent                 Boolean that tells ut that at least one
%                                   sender cannot be real, whatvever d+1
%                                   receivers you use to try to solve the
%                                   minimal trinagulation case.
%
%OBS! ONE scary thing... When searching through the six (1+d+d*(d+1)/2) values that mathces the
%known receivers, for each column: There might not be ANY three
%combinations of measurements of these six that gives a real solution to
%the sender.

send=zeros(1+d,size(D,2));
send(:,indSendKnown)=sendKnown;
nbrInliers=0;
[mKnown]=size(recKnown,2);
allSendersIdent=1;
for jj=setdiff(1:size(D,2),indSendKnown) %Loop over the unused senders and find one for each unused column of D
    bestNbrFits=-Inf; %Number of How many of the 1+d+d*(d+1)/2 mesurements that fits best sender Hypothesis
    clear bestSj;
    indRecMat=nchoosek(1:mKnown,d+1); %each row corresponds to whic receivers to be used
    nbrComplexSj=0; %for debug purposes.
    
    for ii=1:size(indRecMat,1)
        
        
        %the d receivers used   %the index of those receivers
        
        if sum(sum(~isfinite(D(indRecKnown(indRecMat(ii,:)),jj))))~=0 %if thereare elements that are not finite, ex NaN
            continue
        end
        sj=triangulateSendMin(recKnown(:,indRecMat(ii,:)), indRecKnown(indRecMat(ii,:)) ,D,jj);
        
        if ~isreal(sj)
            nbrComplexSj=nbrComplexSj+1;
            continue
        end
        %OBS! THERES THEE LINES MIGHT LOOK FASTER THAN THE FOURTH BUT SETDIFF IS VERY SLOW!
%         indExtraReceiversInD =setdiff(indRecKnown,indRecKnown(indRecMat(ii,:)),'stable'); %Can be sped up by removing ad doing over all?
%         indExtraReceivers=setdiff(1:(1+d+d*(d+1)/2),indRecMat(ii,:),'stable');
%         diff=D(indExtraReceiversInD,jj) - sqrt(sum((recKnown(:,indExtraReceivers) -repmat(sj,1,(1+d+d*(d+1)/2)-(d+1)) ).^2)).'; %Difference of ||ri-sj||-d_ij for said column and receivers. Or rather, the part of the column that has not has its values perfectly fit to sj. Example: IF there are six senders on the partial comun of D we are looking at, then three of the values if d=2 are chosen so that sj fits perfectly to them. THis row of code checks the diff of the other three.
        diff=D(indRecKnown,jj)-sqrt(sum((recKnown-repmat(sj,1,(1+d+d*(d+1)/2))).^2,1)).'; 
        nbrFits=sum(abs(diff)<tol);
        if ~(nbrFits  >=d+1) %%TEST
            keyboard
        end
        
        if nbrFits == mKnown %Maximum number of fits. We can stop ii-loop here here. THis check spped things up.
            bestNbrFits=nbrFits;
            bestSj=sj;
            break
        end
        
        
        if  nbrFits>bestNbrFits
            bestSj=sj;
            bestNbrFits=nbrFits;
        end
    end
    
        %Might be that we have not found a single real sender.
                                %We eant to have found at least 2 extra
                                %(i..e. d+3)
                                %inliers to belive it. Else it might be
                                %that we have an outlier sender model that fits
                                %ine more data point. (might stillbe,b ut
                                %less probable)
    if exist('bestSj','var') & bestNbrFits >= d+1 %exist('bestSj','var') &  WARNING! THIS DOES NOT WORK WITH d=1
        %DEBUG PURPOSES
%         if sum(sum((send-repmat(bestSj,1,size(send,2)))==0))~=0
%             keyboard
%         end
        
        send(:,jj)=bestSj;
        
        nbrInliers=nbrInliers + nbrFits;
    else
        allSendersIdent=0;
        return
    end
    
    
end
end



function sj=triangulateSendMin(rec, indRecTot ,D,indSend)

%SOlve for sj,1 sj,2 ... sj,d
d=size(rec,1)-1;
A=2*(-rec(1:(end-1),2:end)+repmat(rec(1:(end-1),1),1,d)).'; % A looks ok
y=D(indRecTot(2:end),indSend).^2 - D(indRecTot(1),indSend).^2;
y=y+sum(rec(1:d,1).^2)-sum(rec(1:d,2:end).^2,1).';
sjD=A\y; 

sjLast=sqrt(D(indRecTot(1),indSend).^2-sum((sjD-rec(1:(end-1),1)).^2)  ); %THis looks ok as the first residual is 0
%sjLast might become imaginary. If we only look for real soultions we can
%punsih this in several ways. One is to just set it to zero.
sj=[sjD; sjLast];
end





