%A script for producing comparative figures of the solutions between the reconstructed microphone and spekaer
%positions using the proposed algorithm and unsing computer vision.
%The data used is for the second real experiment in the paper "TOA sensor network self-calibration for
% receiver and transmitter spaces with difference in dimension".
%
%This code can be used free of charge for academic research and private use, provided
% that you cite "Simon Burgess, Yubin Kuang, Kalle Åström, TOA sensor network self-calibration for
% receiver and transmitter spaces with difference in dimension, Signal Processing (2014),
% 10.1016/j.sigpro.2014.05.034"  
%
%For questions, mail simonb(äth)maths.lth.se . 

close all,clearvars
 
 addpath('..\.')
 addpath('..\..\..\lowlevelTOAmatching\experimentTOAMeasurements\experimentResults_20140405_TOADIFFDIM\.')
 
load DMatrix_7mics2D_5sounds3D_20140415
load CVReconstructedCoords_First7Mics2D_Last4SOunds3D
 
isMirrorSol=1;
 d=2;
 [ xSol, ySol , x] = solveRecLowerDimTOAv2(D,d)
 setDownPlane=0.05;
 markerSize=10;
 LineWidth=1.1;
 fontsize = 15;
 
 [~,CV_aligned]=procrustes([xSol ySol].', CV_reconstructedPts.')
 CV_aligned=CV_aligned.';
 
 
if isMirrorSol
    CV_aligned=[CV_aligned(1,:); -CV_aligned(2,:); CV_aligned(3,:)];
    xSol=[ xSol(1,:); -xSol(2,:); xSol(3,:)];
    ySol=[ ySol(1,:); -ySol(2,:); ySol(3,:)];
    
end
 gt.centers=CV_aligned(:,1:7);
 gt.coords=CV_aligned(:,8:end);
%  fontsize = 19;
%     set(gca,'fontsize',fontsize);
%     xlabel('x','fontsize',fontsize);
%     ylabel('y','fontsize',fontsize);
%     zlabel('z','fontsize',fontsize);
% tempR=[ri ; ones(1,nbrR)];
% tempS=[sj ; ones(1,nbrS)];
% riT=T*tempR;
% sjT=T*tempS;
% riT=riT(1:end-1,:);
% sjT=sjT(1:end-1,:);
%

    


plot3(gt.centers(1,:),gt.centers(2,:),gt.centers(3,:),'b+','MarkerSize',markerSize,'LineWidth',LineWidth)
hold on
plot3(xSol(1,:),xSol(2,:),xSol(3,:),'r+','MarkerSize',markerSize)
plot3(gt.coords(1,:),gt.coords(2,:),gt.coords(3,:),'bx','MarkerSize',markerSize,'LineWidth',LineWidth)
plot3(ySol(1,:),ySol(2,:),ySol(3,:),'rx','markers',markerSize,'LineWidth',LineWidth)
axis equal
for ii=1:length(ySol(1,:))
    plot3([ySol(1,ii) ySol(1,ii)],[ySol(2,ii) ySol(2,ii)], [ySol(3,ii) 0],'r' )
end

[X,Y]=meshgrid([-0.5:0.1:0.6],[-0.1:0.1:1.2]);
Z=zeros(size(X))-setDownPlane;
mesh(X,Y,Z)


set(gca,'fontsize',fontsize);
xlabel('x (m)','fontsize',fontsize);
ylabel('y (m)','fontsize',fontsize);
zlabel('z (m)','fontsize',fontsize);

axis equal
hold off

set(gca,'XTick',[-0.5 0 0.5])
set(gca,'YTick',[0 0.5 1])
set(gca,'ZTick',[0 0.3 0.6])


RMSE= sqrt(sum(sum(([xSol ySol]-CV_aligned).^2))/size([xSol ySol],2))
% figure
% plot(xSol(1,:),xSol(2,:),'ro')
% hold on
% plot(ySol(1,:),ySol(2,:),'bo')
% axis equal