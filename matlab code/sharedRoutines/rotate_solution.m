function [Q1,Q2,res,transl1,transl2,t1,t2]=rotate_solution(pos1,pos2)
% Just estimate the coordinate transform to take c1 to pos 
% x-axis and c2 into pos xy plane

% STEP 1: translate to 0 for both solutions
transl1=pos1(:,1);
transl2=pos2(:,1);
pos1=pos1-repmat(transl1,1,size(pos1,2));
pos2=pos2-repmat(transl2,1,size(pos2,2));

% STEP 2: rotate so that second coordinate points in pure x direction
[Q1,R1]=qr(pos1(:,2:4));
[Q2,R2]=qr(pos2(:,2:4));

% STEP 3: Compensate for mirroring
Q1=(Q1*[sign(R1(1,1)) 0 0;0 sign(R1(2,2)) 0;0 0 sign(R1(3,3))])';
Q2=(Q2*[sign(R2(1,1)) 0 0;0 sign(R2(2,2)) 0;0 0 sign(R2(3,3))])';

if nargout>2
    t1=Q1*pos1;
    t2=Q2*pos2;
    res=mean(sqrt(sum((t2-t1).^2,1))); %WARNING! This is NOT rmse, but mrse (not used as much)
end

    %own sign function. Needed to not get sign(0)-->0 and then FUCK stuff
    %up. Goddamnit, why is this not fixed yet!
    function a=sign(b)
        if b>=0
            a=1;
        else
            a=-1;
        end
        
    end
end
