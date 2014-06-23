function [D, centers, nvects] = simulateCellDists(param,m,k,centers_in);
%simulateCellDists silumates cell phone relative motions and directions to
%the base stationss
%   
%INPUT
%   m - number of motions
%   k - number of base stations
%OUTPUT:
%   D - measurements of the relative distances
%   centers - m by 3 matrix: cell phone positions of the m motions
%   nvects  - 3 by k matrix: directions to the k base stations, with unit
%   L2-norm

if nargin<1,
    param = 1;
end

if nargin<2,
    m = 4;
end

if nargin<3,
    k = 6;
end
if nargin==4
    centers=centers_in;
else
    centers = randn(3,m);
end
is_on_conic=1;
while is_on_conic
    nvects = randn(3,k);
    nvects(:,3) = cross(nvects(:,1),nvects(:,2));
    nvects = psphere(nvects);
    D = centers'*nvects;
    [s]=svd(D);
    if all(abs(s(1:3))>0.1)
        is_on_conic=0;
    end
end

