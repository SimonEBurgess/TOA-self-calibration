function [D, groundtruth] = simulateSynteticTOA(m,k,d,opts) %d,v,centers, type, magOff,dim)
%Simulates sytnhetic TOA or UTOA measurements, and ground truth
%
%INPUT
%   m       - scalar  number of motions +1, or number of receivers
%   k       - scalar  number of base stations, or number of transmittors
%   d       - scalar  distance to basestations from ORIGO, give 2 vaules for interval
%              If you want to simulate true far field mesurements, set to NaN.
%   opts    - a stuct with the following fields                  
%           v       - scalar>=0 standard deviation of distance measurement
%                     (default = 0).
%           centers - 3 x number receivers-matrix  Centers of the receivers (defualt randn(dim,m)). This overrides
%                     m as number o receivers if given.
%           coords  - dim x number transmitters-matrix. Coordinates for transmitters to be sued. (Only works
%                     for far field as of now. Default is randomized
%                     uniformly over sphere). THis overrides k as number of
%                     transmitters
%           type    - String. If set to 'TDOA', the measurements will be TDOA instead of
%                       TOA. If set to 'UTOA', the measurements will be set to be
%                       UTOA. Default is TOA measuements.
%           magOff  - 2*1 - magnitude of the offsets. (Default is randn*10 for each offset)
%           dim     - scalar  Dimension of the underlying affine space (default 3).
%                   (Does not work with far field)
    
%OUTPUT:
%   D - measurements of the  distances
%   groundtruth.centers - coordinates of recievers
%   groundtruth.angles - angles of antennas (ONLY GIVEN BACK for dim=3, far
%                        field measurements as of right now)
%   groundtruth.distances - mxk matrix of distance measures from antennas
%       locations WITHOUT NOISE, each row a set of measurements
%   groundtruth.coords - coordinates of transmitters, or for true far field
%                        simulations they are the direction from
%                        transmittors to receivers
%   groundtruth.recOffset - offsets of the receivere
%   groundtruth.transOffset -offset for the transmittors
%   groundtruth.noise       -Additive noise.
%   groundtruth.direction
if nargin <2,
    error('Missing arguments\n');
end



% generate centers
if nargin <4
    opts.dummy=1; %dummy variable
end

if ~isfield(opts,'dim')
    opts.dim=3;
end

if (~isfield(opts,'centers'))
    groundtruth.centers = randn(opts.dim,m);
else
    groundtruth.centers=opts.centers;
    m=size(opts.centers,2);
end

if (~isfield(opts,'type'))
   opts.type='TOA';
end

if (~isfield(opts,'v'))
   opts.v=0;
end

if ~isfield(opts,'magOff')
   opts.magOff=[10 ; 10]; %should be in the range of the TOA differences               %and for receivers: Seconds i.e. 5s*340m/s=1000
end



if opts.dim~=floor(opts.dim) & opts.dim <1
    disp('wrong dimension. Must be integer > 0')
    return
end

% generate angles. Warning: angles cannot be picked uniformely if we wat a
% uniform distribution overthe hypersphere. We can do it in 3 dimensons,
% but no general formula that i know of exists for d dimensins (SiBu).
% Hence, groundtruth.angles and  spherical coordinates will be osolete
%groundtruth.angles = rand(opts.dim-1,k).*repmat([2*pi ;ones(opts.dim-2,1)*pi],1,k);


if ~isnan(d) %if we have real distances and not simulated far field
    if length(d)==2 % <-- should ve generate distances?
        dd=rand(1,k)*(d(2)-d(1))+d(1);
    else % otherwise just use d
        dd=ones(1,k)*d;
    end
    % calculate x,y,z for the points
  
    %Cartesian coodinated for the affine space of dimension dim
    groundtruth.coords=zeros(opts.dim,k);
%     if opts.dim==1
%         groundtruth.coords(1,:)=dd;
%     else
%         groundtruth.coords(1,:)=dd.*cos(groundtruth.angles(1,:));
%         groundtruth.coords(2,:)=dd.*sin(groundtruth.angles(1,:)); %Looks ok
%         
%         if opts.dim>=3
%             for ii=3:opts.dim
%                 groundtruth.coords(1:ii-1,:) = groundtruth.coords(1:ii-1,:).*repmat( sin(groundtruth.angles(ii-1,:)) ,ii-1,1);
%                 groundtruth.coords(ii,:)=dd.*cos(groundtruth.angles(ii-1,:));
%             end
%         end
%     end
    groundtruth.coords = randn(opts.dim,k);
    coordNorms=sqrt(sum(groundtruth.coords.^2,1));
    groundtruth.coords=groundtruth.coords ./repmat(coordNorms,opts.dim,1).*repmat(dd,opts.dim,1);
        

    % calculate real distances to the centres
    groundtruth.distances=zeros(m,k);
    
    
    for ii=1:k %<-- should be doable without loop, but not important
        %Warning. Huge loss of precision if distances are a lot bigger than the
        %TOA differences.
        groundtruth.distances(:,ii)=sqrt(sum((groundtruth.centers-repmat(groundtruth.coords(:,ii),1,m)).^2,1))';
    end
    
    groundtruth.direction = zeros(opts.dim,k);
    
    mrec = mean(groundtruth.centers,2);
    
    groundtruth.direction = groundtruth.coords - repmat(mrec,1,k);
    groundtruth.direction = groundtruth.direction./(ones(opts.dim,1)* sqrt(sum( groundtruth.direction.^2)));
    
else %if we want to simulate far field
    %the directons from transmittors to receivers. Only Relveant for TDOA
    %and UTOA
    %only works for dim=2
    
    if (~isfield(opts,'coords')) %if we havent got our coords for transmiter directions, randomize them
        if opts.dim == 3;
            groundtruth.angles = zeros(opts.dim-1,k);
            groundtruth.angles(1,:) = rand(1,k).*2*pi;
            groundtruth.angles(2,:) = acos(2*rand(1,k)-1);
            
            groundtruth.coords(1,:)=sin(groundtruth.angles(2,:)).*cos(groundtruth.angles(1,:));
            groundtruth.coords(2,:)=sin(groundtruth.angles(2,:)).*sin(groundtruth.angles(1,:));
            groundtruth.coords(3,:)=cos(groundtruth.angles(2,:));
        elseif opts.dim == 2;
            groundtruth.angles = zeros(opts.dim-1,k);
            groundtruth.angles(1,:) = rand(1,k).*2*pi;
            
            groundtruth.coords(1,:)=sin(groundtruth.angles(1,:));
            groundtruth.coords(2,:)=cos(groundtruth.angles(1,:));
        else
            groundtruth.coords = randn(opts.dim,k);
            groundtruth.coords = bsxfun(@rdivide,groundtruth.coords, sqrt(sum(groundtruth.coords.^2)));
        end
    else  %if we got our tranmitter coordiantes from options
        groundtruth.coords=opts.coords;
        k=size(opts.coords,2);
        
    end
    
    groundtruth.distances=groundtruth.centers' * groundtruth.coords;
    
    
end



groundtruth.noise=opts.v*randn(m,k);
D=groundtruth.distances + groundtruth.noise;


if strcmp(opts.type,'UTOA')
    %add differences randomly
    groundtruth.transOffset=rand(1,size(D,2))*opts.magOff(1); 
    groundtruth.recOffset=rand(size(D,1),1)*opts.magOff(2);
    D = D + repmat(groundtruth.transOffset, size(D,1),1) + repmat(groundtruth.recOffset,1,size(D,2));
    
elseif strcmp(opts.type,'TDOA')
        groundtruth.transOffset=randn(1,size(D,2))*opts.magOff(1);
        D = D + repmat(groundtruth.transOffset, size(D,1),1);
        
    end
end