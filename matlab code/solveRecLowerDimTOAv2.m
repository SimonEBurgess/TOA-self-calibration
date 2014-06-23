function [ xSol, ySol , isPosDef] = solveRecLowerDimTOAv2(D,d)
%SOLVERECLOWERDIMv2 Solves the problem of determining receiver positions in a lower dimension space (dimension d) and
%transmitters in higher dimension space (dimension d+1) from TOA measurents. Solution is unique up to
%rotation/mirroring of coordinate system and translation, and up to
%mirroring of tranmittors in plane of receivers, giving 2^nbrTransmitters solutions up to translation and rotation/mirroring.
%The algorithm needs at least 1 + d + (d+1)d/2 receivers and d+1 transmitters.
%THis is a minimal algorithm, i.e. it in general solves the sets of equations exactly with finite
%number of soultions, in this case one. (up to rot, transl and mirroring in d-dim).
%
% DISCLIAMER: This code can be used free of charge for academic research and private use, provided
% that you cite "Simon Burgess, Yubin Kuang, Kalle Åström, TOA sensor network self-calibration for
% receiver and transmitter spaces with difference in dimension, Signal Processing (2014),
% 10.1016/j.sigpro.2014.05.034"  
%
%  VERSION 2: The Difference from the first version of the solver is how the S-matrix
%first factorization is made. Here, an svd is made and then a d-dimensional
%invertible matrix B0 with cond=1 is calcualted so that the right hand side
%factorization has ones in the last row. See article or ask SiBu for a more
%thorough explanation. 
%  2012-09-04: Version 2 and version 1 should be equivalent if all assumtions hold, but
%version 2 har better numerical performance (but also 20% slower),
%especially in comparing upperbounds. For 1000 random problems in 3 dimensions, mean
%error was  [meanErrV1 meanErrV2]=[4.4180e-011  5.6413e-013]. The
%difference is mostly explained on that the upper bound is much lower for
%version 2.
%
%
%INPUT:
%   D   m x n       Matrix of TOA measurements. Receivers index
%                   conicides with row intex. m >= (1+d+(d+1)d/2)
%                   and n>=d+1.
%   d   scalar      Dimension of affine space receivers recide in
%   
%
%OUTPUT
%   xSol    (d+1)x m    Coordinates of the receivers in the
%                       higher dimensional space. Each column is a
%                       receiver position.
%   ySol    d+1 x n     Coordinates of transittors in the higher 
%                       dimensional space. Each colum is a position.
%   isPosDef scalar     Boolean that says if real solutions exist or
%                       not.
% Receivers and transmitters are symmetrical for TOA, so feel free to
% exchange them. In this function, receivers are assumed to be in a the lower dimensional subspace,
% e.g. receivers in a plane and transmitters in space. 



%% Constants
[nbrRec, nbrTrans]=size(D);

 if nbrRec < (1+d+(d+1)*d/2)
   disp('Too few receivers for given dimension d of subspace');
   %return
elseif nbrTrans < d+1
    disp('Too few transmitters for given dimension d of subspace. ');
    %return
end

%% Making the squared distance difference amtrix S
S=D.^2;
S=S-repmat(S(1,:),nbrRec,1);
S=S(2:end,:);

%% Making the First factorization: S=X0*Y0
[U,s,V]=svd(S,'econ');
X0=U(:,1:d+1);
Y0=s(1:d+1,1:d+1)*V(:,1:d+1)'; %v(:,1:4)'
%% Second factorization: Making left hand side have ones furhtest down. S=X1*Y1
b0=ones(1,size(Y0,2))/Y0; %making the part so that Y1 will have last row of ones (or closest to in least square sense)
B1=[eye(d) zeros(d,1) ; b0];

B1=gs(B1); %one loop of classical gram schmidt to make it better condition number. Makes it 20% slower, but better overall numerical precision.

X1=X0/B1;
Y1=B1*Y0;

%% SOlving for the unknown transformation A 

C=zeros(nbrRec-1, d+(d+1)*d/2);
C(:,1:d)=(X1(:,1:(d)).^2)/4;  %Coeffcients in front of diagonal

%rest of coefficients
nbrUsedPos=d; %number of used colums
for ii=1: (d-1)
   %Each for loop takes care of coefficients in front of off diagonal elements of (upper triangular part) of row ii of A^t*A
   C(:,(nbrUsedPos+1):(nbrUsedPos +(d-ii))) = 0.5*diag(X1(:,ii))*X1(:,ii+1:d);
   nbrUsedPos=nbrUsedPos+ d-ii;
end

%coefficients in front of b
C(:,end-d+1:end)=-X1(:,1:d);

y=X1(:,end);

params=C\y; %ok so far (at least when d=2)

%ugly hack to reshape matrix from params
AtA=diag(params(1:d));
nbrUsedPos=d;  
for ii=1:d-1
    AtA(ii,ii+1:end)=params(nbrUsedPos+1:nbrUsedPos+d-ii);
    AtA(ii+1:end,ii)=AtA(ii,ii+1:end);  % put it in lower triangular part also!
    nbrUsedPos=nbrUsedPos+d-ii;
end

%factorizing AtA=(A^T)*A into A
[U,s2]=eig(AtA);
Atop=U*sqrt(s2);

try 
    AA = chol(AtA);
    isPosDef = 1;
catch
    isPosDef = 0;
end

A=[Atop params(end-d+1:end) ; zeros(1,d) 1];

%% Constructing solutions
leftSol=X1*A; %left part of factorization of S
rightSol=A\Y1;   %Right part of factorization of S.  So S=leftSol*rightSol

xSol=[zeros(1,d) ;-leftSol(:,1:d)/2]; 
xSol=[xSol zeros(nbrRec,1)].';
ySol=rightSol(1:d,:);%Now we only have the coordinates in the plane of transmitters for y. We need to collect them all *nanananana POKEMON!*

normSquaredYproj=sum(ySol.^2);

zCoords=sqrt(D(1,:).^2-normSquaredYproj);
ySol=[ySol; zCoords];

    function A=gs(A)
       %Does one loop of ~Gram schmidt on the columns of A.
       %
       % Keeps the right colum the same, and then does gram-schmidt on the
       % other columns, right from left, to ortogonalize the matrix.
       % FInally scale each column to have the same length as the rightmost
       % col. (everyting to make B have a lower condition numer, and tuhs
       % more precision when inverted)
       
       
       [m,n]=size(A);
       Q=zeros(m,n);
       R=zeros(n,n);
       len=norm(A(end,:));
       
       for j=n-1:-1:1 %j=1:n                       %loop n times
           %v=A(j,:);
           for i=n:-1:j+1 %1:j-1                 % This loop runs total (n^2-n)/2
               A(j,:)=A(j,:) - A(i,:)*A(j,:)' *A(i,:)/(A(i,:)*A(i,:)');
           end
           A(j,:)=A(j,:)*len/norm(A(j,:));             % O(m) Trash flops (few)
           %Q(:,j)=v/R(j,j);            % m Trash flops (few)
       end
       
    end


end

