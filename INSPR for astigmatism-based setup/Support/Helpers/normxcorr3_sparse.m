function C = normxcorr3_sparse(T, A, shape)
% C = normxcorr3(TEMPLATE, IMAGE, SHAPE)
%
%       TEMPLATE - type double, ndims==3, size <= size of image
%       IMAGE    - type double, ndims==3
%       SHAPE    - one of: 'valid', 'same', 'full'. same as conv2 shape parameter
%                  'full' by default
%
%       C        - values in [-1,1]. size depends on SHAPE
%
% the syntax of this function is identical to Matlab's
% normxcorr2, except that it's been extended to 3D matrices,
% and, the SHAPE parameter has been introduced as a convenience
%
% the SHAPE parameter has the same effect as it does for the CONVN function.
% see the documentation for CONVN for a more detailed explanation
%
% caveat emptor: this function does not perform the argument checking that
% normxcorr2 does. for example, it doesn't ensure that std(T(:))~=0
%
% daniel eaton, 2005, danieljameseaton@gmail.com



if nargin<3
	shape = 'full';
end

if ndims(A)~=3 || ndims(T)~=3
	error('A and T must be 3 dimensional matrices');
end

szT = size(T);
szA = size(A);

if any(szT>szA)
	error('template must be smaller than image');
end

pSzT = prod(szT);

% make the running-sum/integral-images of A and A^2, which are
% used to speed up the computation of the NCC denominator
intImgA = integralImage(A,szT);
intImgA2 = integralImage(A.*A,szT);

szOut = size(intImgA);

% compute the numerator of the NCC
% emulate 3D correlation by rotating templates dimensions
% in 3D frequency-domain correlation is MUCH faster than the spatial-domain
% variety
rotT = flip(flip(flip(T,1),2),3); % this is rot90 in 3d
fftRotT = fftn(rotT,szOut);
clear rotT
fftA = fftn(A,szOut);
clear A;
corrTA = real(ifftn(fftA.*fftRotT));
clear fftA;
clear fftRotT;
num = (corrTA - intImgA*sum(T(:))/pSzT ) / (pSzT-1);
clear corrTA
% compute the denominator of the NCC
denomA = sqrt( ( intImgA2 - (intImgA.^2)/pSzT ) / (pSzT-1) );
clear intImgA
clear intImgA2
denomT = std(T(:));
clear T
denom = denomT*denomA;

% compute the NCC
s = warning('off', 'MATLAB:divideByZero');
C = num ./ denom;
clear denom
clear num
s = warning('on', 'MATLAB:divideByZero');

% replace the NaN (if any) with 0's
zeroInd = find(denomA==0);
C(zeroInd) = 0;
clear denomA
clear zeroInd
switch( lower(shape) )
	case 'full'
	case 'same'
		szTp = fix((szT-1)/2);
		C = C( szTp(1)+1:szTp(1)+szA(1), szTp(2)+1:szTp(2)+szA(2), szTp(3)+1:szTp(3)+szA(3) );
	case 'valid'
		C = C(szT(1):end-szT(1)+1,szT(2):end-szT(2)+1,szT(3):end-szT(3)+1);
	otherwise
		error(sprintf('unknown SHAPE %s, assuming FULL by default', shape));
end

function integralImageA = integralImage(A,szT)
% this is adapted from Matlab's normxcorr2

szA = size(A);
szs=szA+2*szT-1;
B=[];
for ii=1:1:szs(3)
B {ii}= sparse(szs(1),szs(2));
end

for ii=(szT(3)+1):(szT(3)+szA(3))
% B( szT(1)+1:szT(1)+szA(1), szT(2)+1:szT(2)+szA(2), szT(3)+1:szT(3)+szA(3) ) = A;
B{ii}(szT(1)+1:szT(1)+szA(1), szT(2)+1:szT(2)+szA(2)) = A(:,:,ii-szT(3));
end

c=[];
for ii=1:numel(B)
    tmp=cumsum(B{ii},1);
    c{ii}=tmp(1+szT(1):end,:)-tmp(1:end-szT(1),:);
end

c2=[];
for ii=1:1:numel(c)
    tmp=cumsum(c{ii},2);
    c2{ii}=tmp(:,1+szT(2):end)-tmp(:,1:end-szT(2));
end

c2sz1=size(c2{1});
c2sz3=numel(c2);
s=[];
s=cell2mat(c2);
s=reshape(full(s),[c2sz1(1) c2sz1(2) c2sz3]);
s = cumsum(s,3);
integralImageA = s(:,:,1+szT(3):end)-s(:,:,1:end-szT(3));

% s = cumsum(B,1);
% c = s(1+szT(1):end,:,:)-s(1:end-szT(1),:,:);
% s = cumsum(c,2);
% c = s(:,1+szT(2):end,:)-s(:,1:end-szT(2),:);
% s = cumsum(c,3);
% integralImageA = s(:,:,1+szT(3):end)-s(:,:,1:end-szT(3));




