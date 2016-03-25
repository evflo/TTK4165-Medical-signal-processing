function logP=imagelog(P,gain,dyn);

%IMAGELOG   Image a matrix of ultrasound power in log scale
%
%   function logP=imagelog(P,gain,dyn);
%
%   Inputs: P    - Power image
%           gain - gain [dB] (=image brightness)
%           dyn  - dynamic range [dB] (=image contrast)
%
%   Output: logP - Output image, log Power.
%
%   If no output image is requested, the image is displayed
%   in 64 graylevels in the current axes. If the input image is
%   complex valued, only the real part is displayed.

% (c) GE Vingmed Ultrasound 1997-99. All rights reserved
% Written by: HT, 04.07.96
% Last update: 14.08.97, JK

if nargin<3, error('Too few input arguments');end;

if ~isreal(P),
    disp('Warning: Imaginary part of image discarded');
end;

if prod(size(gain))~=1|~isreal(gain)|isstr(gain),
    error('Invalid gain. Must be a real number');
end;

if prod(size(dyn))~=1|~isreal(dyn)|dyn<=0|isstr(dyn),
    error('Invalid dynamic range. Must be real and positive');
end;

minIndex=1;
maxIndex=64;
SmallNonZeroNumber=10^(-maxIndex); %To avoid log of zero

logP=maxIndex*(gain+10*log10(max(SmallNonZeroNumber,real(P))))/dyn;
logP=min(maxIndex, max(minIndex,logP));

if ~nargout,
    image(logP);
    colormap(gray(maxIndex));
    logP=[];
end;