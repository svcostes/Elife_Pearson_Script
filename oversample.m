function out_img =  oversample(img,coef,mode)
% Will increase the size of an image, interpolating from initial image img
% using nearest neighbor interpolation.
% Inverse functin of subsample.j
% img = oversample(img,[w h d],METHOD); Resample image increasing current width to
% w, current height to h and current depth to d
% METHOD: interpolation method. Default is 'nearest'
% vailable methods are:
%  
%       'nearest' - nearest neighbor interpolation
%       'linear'  - linear interpolation
%       'spline'  - spline interpolation
%       'cubic'   - cubic interpolation as long as the data is uniformly
%                   spaced, otherwise the same as 'spline'
%
% Sylvain Costes, September 2009, Lawrence Berkeley National Lab
% see also interp3

% Input matrix size
[w,h,d] = size(img);

if ~exist('mode','var')
    mode = 'nearest';
end

% Output matrix size
W = coef(1);H = coef(2);
try
    D = coef(3);
catch
    D = 1;
end

% Current input matrix grid
X = repmat((1:w),[h 1 d]);
Y = repmat((1:h)',[1 w d]);
Z = repmat(reshape(1:d,1,1,d),[h w 1]);

% Matrix to be interporlated
XI = repmat((1:W)/W*w,[H 1 D]);
YI = repmat((1:H)'/H*h,[1 W D]);
ZI = repmat(reshape(1:D,1,1,D)/D*d,[H W 1]);

% Interpolation
if d>1
    out_img = dip_image(interp3(X,Y,Z,double(img),XI,YI,ZI,mode));
else
    out_img = dip_image(interp2(X,Y,double(img),XI,YI,mode));
end