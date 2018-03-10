function [s1,s2, f1, f2] = peakwidth(im, x,y, dx, dy)
% peakwidth: determine a peak width along X and Y
%
% input:
%   im:     image (rgb)
%   x,y:    peak location
%   dx,dy:  pixel size to extract around peak (optional, default=+/-5)

  if nargin < 4, dx=4; end
  if nargin < 5, dy=4; end

  % fist extract the portion of the image to analyze
  X=(x-dx):(x+dx); X=X(X>0 & X<size(im,1));
  Y=(y-dy):(y+dy); Y=Y(Y>0 & Y<size(im,2));
  
  % then project the image along dimension 1
  [s1, f1] = width1(Y, im(x,Y));
  [s2, f2] = width1(X, im(X,y)); 
end % peakwidth
