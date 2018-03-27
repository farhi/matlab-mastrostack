function [x1,y1, m1, im1, sx, sy, iter, sharpness] = max_and_zero(im1, dx, dy, iter)
  % max_and_zero: search for a maximum intensity point and zero image around it
  %
  % input:
  %   im1: image (gray m*n)
  %   dx,dy: area to zero after search
  % output:
  %   x1,y1: location in image
  %   m1:    intensity
  
  if nargin <2, dx=0.3; end
  if nargin <3, dy=0.3; end
  if nargin <4, iter=1; end
  if isscalar(dx) && dx<1
    dx = round(size(im1,1)*dx);
  end
  if isscalar(dy) && dy<1
    dy = round(size(im1,2)*dy);
  end
  sharpness = 0;
  
  % search for max intensity
  [m1,x1]=max(im1(:)); m1 = double(m1);
  [x1,y1]=ind2sub(size(im1),x1);
  [sx,sy, f1, f2] = peakwidth(im1, x1,y1);
  % blank image around max
  dx1 = round((x1-dx):(x1+dx));
  dy1 = round((y1-dy):(y1+dy));
  dx1=dx1(dx1>=1 & dx1 <=size(im1,1));
  dy1=dy1(dy1>=1 & dy1 <=size(im1,2));
  sharpness = image_sharpness(im1(dx1,dy1));
  im1(dx1,dy1) = 0;
  
  if iter>5, return; end
  
  % recursive call if that guess is not acceptable
  % remove dead pixels (too sharp peaks) and image edges.
  if sx*sy < 4 || x1<5 || x1>size(im1,1)-4 || y1<5 || y1>size(im1,2)-4
    [x1,y1, m1, im1, sx, sy, iter, sharpness] = max_and_zero(im1, dx, dy, iter+1);
    return
  end
  
end % max_and_zero
