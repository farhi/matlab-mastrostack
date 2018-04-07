function points = find_control_points(im, N)
  % find_control_points: find N points separated with tol_trans*10 pixels
  
  points = [];
  if isempty(im), return; end
  
  points.x  = [];
  points.y  = [];
  points.m  = [];
  points.sx = [];
  points.sy = [];
  points.sharpness = [];
  
  % we need N points within the image. Area per point is prod(size(im))/N
  % exclusion distance is sqrt(prod(size(im))/N) (full width)
  dx = sqrt(prod(size(im))/N)/2;  % and we half it

  % find 'max' intensity locations. Every time we have a spot, we 'blank' it
  % around so that other ones are separated by 'tol_trans*10'
  for p=1:N
    [x1, y1, m1, im, sx, sy, iter, sharpness] = max_and_zero(im, dx/2, dx/2);
    if ~x1 || ~y1, continue; end % can not find a new star
    points.x(end+1) = x1;
    points.y(end+1) = y1;
    points.m(end+1) = m1;
    points.sx(end+1)= sy;
    points.sy(end+1)= sx;
    points.sharpness(end+1) = sharpness;
  end

  % sort points with increasing 'y' value
  [points.y,p] = sort(points.y);
  points.x = points.x(p);
  points.m = points.m(p);
  points.sx= points.sx(p);
  points.sy= points.sy(p);
  points.sharpness = points.sharpness(p);
  points.handle = [];
end % find_control_points
