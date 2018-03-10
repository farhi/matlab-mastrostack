function points = find_control_points(im, N, tol_trans)
  % find_control_points: find N points separated with tol_trans*10 pixels
  points.x  = zeros(1, N);
  points.y  = zeros(1, N);
  points.m  = zeros(1, N);
  points.sx = zeros(1, N);
  points.sy = zeros(1, N);

  % find 'max' intensity locations. Every time we have a spot, we 'blank' it
  % around so that other ones are separated by 'tol_trans*10'
  for p=1:N
    [x1, y1, m1, im, sx, sy] = max_and_zero(im, tol_trans, tol_trans);
    if ~x1 || ~y1, continue; end % can not find a new star
    points.x(p) = x1;
    points.y(p) = y1;
    points.m(p) = m1;
    points.sx(p)= sx;
    points.sy(p)= sy;
  end

  % sort points with increasing 'y' value
  [points.y,p] = sort(points.y);
  points.x = points.x(p);
  points.m = points.m(p);
  points.sx= points.sx(p);
  points.sy= points.sy(p);
end % find_control_points
