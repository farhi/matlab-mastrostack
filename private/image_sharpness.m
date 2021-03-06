function sharpness = image_sharpness(im)
  % image_sharpness: estimate sharpness from sum of gradients
  
  if ndims(im) == 3, im = rgb2gray(im); end
  im  = double(im);
  % a blurred image has smooth variations. We sum up diff
  im1 = diff(im,[], 1);
  im  = diff(im,[], 2);
  imax= min(numel(im1),numel(im));
  sharpness = sqrt(sum(im1(1:imax).^2+im(1:imax).^2))/imax;

end
