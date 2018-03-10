function sharpness = image_sharpness(im)
  % image_sharpness: estimate sharpness from sum of gradients
  
  im = double(im);
  % a blurred image has smooth variations. We sum up diff
  im1 = abs(diff(im,[], 1))/numel(im);
  im2 = abs(diff(im,[], 2))/numel(im);
  sharpness = sum(im1(:))+sum(im2(:));
