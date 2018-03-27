function im = imlogscale(im)
% enhance the intensity contrast using a log-scale
  im=imdouble(im);
  im=log(1+im);
  im=im-min(im(:));
  im=im2uint(im/max(im(:)));
end
