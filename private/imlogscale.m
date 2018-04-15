function im = imlogscale(im)
% enhance the intensity contrast using a log-scale
  im=imdouble(im);
  im=log(1e-1+im*100);
  im=im-min(im(:));
  im=im2uint(im/max(im(:)));
end
