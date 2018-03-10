function im = imdouble(im)
  % imdouble: cast an image to a double with values in [0-1]

  if isempty(im), return; end
  switch class(im)
  case 'uint8'
    im = double(im)/2^8;
  case 'uint16'
    im = double(im)/2^16;
  case 'double'
    % do nothing
  otherwise
    im = double(im)/max(im(:));
  end
end % imdouble
