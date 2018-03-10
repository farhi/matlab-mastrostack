function [im, img, exif] = imread_single(source, images, flag)
  % imread_single: read a single image.
  %
  %   source: the image filename, or matrix, or image index
  %   images: the loaded images (structure array)
  %   flag:   when 0, the image is not read when not new
  
  if nargin < 3, flag=1; end

  im = []; img = []; exif = [];
  % new image
  if ischar(source) && exist(source, 'file')
    % check if image exists
    try
      exif     = imfinfo(source);
      im       = imread(source);
    catch ME
      disp([ mfilename ': Skipping ' source ]);
      return
    end
  elseif isnumeric(source) && isscalar(source)
    % source is an image index
    if numel(images) >= source
      img    = images(source);
      source = img.source;
      
      % make sure we can get the image (matrix), when img was retrieved but not image
      if flag && ischar(img.source)
        im = imread(img.source);
      end
    end
  elseif isnumeric(source)
    exif     = [];
    im       = source;
  else
    disp([ mfilename ': argument should be char or matrix. Is ' class(source) ])
    return
  end

  
  
  if isempty(im) return; end

  % create fake exif if not set
  if isempty(exif)
    exif.Filename     = sprintf('Image_%d', numel(images)+1);
    exif.FileModDate  = datestr(now);
    exif.FileSize     = numel(im);
    exif.Format       = 'jpg';
    exif.FormatVersion = '';
    exif.Width        = size(im,2);
    exif.Height       = size(im,1);
    exif.BitDepth     = size(im,3)*8;
    exif.Software     = mfilename;
    if size(im,3) == 1
      exif.ColorType    = 'grayscale';
    else
      exif.ColorType    = 'truecolor';
    end
  end
  
  % create the img structure holding information and EXIF
  img.index       = numel(images)+1;
  [~,img.id]      = fileparts(exif.Filename);
  img.source      = source;
  img.image       = im;
  img.image_size  = size(im);
  img.image_sum   = sum(double(im(:)));
  img.exif        = exif;
  img.thumbnail   = im;
  img.sharpness   = image_sharpness(im);
  
  if numel(img.thumbnail) > 1e6
    ratio         = ceil(sqrt(numel(img.thumbnail)/1e6));
    img.thumbnail = img.thumbnail(1:ratio:end, 1:ratio:end,:);
  end
    
  % these will be set later
  img.type        = [];
  img.points      = [];
  img.rotation    = [];
  img.translation = [];
  img.angle       = [];
  img.scaling     = [];
  img.reference   = [];
  img.ra          = [];
  img.dec         = [];
  img.focal_length= [];
  img.sensor_width= [];
  img.sensor_height=[];
  img.exposureTime =[];
  
end % imread_single
