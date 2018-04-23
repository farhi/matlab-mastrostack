function filename = write_stacked(im, filename, info)
  % we create a cell with PNG options for imwrite
    % clean the EXIF
    if nargin < 3, info = struct(); end
    f = fieldnames(info);
    for index=1:numel(f)
      if isfield(info, f{index})
        this = info.(f{index});
        flag = ~isempty(this) && ~isstruct(this) && (ischar(this) || (isnumeric(this) && isscalar(this)));
        if ~flag
          info = rmfield(info, f{index});
        end
      end
    end
    % fix DateTime
    if isfield(info,'DateTime')
      try
        d = datenum(info.DateTime);
      catch
        index = find(info.DateTime == ':');
        if numel(index) == 4
          info.DateTime(index(1:2)) = '-';
        end
      end
    else info.DateTime = datestr(now);
    end
    desc = '';
    try
      desc = evalc('disp(info)');
    end
    if isempty(desc)
      try
        desc = disp(info);
      end
    end
    if ~isfield(info,'Software'),    info.Software    = mfilename; end
    if ~isfield(info,'XResolution'), info.XResolution = 350; end
    if ~isfield(info,'YResolution'), info.YResolution = 350; end
    if ~isfield(info,'Model'),       info.Model       = version; end
    if ~isfield(info,'Make'),        info.Make       = 'Matlab'; end
    info.Filename     = filename;
    info.FileModDate  = datestr(now);
    info.Software     = mfilename;
    
    info.FileSize     = numel(im);
    info.Width        = size(im,2);
    info.Height       = size(im,1);
    info.BitDepth     = size(im,3)*8;
    if size(im,3) == 1
      exif.ColorType    = 'grayscale';
    else
      exif.ColorType    = 'truecolor';
    end
    args={ ...
      'XResolution', info.XResolution , ...
      'YResolution', info.YResolution , ...
      'Title', [ mfilename ' stacked images ']  , ...
      'Author', 'E. Farhi' , ...
      'Description', desc , ...
      'CreationTime', info.DateTime , ...
      'Software', info.Software , ...
      'Source', [ info.Model ' ' info.Make ] , ...
      'Comment', desc };
    % writing stacked image
    disp([ mfilename ': writing ' filename ])
    warning('off','MATLAB:writepng:changedTextChunk')
    try
      imwrite(im, filename, args{:});
    catch
      imwrite(im, filename);
    end
