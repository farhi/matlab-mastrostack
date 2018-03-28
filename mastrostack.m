classdef mastrostack < handle

  properties
    toleranceTranslation  = 0.01; % in percent
    toleranceRotation     = 3;    % in deg
    reference             = 0;    % by default the first 'light' image
    figure                = [];
    images = [];
    %      id
    %      source (or matrix if directly given)
    %      image_size
    %      image_sum
    %      exif
    %      type (light, dark, flat)
    %      points
    %      rotation
    %      translation
    %      scaling
    %      thumbnail
    %      ra
    %      dec
    %      focal_length
    %      sensor_width
    %      sensor_height
    %      exposureTime
    %      sharpness
    
    % the actual image matrix are not stored (except when given as such)
    % as we expect to handle 100's of images.
    dark    = 0;  % stored as a double in [0-1]
    flat    = 0;  % stored as a double in [0-1]
    light   = 0;  % stored as a double in [0-1]
    lightN  = 0;  % stored as a uint16
    currentImage = [];  % index of the current image on figure
    dndcontrol = [];
  end % properties
  
  methods
  
    % I/O stuff ================================================================
    function self = mastrostack(source)

      % load catalogs: objects, stars
      disp([ mfilename ': Welcome !' ]);
      
      if nargin % import images
        imread(self, source, 0);
      end
      
      plot(self);
      
    end % mastrostack
    
    function [im, img] = imread(self, source, flag)
      % imread: read an image, and store its information
      %
      %  [im, img] = imread(self, filename or matrix or index)
      %    returns the image matrix and information
      %
      %  The first 'light' loaded image defines the default reference image.
      
      % flag: when 0, do not force to read matrix images
      if nargin < 3, flag=1; end
      
      % first check if image(s) are already loaded
      [found, src] = exist(self, source);
      if ~iscell(src) && ~isstruct(src), src = { src }; end
      img   = []; im = {};
      if isempty(src), return; end
      
      for index=1:numel(found)
      
        % read new image, or retrieve it
        if isstruct(src)
          this_src = src(index);
        else
          this_src = src{index};
        end
        [this_im, this_img] = imread_single(this_src, self.images, flag);
        if ischar(this_src)
          if strfind(lower(this_src), 'dark'), this_img.type = 'dark';
          elseif strfind(lower(this_src), 'flat'), this_img.type = 'flat';
          end
        end
        
        if isempty(this_img) continue; end  % invalid image
        
        % store if this is a new image
        if isnan(found(index))
          this_img.index = numel(self.images)+1;
          if isempty(self.images)
            self.images    = this_img;
          else 
            self.images(this_img.index) = this_img; 
          end
        else % update
          this_img.index = found(index);
          self.images(this_img.index) = this_img; 
        end

        % add to return arguments
        if isempty(img), img        = this_img;
        else             img(end+1) = this_img; end
        if flag, im{end+1} = this_im; end
      end % for
      
      if iscell(im) && numel(im) == 1
        im = im{1};
      end
      
    end % imread
    
    function load(self, source)
      imread(self, source, 0);
    end % load
    
    function [found, src] = exist(self, source)
      % exist: check if an image is already loaded
      %
      %   found = exist(self, source)
      %     returns the index of any matching image, or nan (not found)

      % handle array of images
      found = []; src = {};
      if ischar(source)
        source = resolvefiles(source);
      end
      if iscellstr(source) && numel(source) == 1
        source = source{1}; 
      end
      if (iscell(source) && numel(source) > 1) || (isnumeric(source) && ~isscalar(source) ...
        && numel(source) == max(size(source))) || (isstruct(source) && numel(source) > 1)
        for index=1:numel(source)
          if iscell(source), this=source{index};
          else this=source(index); end
          [found(end+1), src{end+1}] = exist(self, this);
        end
        return
      end

      found = nan; src = source;
      if isnumeric(source) && isscalar(source) && source <= numel(self.images)
        found = source;
      else
        for index=1:numel(self.images)
          this = self.images(index);  % a struct
          if isempty(this), continue; end
          if ischar(source) && ~isempty(this.source) && ischar(this.source)
            % found: same name
            if strcmp(this.source, source) || strcmp(this.id, source)
              found = index;
              break; 
            end
          elseif isnumeric(source)
            % found: equal sum and size
            if all(size(this.image_size) == size(source)) ...
              && this.image_sum && sum(source(:)) == this.image_sum
              found = index;
              break; 
            end
          elseif isstruct(source)
            if any(strcmp(this.source, { source.source })) ...
            || any(strcmp(this.id, { source.id }))
              found = index;
              break; 
            end
          end
        end % for
      end
      
    end % exist
    
    % compute stuff ============================================================
    
    function im = getdark(self)
      % getdark: compute the mean Dark frame, using labelled images
      
      if any(~isfinite(self.dark)), self.dark=0; end
      if isempty(self.dark) || isscalar(self.dark)
        disp([ mfilename ': computing master Dark' ])
        sumDarks = 0;
        nbDarks  = 0;
        for index=1:numel(self.images)
          img = self.images(index);
          if ~strcmp(img.type, 'dark'), continue; end
          % read image (re-read if needed), convert to double [0-1]
          sumDarks  = sumDarks+ imdouble(imread(self, index));
          nbDarks   = nbDarks + 1;
        end
        if nbDarks>0
          self.dark = sumDarks/nbDarks;
          disp([ mfilename ':   used ' num2str(nbDarks) ' images.' ])
        else self.dark=0;
        end
      end
      im = self.dark;
      
    end % getdark
    
    function im = getflat(self)
      % getflat: compute the mean Flat frame, using labelled images
      
      if any(~isfinite(self.flat)), self.flat=0; end
      if isempty(self.flat) || isscalar(self.flat)
        disp([ mfilename ': computing master Flat' ])
        sumFlats = 0;
        nbFlats  = 0;
        for index=1:numel(self.images)
          img = self.images(index);
          if ~strcmp(img.type, 'flat'), continue; end
          % read image (re-read if needed), convert to double [0-1]
          im = imdouble(rgb2gray(imread(self, index)));
          im = im/max(im(:));
          sumFlats  = sumFlats+ im;
          nbFlats   = nbFlats + 1;
        end
        if nbFlats>0
          self.flat = sumFlats/nbFlats;
          index     = isfinite(self.flat) & self.flat > 0;
          self.flat = self.flat/max(max(self.flat(index)));
          index     = ~isfinite(self.flat) | self.flat <= 0;
          self.flat(index) = 1;
          disp([ mfilename ':   used ' num2str(nbFlats) ' images.' ])
        else self.flat=0; end
      end
      im = self.flat;
      
    end % getflat

    function im = correct(self, im, cl)
      % correct: correct an image for vigneting and background
      if nargin < 2, return; end
      if nargin < 3
        cl = class(im); % we shall then cast back
      end
      
      if isempty(im), return; end
      im = imdouble(im);
      
      % correct for flat field to compensate vigneting (divide)
      if ~isempty(self.flat) && ~isscalar(self.flat)
        for l=1:size(im,3)
          im(:,:,l) = im(:,:,l)./self.flat;
        end
      end
      
      % correct for read-out/sensor noise (subtract)
      if ~isempty(self.dark) && ~isscalar(self.dark)
        im = im - self.dark; % dark is a double
      end

      im = im2uint(im, cl);
      
    end % correct
    
    function plot(self, img, im)
    
      if nargin < 2
        % we search for the first 'light' image
        img = [];
        for index=1:numel(self.images)
          if strcmp(self.images(index).type, 'light')
            img = self.images(index);
            break
          end
        end
        if isempty(img), img = 1; end
      end
      
      fig = findall(0, 'Tag', 'mastrostack');
      if isempty(fig) % build the figure
        fig = figure('Name','Ma(e)stroStack','Tag','mastrostack', ...
          'MenuBar','none', 'Toolbar','figure',...
          'CloseRequestFcn',      {@MenuCallback, self }, ...
          'WindowButtonDownFcn',  {@ButtonDownCallback, self}, ...
          'KeyPressFcn',          {@KeyPressCallback, self }, ...
          'WindowScrollWheelFcn', {@ScrollWheelCallback,self});
        % add menu items and mouse/keyboard callbacks
        m = uimenu(fig, 'Label', 'File');
        uimenu(m, 'Label', 'Open',        ...
          'Callback', {@MenuCallback, self },'Accelerator','o');
        uimenu(m, 'Label', 'Open light images',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Open dark images',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Open flat images',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Save',        ...
          'Callback', 'filemenufcn(gcbf,''FileSave'')','Accelerator','s');
        uimenu(m, 'Label', 'Save As...',        ...
          'Callback', 'filemenufcn(gcbf,''FileSaveAs'')');
        uimenu(m, 'Label', 'Print',        ...
          'Callback', 'printdlg(gcbf)');
        uimenu(m, 'Label', 'Close',        ...
          'Callback', {@MenuCallback, self }, ...
          'Accelerator','w', 'Separator','on');
          
        m = uimenu(fig, 'Label', 'Image');
        uimenu(m, 'Label', 'Mark as Light image...',       ...
          'Callback', {@MenuCallback, self },'Accelerator','l');
        uimenu(m, 'Label', 'Mark as Dark image...',        ...
          'Callback', {@MenuCallback, self },'Accelerator','d');
        uimenu(m, 'Label', 'Mark as Flat image...',        ...
          'Callback', {@MenuCallback, self },'Accelerator','f');
        uimenu(m, 'Label', 'Mark as "skip"...',        ...
          'Callback', {@MenuCallback, self },'Accelerator','i');
        uimenu(m, 'Label', 'Mark as Reference image',        ...
          'Callback', {@MenuCallback, self },'Accelerator','r');
        uimenu(m, 'Label', 'Show in log scale',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Next image', 'Separator','on',       ...
          'Callback', {@MenuCallback, self },'Accelerator','n');
        uimenu(m, 'Label', 'Previous image',        ...
          'Callback', {@MenuCallback, self },'Accelerator','p');
        uimenu(m, 'Label', 'Goto image...',    ...
          'Callback', {@MenuCallback, self });
          uimenu(m, 'Label', 'Clear image(s)...', 'Separator','on',    ...
          'Callback', {@MenuCallback, self });
        
        m = uimenu(fig, 'Label', 'Compute');
        uimenu(m, 'Label', 'Set tolerances',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Automatic control points',        ...
          'Callback', {@MenuCallback, self },'Accelerator','a');
        uimenu(m, 'Label', 'Clear all control points',        ...
          'Callback', {@MenuCallback, self },'Accelerator','c');
        uimenu(m, 'Label', 'Show sharpness',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Compute master Dark',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Compute master Flat',        ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Align', 'Separator','on',       ...
          'Callback', {@MenuCallback, self });
        uimenu(m, 'Label', 'Stack',       ...
          'Callback', {@MenuCallback, self });
        % install mouse event handler on axis
        set(gca, 'ButtonDownFcn', {@ButtonDownCallback, self }); % will get its CurrentPoint
        set(fig, 'KeyPressFcn',   {@MenuCallback, self });       % will get its CurrentCharacter
        set(fig, 'WindowScrollWheelFcn',   {@ScrollWheelCallback, self });
        set(fig, 'UserData',self);
        
        t = {'\bf Welcome to {\color{magenta}Ma(e)stroStack} !', ...
          ' ', ...
          'This application allows to stack amateur astro-photography images', ...
          ' ', ...
          '1- import any {\color{blue}Dark} images (taken with the cap on the scope)', ...
          '2- import the {\color{blue}Flat} images (taken viewing a uniform surface)', ...
          '3- import the {\color{blue}Light} images (those showing stars and deep sky objects)', ...
          '4- compute master Dark and Flat', ...
          '5- select the {\color{red}Reference} image (on which all images will be added)', ...
          '6- stack ! (the alignment will then be done at the same time)', ...
          ' ', ...
          'You can look at images to mark those to ignore/skip (e.g. blurred)', ...
          'To use them back set their type to e.g. Light', ...
          ' ', ...
          '{\color{red}Use keys and menu items:}', ...
          '  L: mark current image as Light',...
          '  D: mark current image as Dark',...
          '  F: mark current image as Flat',...
          '  I: mark current image as Skip/Ignore',...
          '  R: mark current image as Reference',...
          '  N or ->: see next image', ...
          '  P or <-: see previous image', ...
          '  A: select automatically control points for alignment', ...
          '  click: add a control point in the image for alignment', ...
          ' ' };
        t = text(0.1,0.5, t, 'Units','normalized');
        
      else
        set(0, 'CurrentFigure',fig); % activate but no raise
      end
      self.figure = fig;
      
      if nargin >= 3 && ischar(im)
        option = im;
        im     = [];
      else option = '';
      end

      % get the image
      if nargin < 3 || isempty(im) || ~isnumeric(im)
        [im, img] = imread(self, img, 1);
      else
        [~, img] = imread(self, img, 0);
      end
      if isempty(img), return; end
      self.currentImage = img.index;
      t = [ self.images(self.currentImage).id ...
        ' ' self.images(self.currentImage).type ...
        ' sharp=' num2str(self.images(self.currentImage).sharpness) ];
      if strcmp(self.images(self.currentImage).type, 'light') && ...
        ((~isempty(self.flat) && ~isscalar(self.flat)) || ...
         (~isempty(self.dark) && ~isscalar(self.dark)))
        im = correct(self, im);
        t = [ t ' (/flat, -dark)' ];
      end
      t = [ t ' [' num2str(self.currentImage) '/' num2str(numel(self.images)) ']' ];
      
      if strcmp(option, 'log')
        im = imlogscale(im);
      end
      h = image(im);
      title(t);

      % display control points
      hold on
      self.images(self.currentImage).points.handle = scatter(...
        self.images(self.currentImage).points.y, ...
        self.images(self.currentImage).points.x, 100, 'g', 'x');
      hold off
      
    end
    
    function im = stack(self)
      % stack: stack all images on the reference
      
      % compute dark and flat (if not done yet)
      getdark(self);
      getflat(self);
      self.light = 0.0;
      self.lightN= uint8(0);
      ExposureTime = 0;
      
      % we first check for the reference and its control points
      for index=1:numel(self.images)
        this_img = self.images(index);
    
        if isempty(this_img.type) || strcmp('light', this_img.type)
          if ~self.reference || self.reference == index
            self.reference = index;
            this_img       = self.images(self.reference);
            self.images(self.reference).type = 'reference';
            disp([ mfilename ':   using reference as ' this_img.id ' (' num2str(index) ')' ])
            if isempty(this_img.points) || ~isstruct(this_img.points) || ...
              ~isfield(this_img.points, 'x') || numel(this_img.points.x) < 2
              % compute control points
              disp([ mfilename ':   computing control points for ' this_img.id ' (' num2str(index) ')' ])
              im = correct(self, imread(self, this_img)); % read and correct image
              this_img = cpselect(self, this_img, im);
            end
            if isfield(this_img.exif, 'ExposureTime')
              ExposureTime = ExposureTime + this_img.exif.ExposureTime;
            end
            break
          end
        end
      end
      if self.reference < 1 || self.reference > numel(self.images), return; end
      
      delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
      wb  = waitbar(0, [ mfilename ': Stacking images (close to abort)...' ]); 
      set(wb, 'Tag', [ mfilename '_waitbar' ]);
      t0 = clock; stackedimages = 0;
      
      for index=1:numel(self.images)
        this_img = self.images(index);
        
        if isempty(this_img.type) || strcmp('light', this_img.type)
          % check for the reference
          if ~self.reference || self.reference == index
            continue; 
          end
          
          % stack light images
          im = correct(self, imread(self, this_img)); % read and correct image
          if isempty(im), continue; end
          
          % get control points
          points = this_img.points;
          if isempty(points) || numel(points.x) < 2
            disp([ mfilename ':   computing control points for ' this_img.id ' (' num2str(index) ')' ])
            this_img = cpselect(self, this_img, im);
          end
          
          % display
          plot(self, this_img, im);
          
          % compute the affine transformation wrt reference (using control points)
          [ret_t, ret_R] = imdiff(self, self.images(self.reference), this_img);
          if isempty (ret_t), continue; end
          
          % we read the image and transform it
          [im,M]  = imaffine(imdouble(im), ret_R, ret_t);
          self.lightN= self.lightN+M;
          clear M
          self.light = self.light + im; % add rotated/translated image on the reference
          clear im
          
          % sum-up exposure
          if isfield(this_img.exif, 'ExposureTime')
            ExposureTime = ExposureTime + this_img.exif.ExposureTime;
          end
          
          % compute ETA
          stackedimages = stackedimages+1;
          dt_from0     = etime(clock, t0);
          dt_per_image = dt_from0/stackedimages;
          % remaining images: numel(s)-index
          eta    = dt_per_image*(numel(self.images)-index+1);
          ending = addtodate(now, ceil(eta), 'second');
          ending = [ 'Ending ' datestr(ending) ];
          eta    = sprintf('ETA %i [s]. %s', round(eta), ending);
          disp([ mfilename ': ' ...
              num2str(index) '/' num2str(numel(self.images)) ': ' this_img.id '. ' eta]);
          
          % update waitbar and ETA display
          if ishandle(wb)
            waitbar(index/numel(self.images), wb, [ 'Stack: ' this_img.id ' (close to abort)...' ]);
            try;
            set(wb, 'Name', [ num2str(index) '/' num2str(numel(self.images)) ' ' eta ]);
            end
          else
            disp('Stack: Aborting (user closed the waitbar).')
            break;
          end
        end
      end % for
      disp([ mfilename ': Elapsed time ' num2str(etime(clock, t0)) ' [s]' ])
      disp([ mfilename ': Total exposure on stacked image: ' num2str(ExposureTime) ' [s]' ]);
      delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
      for index=1:size(self.light,3)
        self.light(:,:,index) = self.light(:,:,index) ./ double(self.lightN);
      end
      self.light= im2uint(self.light, 'uint16');
      im        = self.light;
      % display final image
      image(self.light); 
      title([ mfilename ': Stacked image. Exposure=' num2str(ExposureTime) ' [s]' ]);
      % write stacked file
      this_img = self.images(self.reference);
      this_img.exif.ExposureTime = ExposureTime;
      write_stacked(self.light, [ this_img.id '.png' ], this_img.exif);
      
    end % stack
    
    function [ret_t, ret_R] = diff(self, img1)
      [~,img1] = imread(self, img1, 0);
      [~,img2] = imread(self, self.reference, 0);
      if numel(img1.points.x) < 2
        img1 = cpselect(self, img1);
      end
      if numel(img2.points.x) < 2
        img2 = cpselect(self, img2);
      end
      ret_t= []; ret_R = []
      if isempty(img2), return; end
      [ret_t, ret_R] = imdiff(self, img2, img1);
    end % diff
    
    function label(self, lab, img)
      if nargin < 2, return; end
      if nargin < 3, img = 1:numel(self.images); end
      
      if numel(img) == 1
        % mark a single image
        if img > 0 && img <= numel(self.images)
          self.images(img).type = lab;
          disp([ mfilename ': ' self.images(img).id ' marked as ' upper(lab) ]);
        end
      else
        selection = listdlg(self, [ 'select images to mark as ' upper(lab) ]);
        if isempty(selection), return; end
        for index=selection(:)'
          label(self, lab, img(index));
        end
      end
    end % label
    
    function selection = listdlg(self, name, SelectionMode)
      % listdlg: display a dialog list to select images
      if nargin < 2, name          = 'select image'; end
      if nargin < 3, SelectionMode = 'multiple'; end
      liststring = {};
      for index=1:numel(self.images)
        this = self.images(index);
        if ischar(this.source)
          liststring{end+1} = this.source;
        else
          liststring{end+1} = this.id;
        end
        % we add type sharpness
        liststring{end} = [ liststring{end} ' [' num2str(index) '] ' upper(this.type) ...
          ' sharp=' num2str(this.sharpness) ];
        if isfield(this.points,'sx') && ~isempty(this.points.sx) && ~isempty(this.points.sy)
          liststring{end} = [ liststring{end} ...
          ' width=' num2str(this.width) ];
        end
      end
      if isempty(liststring), selection=[]; return; end
      [selection, ok] = listdlg('ListString', liststring, ...
        'ListSize', [480 640], ...
        'SelectionMode', SelectionMode, ...
        'Name', [ mfilename ': ' name ]);
    end % listdlg
    
    function this_img = cpselect(self, img, im)
      % cpselect: automatically set control points
      if nargin < 2, img=1:numel(self.images); end
      
      delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
      if numel(img) > 1
        wb  = waitbar(0, [ mfilename ': Aligning images (close to abort)...' ]); 
        set(wb, 'Tag', [ mfilename '_waitbar' ]);
      else wb = [];
      end
      t0 = clock; 
      
      for index=1:numel(img)
        this_img = img(index);
        % get the image
        if nargin < 3 || isempty(im) || ~isnumeric(im)
          [im, this_img] = imread(self, this_img, 1);
        else
          [~, this_img]  = imread(self, this_img, 0);
        end
        
        % update waitbar and ETA display
        if numel(img) > 1 && ~isempty(wb)
          if ~ishandle(wb)
            disp('Align: Aborting (user closed the waitbar).')
            break;
          end
          % compute ETA
          dt_from0     = etime(clock, t0);
          dt_per_image = dt_from0/index;
          % remaining images: numel(s)-index
          eta    = dt_per_image*(numel(img)-index+1);
          ending = addtodate(now, ceil(eta), 'second');
          ending = [ 'Ending ' datestr(ending) ];
          eta    = sprintf('ETA %i [s]. %s', round(eta), ending);
          waitbar(index/numel(img), wb, [ 'Align: ' this_img.id ' (close to abort)...' ]);
          try
          set(wb, 'Name', [ num2str(index) '/' num2str(numel(img)) ' ' eta ]);
          end
        end

        this_img.points = ...
            find_control_points(rgb2gray(im), 20, self.toleranceTranslation*10);
        this_img.width  = ...
            sqrt(sum(this_img.points.sx.^2.*this_img.points.sy.^2)) ...
            /numel(this_img.points.sx);
        this_img.sharpness  = ...
            sqrt(sum(this_img.points.sharpness.^2)) ...
            /numel(this_img.points.sharpness);
        self.images(this_img.index) = this_img;
      end
      
      delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
      
    end % cpselect
  
  end % methods
  
end % mastrostack

% -----------------------------------------------------------------------------
function ButtonDownCallback(src, evnt, self)
  % ButtonDownCallback: callback when user clicks on the image

  fig = self.figure;
  if ~ishandle(fig), return; end
  
  if strcmp(get(self.figure, 'SelectionType'),'normal') % left click adds a control point
    
    xy = get(gca, 'CurrentPoint'); 
    y = xy(1,1); x = xy(1,2);

    % get the image data
    im = findall(gca, 'Type','image');
    im = get(im, 'CData'); im = rgb2gray(im);
    dx = self.toleranceTranslation*size(im,1)*2;
    dy = self.toleranceTranslation*size(im,2)*2;

    % extract local image
    dx1 = round((x-dx):(x+dx));
    dy1 = round((y-dy):(y+dy));
    try
      im1 = im(dx1, dy1); clear im
    catch
      return
    end

    % get closest max
    [m0,x0]=max(im1(:));
    [x,y]=ind2sub(size(im1),x0);
    [sx,sy, fx, fy] = peakwidth(im1, x,y);  % refine around max
    x=x+min(dx1); y=y+min(dy1);
    self.images(self.currentImage).points.x(end+1) = x;
    self.images(self.currentImage).points.y(end+1) = y;
    self.images(self.currentImage).points.m(end+1) = m0;
    self.images(self.currentImage).points.sx(end+1) = sx;
    self.images(self.currentImage).points.sy(end+1) = sy;
  elseif strcmp(get(self.figure, 'SelectionType'),'alt') % right -> remove last
    if numel(self.images(self.currentImage).points.x)
      self.images(self.currentImage).points.x(end) = [];
      self.images(self.currentImage).points.y(end) = [];
      self.images(self.currentImage).points.m(end) = [];
      self.images(self.currentImage).points.sx(end) = [];
      self.images(self.currentImage).points.sy(end) = [];
    end
  end
  
  % plot control points
  if ishandle(self.images(self.currentImage).points.handle)
    delete(self.images(self.currentImage).points.handle);
  end
  hold on
  self.images(self.currentImage).points.handle = scatter(...
    self.images(self.currentImage).points.y, ...
    self.images(self.currentImage).points.x, 100, 'g', 'x');
  hold off
  
end % ButtonDownCallback

function ScrollWheelCallback(src, evnt, self)
  % ScrollWheelCallback: callback to change image (back/forth) with mouse wheel
  
  % we just replot next/previous image
  index=self.currentImage+sign(evnt.VerticalScrollCount);
  index=max(index,1); index=min(index,numel(self.images));
  if index ~= self.currentImage
    plot(self, index);
  end

end % ScrollWheelCallback

function MenuCallback(src, evnt, self)

  if nargin < 3, self=[]; end
  
  if ~isempty(self) && isempty(self.dndcontrol)
    % initialize Java hooks and activate Drag-n-Drop from external source (files, text)
    if ~(exist('MLDropTarget','class') == 8)
      dndcontrol.initJava;
    end
    target = uicontrol('style','pushbutton','String','Drop Files Here', ...
      'Units','normalized', 'Position', [ 0 0 0.3 0.07 ], ...
      'callback', {@MenuCallback, self });
    self.dndcontrol = dndcontrol(target,@MenuCallback,@MenuCallback);
  end

  if ishandle(src) && isempty(evnt) % menu callback
    if strcmp(get(src,'type'), 'figure')
      action = 'Close';
    elseif strcmp(get(src,'type'), 'uicontrol')
      self.load('');
      return
    elseif strcmp(get(src,'type'), 'uimenu')
      action = get(src,'Label');
    end
  elseif isfield(evnt, 'Key')
    action = upper(evnt.Key);
  elseif isfield(evnt, 'DropType')  % drag-n-drop from dndcontrol
    self = get(gcf, 'UserData');
    switch evnt.DropType
    case 'file'
        for n = 1:numel(evnt.Data)
            self.load(deblank(evnt.Data{n}));
        end
    case 'string'
        lines = textscan(evnt.Data, '%s','Delimiter',sprintf('\r\n'));
        if isempty(lines), retuen; end
        lines = lines{1};
        for n = 1:numel(lines)
          this = deblank(lines{n});
          if ~isempty(this), self.load(this); end
        end
    end
    return
  end

  switch action
  case {'LEFTARROW','PAGEUP','UPARROW','P','Previous image'}   % previous
    if isempty(self.currentImage) || self.currentImage < 2
      self.currentImage = 2;
    end
    plot(self, self.currentImage-1);
  case {'RIGHTARROW','PAGEDOWN','DOWNARROW','N','Next image'}  % next
    if isempty(self.currentImage) || self.currentImage>=numel(self.images)
      self.currentImage = numel(self.images)-1;
    end
    plot(self, self.currentImage+1);
  case {'L'}  % Light
    label(self, 'light', self.currentImage);
  case {'D'}  % Dark
    label(self, 'dark', self.currentImage);
  case {'F'}  % Flat
    label(self, 'flat', self.currentImage);
  case {'I'}
    label(self, 'skip', self.currentImage);
    if self.currentImage<numel(self.images)
      plot(self, self.currentImage+1);
    end
  case {'Mark as Light image...'}
    label(self, 'light');
  case {'Mark as Dark image...'}  % Dark
    label(self, 'dark');
  case {'Mark as Flat image...'}  % Flat
    label(self, 'flat');
  case {'Mark as "skip"...'}
    label(self, 'skip');
  case {'R','Mark as Reference image'}  % Reference
    if ~isempty(self.currentImage) && ...
      numel(self.images) >= self.currentImage && self.currentImage > 0
      if self.reference && self.reference ~= self.currentImage;
        self.images(self.reference).type = 'light';
        disp([ mfilename ': image ' self.images(self.reference).id ' set as LIGHT (was REFERENCE)' ]);
      end
      self.images(self.currentImage).type = 'reference';
      self.reference = self.currentImage;
      disp([ mfilename ': image ' self.images(self.currentImage).id ' set as REFERENCE' ]);
      title([ self.images(self.currentImage).id ...
          ' ' self.images(self.currentImage).type ]);
    end
  case {'C','Clear all control points'}  % clear points
    if ~isempty(self.currentImage) && ...
      numel(self.images) >= self.currentImage && self.currentImage > 0
      self.images(self.currentImage).points.x = [];
      self.images(self.currentImage).points.y = [];
      self.images(self.currentImage).points.sx= [];
      self.images(self.currentImage).points.sy= [];
      self.images(self.currentImage).points.m = [];
      if ishandle(self.images(self.currentImage).points.handle)
        delete(self.images(self.currentImage).points.handle);
      end
    end
  case {'Automatic control points','A'}
    % search points=(x,y,m,sx,sy) and store
    if ~isempty(self.currentImage) && ...
      numel(self.images) >= self.currentImage && self.currentImage > 0
      % get the image data
      im = findall(gca, 'Type','image');
      im = get(im, 'CData');
      if ~isempty(im)
        cpselect(self, self.currentImage, im);
        plot(self, self.currentImage, im);
      end
    end
  case {'Open','O'}
    [~,img] = imread(self,'',0);
  case {'Open light images','O'}
    [~,img] = imread(self,'',0);
    % label them as light
    for index=1:numel(img)
      this_img = img(index);
      self.images(this_img.index).type = 'light';
    end
  case {'Open dark images'}
    [~,img] = imread(self,'',0);
    % label them as dark
    for index=1:numel(img)
      this_img = img(index);
      self.images(this_img.index).type = 'dark';
    end
  case {'Open flat images'}
    [~,img] = imread(self,'',0);
    % label them as dark
    for index=1:numel(img)
      this_img = img(index);
      self.images(this_img.index).type = 'flat';
    end
  case {'Show in log scale'}
    if ~isempty(self.currentImage) && ...
      numel(self.images) >= self.currentImage && self.currentImage > 0
      plot(self, self.currentImage, 'log');
    end
  case 'Compute master Dark'
    self.dark = [];
    image(im2uint(self.getdark));
    title('Master Dark')
  case 'Compute master Flat'
    self.flat = [];
    im = self.getflat;
    image(im2uint(im));
    title([ 'Master Flat ' mat2str([ min(im(:)) max(im(:)) ]) ]);
  case 'Stack'
    stack(self);
  case 'Align'
    cpselect(self);
  case 'RETURN'
    if isempty(self.currentImage) || self.currentImage < 1 ...
      || self.currentImage>numel(self.images)
      self.currentImage = 1;
    end
    plot(self, self.currentImage);
  case 'Goto image...'
    selection = listdlg(self, 'select image to go to','single');
    if isempty(selection), return; end
    plot(self, selection);
  case 'Clear image(s)...'
    selection = listdlg(self, 'select image(s) to clear');
    if isempty(selection), return; end
    self.images(selection) = [];
    self.reference    = 0;
    for index=1:numel(self.images)
      self.images(index).index = index;
      if strcmp(self.images(index).type, 'reference')
        self.reference = index;
      end    
    end
    self.currentImage = [];
  case 'Close'
    delete(self.dndcontrol);
    self.dndcontrol = [];
    delete(gcf);
  case 'Show sharpness'
    figure; plot([ self.images.sharpness ]./[ self.images.width ])
    xlabel('Image index')
    ylabel('Sharpness (higher is better)');
    titme('Sharpness: identify low sharpness with the Marker tool, then mark them as Skip');
  case 'Set tolerances'
    prompt={'\bf {\color{blue}Translation} (0-1, e.g 0.01):','\bf {\color{blue}Rotation} [deg, e.g. 3]:'};
    name=[ mfilename ': Tolerances' ];
    numlines=1;
    defaultanswer={ num2str(self.toleranceTranslation), ...
                    num2str(self.toleranceRotation) };
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    if isempty(answer), return; end
    if isfinite(str2double(answer{1}))
      self.toleranceTranslation = str2double(answer{1});
    end
    if isfinite(str2double(answer{2}))
      self.toleranceRotation = str2double(answer{2});
    end
  otherwise
    disp([ mfilename ': unknown action ' action ])
  end
end % MenuCallback
