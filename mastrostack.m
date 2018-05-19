classdef mastrostack < handle
  % mastrostack: a class to automatically align and stack astro-photography images
  %
  %  Purpose
  %
  %  This class gets a list of images, and automatically determines bright stars as
  %  control points. These are followed along pictures, and used to build an affine
  %  transformation at constant scale (e.g. a rotation and translation). All images
  %  are then stacked. The images can be given as file names (may include wildcards),
  %  or matrices from e.g. imread, and support both RGB and gray images. As stars are
  %  used for the alignment, this method is suited for deep sky images, but not for
  %  planetary imaging.
  %
  %  This function does not make use of the phase correlation technique, but operates
  %  directly on the images. It assumes that at least two bright stars remain on each
  %  picture, but tolerate translations and rotations, as well as stars/control
  %  points that appear/disappear on the field of view. The procedure also ignores
  %  very sharp peaks, assuming they originate from dead pixels (and would then be
  %  static over images, leading to traces).
  %
  %  It is highly recommended to specify a dark frame filename, which will be
  %  subtracted to all images, in order to e.g. remove most hot/dead pixels. To get
  %  such a 'dark', use your camera alone and shot once with the cap on, same
  %  settings as the other pictures (duration, ISO). Image should be full black.
  %
  %  You may as well specify a flat frame filename, which will be divided to all
  %  images, in order to counter vignetting (darker borders). To get such a 'flat',
  %  shot once with the scope pointing at a uniform view (sky, white wall). Adapt the
  %  duration so that you get a rather gray image (not full black or white).
  %
  % Syntax
  % -----------
  %   ma = mastrostack;
  %     start the user interface, without loading images (can be done afterwards)
  %   ma = mastrostack(images)
  %     loads images without setting their type
  %   ma = mastrostack(light, dark)
  %   ma = mastrostack(light, dark, flat)
  %     loads light, dark (background) and flat (scope response) images, and label them.
  %    
  %  Importing images
  %  ----------------
  %  Start with:
  %
  %    ma=mastrostack;
  %    
  %  Then press the Return key on the main interface. A Drop Files Here button
  %  appears in the lower left side. Drag and drop your Dark, Flat and Light images
  %  there. Images having 'dark' or 'flat' in their path/file name are marked as such
  %  automatically. You may alternatively use the File menu items.
  %
  %  Supported image formats include JPG, PNG, TIFF, FITS. If you have installed
  %  readraw, you may as well directly import RAW camera images. This is highly
  %  recommended, as it retains much more information from the camera shot than the
  %  generated JPEG images, which prooves to be essential for subtracting the Dark
  %  image (background), and revealing faint objects.
  %
  %  Preparing the Stacking
  %  ----------------------
  %  After importing the files, you may label them using the 'Image/Mark as...'
  %  menu items. You can navigate within images with the Image/Goto menu item, and
  %  the arrow keys, or the mouse wheel. 'Bad' images can be skipped (ignored). To
  %  use them back, set their type to 'light'. You should then compute the master
  %  Dark and Flat images (Compute menu).
  %
  %  It is recommended to zoom onto specific features (e.g. a set of stars) to check
  %  visually for their sharpness. Deselect the Zoom tool, and scan through images
  %  using the left arrow key, and press the 'I' key to mark images to be ignored,
  %  such as those blurred. To reset the plot, press the Return key.
  %
  %  You can select the Reference image, which will be used as template for stacking.
  %  If not defined, the first image in the list will be used as such when stacking.
  %
  %  Optionally, use the Compute/Align item to compute the images control points
  %  (stars) which also corrects for background and scope response when dark and flat
  %  are defined. This procedure computes a metric to quantify the sharpness. You may
  %  then select the 'Compute/Select on sharpness' menu item which metric is highest for
  %  clearer images, and low for blurred/moved ones. The axis is the image index. To
  %  select best images, use the left-click to set label as 'skip', right-click to 
  %  set label as 'light', shift-click to open that image and any ESC to
  %  abort.
  %
  %  Stacking
  %  --------
  %  When ready, use the Compute/Stack menu item. If the Alignment has not been
  %  executed previously, it is achieved for each image. The final image is then
  %  shown and written to disk. Use e.g. Lightroom, RawTherapee, DarkTable to enhance
  %  contrast.
  %
  %  Notes
  %  -----
  %  If you have difficulties in stacking (some images do not have enough control
  %  points), relax e.g. the translation tolerance, using the menu item 'Compute/Set
  %  tolerances'. You can also increase the number of control points. In case the
  %  main interface is closed, get it back with
  %
  %    plot(ma)
  %    
  %  Using commands (scripting)
  %  --------------------------
  %  Using commands allow to script/automate the procedure.
  %
  %     % create Ma(e)stroStack and import images
  %     ma=mastrostack('path/to/images/*.JPG','path/to/darks/*.JPG','path/to/flats/*.JPG');
  %
  %     % stack. The first 'light' image will be used as Reference for stacking
  %     stack(ma);
  %
  %
  % Methods
  % -------
  % about           display a dialg box
  % correct         correct and image for dark (background) and flat (vigneting)
  % cpselect        ALIGN images on reference
  % diff            compute difference of an image with reference
  % delete          delete the mastrostack, and clear memory
  % exist           check if an image has already been loaded            
  % getdark         compute the master Dark
  % getflat         compute the master Flat
  % help            open the help page
  % imread          load an image and return its matrix and information
  % label           label an image as light, dark, flat or skip
  % listdlg         display a selector for images
  % load            load an image and return its information
  % mastrostack     create the Ma(e)strostack session
  % plot            plot the user interface and images
  % stack           STACK light images, correcting with dark and flat.
  %
  % (c) E. Farhi, 2018. GPL2.

  properties
    toleranceTranslation  = 0.01; % in percent
    toleranceRotation     = 1;    % in deg
    deadPixelArea         = 9;
    reference             = 0;    % by default the first 'light' image
    figure                = [];
    images = [];
    %      id
    %      source (or matrix if directly given)
    %      image_size
    %      image_sum
    %      exif
    %      type (light, dark, flat, skip, none)
    %      points
    %      rotation
    %      translation
    %      sharpness
    %      width
    %      thumbnail
    
    % the actual image matrices are not stored (except when given as such)
    % as we expect to handle 100's of images.
    dark    = 0;  % stored as a double in [0-1] uint16
    flat    = 0;  % stored as a double in [0-1] float
    light   = 0;  % stored as a double in [0-1]
    lightN  = 0;  % stored as a uint16
    nbControlPoints = 50;
    currentImage = 0;  % index of the current image on figure
    dndcontrol = [];
  end % properties
  
  methods
  
    % I/O stuff ================================================================
    function self = mastrostack(source, dark, flat)
      % mastrostack: create the Ma(e)strostack session, load files and display
      %   the interface.
      %
      % self = mastrostack(images)
      %   import given images, specified as filenames (supports directories, wildcards)
      % self = mastrostack(light, dark, flat)
      %   import specified 'light','dark','flat' images.

      % load catalogs: objects, stars
      disp([ mfilename ': Welcome !' ]);
      
      if nargin == 1 % import images
        img = load(self, source);
      elseif nargin > 1
        img = load(self, source);
        label(self, 'light', img);
        img = load(self, dark);
        label(self, 'dark', img);
        if narhin > 2
          img = load(self, flat);
          label(self, 'flat', img);
        end
      end
      
      plot(self);
      
    end % mastrostack
    
    function [im, img] = imread(self, source, flag)
      % imread: read an image, and store its information
      %
      % Supported file formats include JPG, PNG, BMP, GIF, TIFF, FITS
      %
      %  [im, img] = imread(self, filename or matrix or index)
      %    returns the image matrix and information
      %  [~, img] = imread(self, filename or matrix or index, 0)
      %    only returns the image information, avoiding to read the file
      %
      % see also: imread
      
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
        
        % set image index
        if isnan(found(index))
          this_img.index = numel(self.images)+1;
        else % update
          this_img.index = found(index);
        end
        % store
        if isempty(self.images) && this_img.index > 0
          self.images    = this_img;
        elseif this_img.index > 0
          self.images(this_img.index) = this_img; 
        else
          return
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
    
    function img = load(self, source)
      % load: read an image, and store its information
      %
      % Supported file formats include JPG, PNG, BMP, GIF, TIFF, FITS
      %
      %  img = load(self, filename or matrix or index)
      %    returns the image matrix and information.
      [~,img] = imread(self, source, 0);
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
    
    function save(self, filename)
      % save: save the maestrostack session
      if nargin < 2, filename = 'mastrostack_session.mat'; end
      builtin(@save, 'self', filename);
      disp([ mfilename ': saved as ' filename ]);
    end 
    
    % compute stuff ============================================================
    
    function im = getdark(self)
      % getdark: compute the mean Dark frame, using labelled images
      %
      % im = getdark(self)
      %   return the master flat image as the mean of all 'dark' labelled images
      
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
          self.dark = im2uint(sumDarks/nbDarks, 'uint16');
          clear sumDarks
          disp([ mfilename ':   used ' num2str(nbDarks) ' images.' ])
        else self.dark=0;
        end
      end
      if nargout
        im = self.dark;
      end
      
    end % getdark
    
    function im = getflat(self)
      % getflat: compute the mean Flat frame, using labelled images
      %
      % im = getflat(self)
      %   return the master flat image as the mean of all 'flat' labelled images
      
      if any(~isfinite(self.flat)), self.flat=0; end
      if isempty(self.flat) || isscalar(self.flat)
        disp([ mfilename ': computing master Flat' ])
        sumFlats = 0;
        nbFlats  = 0;
        for index=1:numel(self.images)
          img = self.images(index);
          if ~strcmp(img.type, 'flat'), continue; end
          % read image (re-read if needed), convert to uint16
          im = im2uint(rgb2gray(imread(self, index)), 'uint16');
          % remove the dark (uint16) if any is available
          if ~isempty(self.dark) && ~isscalar(self.dark)
            im = im - rgb2gray(self.dark);
          end
          im = imdouble(im);
          im = im/max(im(:));
          sumFlats  = sumFlats+ im;
          clear im
          nbFlats   = nbFlats + 1;
        end
        if nbFlats>0
          self.flat = single(sumFlats/nbFlats);
          clear sumFlats
          index     = isfinite(self.flat) & self.flat > 0;
          self.flat = self.flat/max(max(self.flat(index)));
          index     = ~isfinite(self.flat) | self.flat <= 0;
          self.flat(index) = 1;
          clear index
          disp([ mfilename ':   used ' num2str(nbFlats) ' images.' ])
        else self.flat=0; end
      end
      if nargout
        im = self.flat;
      end
      
    end % getflat

    function im = correct(self, im, cl)
      % correct: correct an image for vigneting and background
      %
      % correct(self, im)
      %   divide the image by the flat, and subtract the dark
      %   'im' is a MxN or MxNx3 matrix, from e.g. imread
      if nargin < 2, return; end
      if nargin < 3
        cl = class(im); % we shall then cast back
      end
      
      if isempty(im), return; end
      im = im2uint(im, 'uint16');
      
      % correct for read-out/sensor noise (subtract)
      if ~isempty(self.dark) && ~isscalar(self.dark) && ndims(self.dark) == 3
        im = im - self.dark; % dark is a uint16
      end
      
      % correct for flat field to compensate vigneting (divide)
      if ~isempty(self.flat) && ~isscalar(self.flat)
        for l=1:size(im,3)
          im(:,:,l) = im2uint(imdouble(im(:,:,l))./double(self.flat),'uint16');
        end
      end

      im = im2uint(im, cl);
      
    end % correct
    
    function t=about(self)
      % about: display an About dialog box
      %
      % t=about(self)
      %   returns the 'about' string
      
      t = {'\bf Welcome to {\color{magenta}Ma(e)stroStack} ! (c) E.Farhi', ...
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
      if nargout == 0
        CreateMode.WindowStyle = 'modal';
        CreateMode.Interpreter='tex';
        msgbox(t, ...
          'Ma(e)stroStack: About', CreateMode);
      end
    end % about
    
    function url=help(self)
      % help(self): open the Help page
      url = fullfile('file:///',fileparts(which(mfilename)),'doc','mastrostack.html');
      open_system_browser(url);
    end
    
    function h=plot(self, img, im)
      % plot: plot images
      %
      % h=plot(self)
      %   plot the current or first 'light' image
      % h=plot(self, image)
      %   same as above with specified image (can be given as name or index)
      if nargin < 2
        % we search for the first 'light' image
        img = [];
        for index=1:numel(self.images)
          if strcmp(self.images(index).type, 'light') || isempty(self.images(index).type)
            img = self.images(index);
            break
          end
        end
        if isempty(img), img = 1; end
      end
      
      fig = findall(0, 'Tag', 'mastrostack');
      if isempty(fig) % build the figure
        fig = build_interface(self);
        
        t = text(0.1,0.5, about(self), 'Units','normalized');
        xl=[]; yl=[];
      else
        set(0, 'CurrentFigure',fig); % activate but no raise
        % get current xlim/ylim (to retain any zoom)
        xl = xlim; yl=ylim;
      end
      self.figure = fig;
      
      if nargin >= 3 && ischar(im)
        option = im;
        im     = [];
      else option = '';
      end

      % get the image
      if nargin < 3 || isempty(im) || ~isnumeric(im)
        [im, img] = imread(self, img);
      else
        [~, img] = imread(self, img, 0);
      end
      if isempty(img) || isempty(im), return; end
      self.currentImage = img.index;
      t = [ self.images(self.currentImage).id ...
        ' ' self.images(self.currentImage).type ];
      if self.images(self.currentImage).sharpness > 0
        t = [ t ' sharp=' num2str(self.images(self.currentImage).sharpness) ];
      end
      if self.images(self.currentImage).width > 0
        t = [ t ' width=' num2str(self.images(self.currentImage).width) ];
      end
      if self.images(self.currentImage).intensity > 0
        t = [ t ' int=' num2str(self.images(self.currentImage).intensity) ];
      end
      if strcmp(self.images(self.currentImage).type, 'light') && ...
        ((~isempty(self.flat) && ~isscalar(self.flat)) || ...
         (~isempty(self.dark) && ~isscalar(self.dark)))
        im = correct(self, im);
        t = [ t ' (-dark, /flat)' ];
      end
      t = [ t ' [' num2str(self.currentImage) '/' num2str(numel(self.images)) ']' ];
      
      im = im2uint(im,'uint8');
      if strcmp(option, 'log')
        im = imlogscale(im);
      end
      h = image(im);
      if ~isempty(xl)
        xlim(xl); ylim(yl);
      end
      title(t);

      % display control points
      hold on
      self.images(self.currentImage).points.handle = scatter(...
        self.images(self.currentImage).points.y, ...
        self.images(self.currentImage).points.x, 100, 'g', 'o');
      hold off
      
    end % plot
    
    function [im, filename] = stack(self, filename)
      % stack: stack all images on the reference
      %
      % im = stack(self)
      %   return the stacked image, using the master dark and flat for correction,
      %   aligning all images on the reference, and adding them with sub-pixel
      %   resolution. The resulting image is written to disk at the reference
      %   location with name '_stacked.png'. 
      %   You can stop the stacking by closing the progress bar. 
      %   If some images are ignored, you may increase e.g. toleranceTranslation.
      % im = stack(self, filename)
      %   same as above, and specifies an output filename.
      %
      % The stacking does sum of (light - dark)/(flat - dark)
      
      % get dark and flat
      getdark(self);
      getflat(self);

      self.light = 0.0;
      self.lightN= uint8(0);
      ExposureTime = 0;
      if nargin < 2, filename = ''; end
      
      nbtostack = 0;
      
      % we first check for the reference and its control points
      for index=1:numel(self.images)
        this_img = self.images(index);
    
        if isempty(this_img.type) || strcmp('light', this_img.type)
          if ~self.reference || self.reference == index
            self.reference = index;
            this_img       = self.images(self.reference);
            self.images(self.reference).type = 'reference';
            disp([ mfilename ':   using reference as ' this_img.id ' (' num2str(index) ')' ])
            im = correct(self, imread(self, this_img)); % read and correct image
            if isempty(this_img.points) || ~isstruct(this_img.points) || ...
              ~isfield(this_img.points, 'x') || numel(this_img.points.x) < 2
              % compute control points
              disp([ mfilename ':   computing control points for ' this_img.id ' (' num2str(index) ') REF' ])
              [this_img, self] = cpselect(self, this_img, im);
              % self.images(self.reference) = this_img;
            end
            self.light  = imdouble(im);
            clear im
            self.lightN = self.lightN+1;
            if isfield(this_img.exif, 'ExposureTime')
              ExposureTime = ExposureTime + this_img.exif.ExposureTime;
            end
          end
          nbtostack = nbtostack+1;
        end
      end
      if self.reference < 1 || self.reference > numel(self.images), return; end
      
      delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
      wb  = waitbar(0, [ mfilename ': Stacking ' num2str(nbtostack) ' images (close to abort)...' ]); 
      set(wb, 'Tag', [ mfilename '_waitbar' ]);
      t0 = clock; stackedimages = 0;
      disp([ mfilename ':   Stacking ' num2str(nbtostack) ' images... ' datestr(now) ])
      
      for index=1:numel(self.images)
        this_img = self.images(index);
        
        if isempty(this_img.type) || strcmp('light', this_img.type)
          stackedimages= stackedimages+1;
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
            [this_img, self] = cpselect(self, this_img, im);
            % self.images(index) = this_img;
          end
          
          % display
          plot(self, this_img, im);
          
          % compute the affine transformation wrt reference (using control points)
          [ret_t, ret_R] = diff(self, this_img);
          if isempty (ret_t), continue; end
          
          % we read the image and transform it
          [im,M]     = imaffine(im, ret_R, ret_t);
          self.lightN= self.lightN+M;
          clear M
          self.light = self.light + imdouble(im); % add rotated/translated image on the reference
          clear im
          
          % sum-up exposure
          if isfield(this_img.exif, 'ExposureTime')
            ExposureTime = ExposureTime + this_img.exif.ExposureTime;
          end
          
          % compute ETA
          dt_from0     = etime(clock, t0);
          dt_per_image = dt_from0/stackedimages;
          % remaining images: numel(s)-index
          eta    = dt_per_image*(nbtostack-stackedimages);
          ending = addtodate(now, ceil(eta), 'second');
          ending = [ 'Ending ' datestr(ending) ];
          eta    = sprintf('ETA %i [s]. %s', round(eta), ending);
          disp([ mfilename ': ' ...
              num2str(stackedimages) '/' num2str(nbtostack) ': ' ...
              this_img.id ' [' num2str(this_img.index) ']. ' eta]);
          
          % update waitbar and ETA display
          if ishandle(wb)
            waitbar(stackedimages/nbtostack, wb, [ 'Stack: ' this_img.id ' (close to abort)...' ]);
            try
            set(wb, 'Name', [ num2str(stackedimages) '/' num2str(nbtostack) ' ' eta ]);
            end
          else
            disp('Stack: Aborting (user closed the waitbar).')
            break;
          end
        end
      end % for
      
      disp([ mfilename ': Total exposure on stacked image: ' num2str(ExposureTime) ' [s]' ]);
      delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
      for index=1:size(self.light,3)
        self.light(:,:,index) = self.light(:,:,index) ./ double(self.lightN);
      end
      self.light= im2uint(self.light, 'uint16');
      
      % display final image
      image(self.light); 
      title([ mfilename ': Stacked image. Exposure=' num2str(ExposureTime) ' [s]' ]);
      % write stacked file
      this_img = self.images(self.reference);
      this_img.exif.ExposureTime = ExposureTime;
      if isempty(filename)
        if ischar(this_img.source)
          p = fileparts(this_img.source);
        else p = pwd; end
        filename = fullfile(p, [ this_img.id '_stacked.png' ]);
      end
      write_stacked(self.light, filename, this_img.exif);
      disp([ mfilename ': Stacking DONE. Elapsed time ' num2str(etime(clock, t0)) ' [s] ' datestr(now) ])
      if nargout
        im        = self.light;
      end
    end % stack
    
    function [ret_t, ret_R, self] = diff(self, img1)
      % diff: compute the translation and rotation wrt reference image
      %
      % [t,R] = diff(self, image)
      %   return the translation and rotation of image wrt reference
      [~,img1] = imread(self, img1, 0);
      [~,img2] = imread(self, self.reference, 0);
      ret_t= []; ret_R = [];
      if isempty(img1) || isempty(img2), return; end
      if numel(img1.points.x) < 2
        [img1, self] = cpselect(self, img1);
      end
      if numel(img2.points.x) < 2
        [img2, self] = cpselect(self, img2);
      end
      
      if isempty(img2), return; end
      [ret_t, ret_R] = imdiff(self, img2, img1);
    end % diff
    
    function label(self, lab, img)
      % label: set the label on given/all data sets
      %
      % supported labels are '' (none), 'light','dark','flat',skip'
      % 
      % label(self, label)
      %   set the same label on all images
      % label(self, label, images)
      %   set the same label on on given images
      %   'images' can be given as name or index
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
      %
      % listdlg(self)
      %   display a dialog to select some images. return the selection indices
      % listdlg(self, title)
      %   display a dialog to select some images with given dialog title
      % listdlg(self, title, 'single' or 'multiple')
      %   display a dialog with title, and single or multiple selection ability

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
    
    function [this_img, self] = align(self, img, im)
      [this_img, self] = cpselect(self, img, im);
    end
    
    function [nimg, self] = cpselect(self, img, im)
      % cpselect: automatically set control points
      %
      % cpselect(self)
      %   set all control points and compute sharpness
      % cpselect(self, images)
      %   set control points and compute sharpness on given images
      %   'images' can be given as name or index
      
      if nargin < 2, img=1:numel(self.images); end
      
      if numel(img) > 1
        delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
        wb  = waitbar(0, [ mfilename ': Aligning images (close to abort)...' ]); 
        set(wb, 'Tag', [ mfilename '_waitbar' ]);
        disp([ mfilename ':   Aligning... ' datestr(now) ])
      else wb = [];
      end
      t0 = clock; 
      
      nimg = [];
      
      for index=1:numel(img)
        this_img = img(index);
        
        % get the image
        [~, this_img]  = imread(self, this_img, 0);
        if isempty(this_img.points.x) || nargin < 2
        
          if nargin < 3 || isempty(im) || ~isnumeric(im)
            [im, this_img] = imread(self, this_img, 1);
          end
          im = rgb2gray(im);
          
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
              find_control_points(im, self.nbControlPoints, self.deadPixelArea);
          this_img.width  = ...
              sqrt(sum(this_img.points.sx.^2.*this_img.points.sy.^2)) ...
              /numel(this_img.points.sx);
          this_img.sharpness  = ...
              sqrt(sum(this_img.points.sharpness.^2)) ...
              /numel(this_img.points.sharpness);
          this_img.intensity  = ...
              sum(this_img.points.m) ...
              /numel(this_img.points.m);
        end
        
        if isempty(nimg), nimg = this_img;
        else nimg = [ nimg this_img]; end
        self.images(this_img.index) = this_img;
      end % for
      if numel(img) > 1 && ~isempty(wb)
        delete(findall(0, 'Tag', [ mfilename '_waitbar' ]));
        disp([ mfilename ':   Alignment done ' datestr(now) ])
      end
      
    end % cpselect
  
  end % methods
  
end % mastrostack

% -----------------------------------------------------------------------------
function ButtonDownCallback(src, evnt, self)
  % ButtonDownCallback: callback when user clicks on the image

  fig = self.figure;
  if ~ishandle(fig), return; end
  
  if isempty(self.currentImage) || ...
      self.currentImage < 1 || self.currentImage > numel(self.images)
      return
  end
  
  if strcmp(get(self.figure, 'SelectionType'),'normal') % left click adds a control point
    
    xy = get(gca, 'CurrentPoint'); 
    y  = xy(1,1); x = xy(1,2);

    % get the image data
    im = findall(gca, 'Type','image');
    im = get(im, 'CData'); im = rgb2gray(im);
    dx = max(self.toleranceTranslation*size(im));

    % get closest max in sub-image
    [s, f, m] = peakwidth(im, [x y], dx);  % refine around max
    if numel(f) < 2, return; end
    self.images(self.currentImage).points.x(end+1)  = f(1);
    self.images(self.currentImage).points.y(end+1)  = f(2);
    self.images(self.currentImage).points.m(end+1)  = m;
    self.images(self.currentImage).points.sx(end+1) = s(1);
    self.images(self.currentImage).points.sy(end+1) = s(2);
  elseif strcmp(get(self.figure, 'SelectionType'),'alt') % right -> remove last
    if numel(self.images(self.currentImage).points.x)
      self.images(self.currentImage).points.x(end)  = [];
      self.images(self.currentImage).points.y(end)  = [];
      self.images(self.currentImage).points.m(end)  = [];
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
    self.images(self.currentImage).points.x, 100, 'g', 'o');
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
  % MenuCallback: the main callback for menu and keyboard actions
  if nargin < 3, self=[]; end
  if isobject(evnt)
    publicProperties = fieldnames(evnt);
    myStruct = struct();
    for iField = 1:numel(publicProperties)
        myStruct.(publicProperties{iField}) = evnt.(publicProperties{iField}); 
    end
    evnt = myStruct;
  end
  
  
  if ~isempty(self) && isempty(self.dndcontrol)
    % initialize Java hooks and activate Drag-n-Drop from external source (files, text)
    if ~(exist('MLDropTarget','class') == 8)
      dndcontrol.initJava;
    end
    target = uicontrol('style','pushbutton','String','Drop Files Here', ...
      'Units','normalized', 'Position', [ 0 0 0.3 0.07 ], ...
      'ToolTipString','Drop files or click me', ...
      'Tag','mastrostack_DropMeHere', ...
      'callback', {@MenuCallback, self });
    dndcontrol(target,@MenuCallback,@MenuCallback);
    self.dndcontrol = target;
  end

  if isfield(evnt, 'Key')
    action = upper(evnt.Key);
  elseif ishandle(src) % menu callbacksrc
    if strcmp(get(src,'Type'), 'figure')    % clicked Figure Close button
      action = 'Close';
    elseif strcmp(get(src,'Type'), 'uicontrol') % clicked Drop button
      self.load('');
      return
    elseif strcmp(get(src,'Type'), 'uimenu')    % a menu item
      action = get(src,'Label');
    end
  elseif isfield(evnt, 'DropType')  % drag-n-drop from dndcontrol
    self = get(gcf, 'UserData');
    switch evnt.DropType
    case 'file'
      lines = deblank(evnt.Data);
    case 'string'
      lines = textscan(evnt.Data, '%s','Delimiter',sprintf('\r\n'));
      if isempty(lines), return; end
      lines = lines{1};
    end
    
    % now import files
    for n = 1:numel(lines)
      this = deblank(lines{n});
      if ~isempty(this), self.load(this); end
    end
    return
  end

  switch action
  case {'LEFTARROW','PAGEUP','UPARROW','P','Previous image'}   % previous
    plot(self, max(1, self.currentImage-1));
  case {'RIGHTARROW','PAGEDOWN','DOWNARROW','N','Next image',' ','SPACE'}  % next
    plot(self, min(numel(self.images), self.currentImage+1));
  case {'L'}  % Light
    label(self, 'light', self.currentImage);
  case {'D'}  % Dark
    label(self, 'dark', self.currentImage);
  case {'F'}  % Flat
    label(self, 'flat', self.currentImage);
  case {'I','S'}
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
  case {'C','Clear all control points','Clear control points (this image)'}  % clear points
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
    image(self.getdark);
    title('Master Dark')
  case 'Compute master Flat'
    self.flat = [];
    im = im2uint(self.getflat,'uint16');
    image(im);
    title([ 'Master Flat ' mat2str([ min(im(:)) max(im(:)) ]) ]);
  case {'Stack', 'Stack (and Align if needed)' }
    stack(self);
  case 'Align'
    cpselect(self);
  case {'Re-Plot','RETURN'}
    if isempty(self.currentImage) || self.currentImage < 1 ...
      || self.currentImage>numel(self.images)
      self.currentImage = 1;
    end
    plot(self, self.currentImage);
    axis tight
  case {'G','Goto image...'}
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
    self.currentImage = 0;
  case 'Clear skipped image(s)'
    disp([ mfilename ': Removing skipped images:' ]);
    for index=numel(self.images):-1:1
      if strcmp(self.images(index).index, 'skip')
        if ischar(self.images(index).source)
          disp(self.images(index).source)
        end
        self.images(index) = [];
      end
    end
  case 'Close'
    self.dndcontrol = [];
    delete(gcf);
  case 'Set tolerances'
    prompt={'\bf {\color{blue}tolerance on Translation} (in % of image size, 0-1, e.g 0.01):', ...
      '\bf {\color{blue}tolerance on Rotation} [in deg, e.g. 1]:', ...
      '\bf {\color{blue}Nb of Control Points} [e.g. 30-50]:', ...
      '\bf {\color{blue}Dead Pixel Area} [in pixel^2 e.g. 6-9, 0 not to ignore dead pixels]:' };
    name=[ mfilename ': Settings' ];
    numlines=1;
    defaultanswer={ num2str(self.toleranceTranslation), ...
                    num2str(self.toleranceRotation), ...
                    num2str(self.nbControlPoints), ...
                    num2str(self.deadPixelArea) };
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
    if isfinite(str2double(answer{3}))
      self.nbControlPoints = ceil(str2double(answer{3}));
    end
    if isfinite(str2double(answer{4}))
      self.deadPixelArea = str2double(answer{4});
    end
  case 'Select images on sharpness...'
    select_on_sharpness(self);
  case 'About'
    about(self);
  case 'Help'
    help(self);
  case {'Save stacked image','Save master Dark', ...
        'Save master Flat',  'Save current image' }
    exif = []; im = [];
    if strfind(action, 'stacked')
      if self.reference
        this_img = self.images(self.reference);
        exif     = this_img.exif;
      end
      im       = self.light;
    elseif strfind(action, 'Dark')
      im       = self.getdark;
    elseif strfind(action, 'Flat')
      im       = im2uint(self.getflat,'uint16');
    elseif strfind(action, 'current')
      im = findall(gca, 'Type','image');
      im = get(im, 'CData');
    end
    if isempty(im), return; end
    
    [filename, pathname, filterindex] = uiputfile( ...
      {'*.png;*.jpg;*.tiff', 'All image Files (PNG, JPG, TIFF)';
        '*.png',  'PNG image (*.png)'; ...
        '*.jpg',  'JPEG image (*.jpg)'; ...
        '*.tiff', 'TIFF image (*.tiff)'; ...
        '*.*',  'All Files (*.*)'}, ...
        action);
    if isequal(filename,0) || isequal(pathname,0)
      return
    end
    filename = fullfile(pathname, filename);
    
    exif.Comment = action;
    write_stacked(im, filename, exif);
    disp([ mfilename ': ' action ' as ' filename ])
  otherwise
    disp([ mfilename ': unknown action ' action ])
  end
end % MenuCallback

function fig = build_interface(self)
  % build_interface: build the user interface (when not there yet)
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
    uimenu(m, 'Label', 'Save stacked image',        ...
      'Callback', {@MenuCallback, self },'Accelerator','s');
    uimenu(m, 'Label', 'Save master Dark',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Save master Flat',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Save current image',        ...
      'Callback', {@MenuCallback, self });
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
    uimenu(m, 'Label', 'Re-Plot',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Show in log scale',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Goto image...', 'Separator','on',    ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Next image',       ...
      'Callback', {@MenuCallback, self },'Accelerator','n');
    uimenu(m, 'Label', 'Previous image',        ...
      'Callback', {@MenuCallback, self },'Accelerator','p');
    
    uimenu(m, 'Label', 'Clear image(s)...', 'Separator','on',    ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Clear skipped image(s)',  ...
      'Callback', {@MenuCallback, self });
    
    m = uimenu(fig, 'Label', 'Compute');
    uimenu(m, 'Label', 'Set tolerances',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Automatic control points',        ...
      'Callback', {@MenuCallback, self },'Accelerator','a');
    uimenu(m, 'Label', 'Clear control points (this image)',        ...
      'Callback', {@MenuCallback, self },'Accelerator','c');
    uimenu(m, 'Label', 'Select images on sharpness...',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Compute master Dark', 'Separator','on',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Compute master Flat',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Align', 'Separator','on',       ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'Stack (and Align if needed)',       ...
      'Callback', {@MenuCallback, self });
      
    m = uimenu(fig, 'Label', 'Help');
    uimenu(m, 'Label', 'Help',        ...
      'Callback', {@MenuCallback, self });
    uimenu(m, 'Label', 'About',        ...
      'Callback', {@MenuCallback, self });
    
    % install mouse event handler on axis and figure
    set(gca, 'ButtonDownFcn', {@ButtonDownCallback, self }); % will get its CurrentPoint
    set(fig, 'KeyPressFcn',   {@MenuCallback, self });       % will get its CurrentCharacter
    set(fig, 'WindowScrollWheelFcn',   {@ScrollWheelCallback, self });
    set(fig, 'UserData',self);
end % build_interface

function select_on_sharpness(self)

  persistent first_use

  fig=figure('Name', [ mfilename ': Sharpness' ]); 
  if isempty(self.images), close(fig); return; end
  
  [h, x, y, xs, ys] = plot_sharpness(self.images, self.currentImage);
  
  t = { [ 'mastrostack: Sharpness selection' ], ...
    'You may select a rectangle area with the mouse, using:', ...
    '  LEFT       button   selected images are set to SKIP/IGNORE', ...
    '  RIGHT      button   selected images are set to LIGHT', ...
    '  SHIFT-LEFT button   open first image in selection', ...
    '  LEFT/UP    arrow    go to previous image', ...
    '  DOWN/RIGHT arrow    go to next image', ...
    '  G          key      goto image (selection from dialogue)', ...
    '  I/S        key      current image is set to IGNORE/SKIP', ...
    '  L          key      current image is set to LIGHT', ...
    '  X/Q/ESC    key      abort' };
  
  if isempty(first_use)
    first_use = false;
    fprintf(1, '%s\n', t{:});
    % helpdlg(t, [ mfilename ': Sharpness selection Help' ]);
  end
  
  % loop until user press key
  while true
    try
      k = waitforbuttonpress;
    catch
      return % no more active window
    end
    % determines what has been pressed
    key = get(gcf, 'CurrentCharacter');
    p1  = get(gca, 'CurrentPoint');
    but = get(gcf, 'SelectionType'); % = 'normal': left or 'alt': right
    
    % key pressed: test if we exit the loop
    if ~isempty(key) || k == 1
      if key == 28 || key == 30 || key == double('p')     % left/up  arrow
        plot(self, max(1, self.currentImage-1));
        figure(fig);
      elseif key == 29 || key == 31 || key == double('n') || key == 32 % right/down arrow
        plot(self, min(numel(self.images), self.currentImage+1));
        figure(fig);
      elseif key == double('i') || key == double('s') % Ignore/Skip
        self.images(self.currentImage).type = 'skip';
      elseif key == double('l') % Light
        self.images(self.currentImage).type = 'light';
      elseif key == double('g') % Goto image
        selection = listdlg(self, 'select image to go to','single');
        if ~isempty(selection)
          plot(self, selection);
          figure(fig);
        end
      elseif key == 27 || key == double('x') || key == double('q') % ESC/quit
        close(fig)
        return
      end
    end

    switch but
    case 'alt'              % ctrl-click or right button
      title('Setting images to Light')
      x = xs; y = ys;
    case 'normal'           % left-button
      title('Setting images to Skip')
    case {'open','extend'}  % shift-click or middle-button
      title('Opening first image')
    end
    % drag rectangle and get area [x y width height]
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units 
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    p2 = max(point1,point2);
    
    % determine which images are in the rectangle
    selection = find(x >= p1(1) & x <= p2(1) ...
                  &  y >= p1(2) & y <= p2(2));
    selection = x(selection); % image indices
    
    % (de-)select and replot
    for index=selection(:)'
      if strcmp(but, 'normal')
        self.images(index).type='skip';
      elseif strcmp(but, 'alt')
        self.images(index).type='light';
      else
        plot(self, index);
        figure(fig);
        break
      end
    end
    figure(fig);
    cla;
    [h, x, y, xs, ys] = plot_sharpness(self.images, self.currentImage);
    
  end % while
end % select_on_sharpness
