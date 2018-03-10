classdef mastrostack < handle

  properties
    nbControlPoints       = 20;
    toleranceTranslation  = 0.01; % in percent
    toleranceRotation     = 3;    % in deg
    thresholdDark         = 2;    % over 256: lower  is Dark (mostly dark)
    thresholdFlat         = 100;  % over 256: higher is Flat (mostly bright)
    reference             = 0;    % by default the first 'light' image
  %end % properties
  
  %properties(Access=private)
    catalogs = [];  % star catalogs for field identification (RA/DEC)
                    % obtained from SkyChart
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
    dark    = 0;  % stored as an uint16
    flat    = 0;  % stored as a double
    light   = 0;  % stored as a uint16
    lightN  = 0;  % stored as a uint16
  end % properties
  
  methods
  
    % I/O stuff ================================================================
    function ar = mastrostack

      % load catalogs: objects, stars
      disp([ mfilename ': Welcome ! Loading Catalogs:' ]);
      self.catalogs = load(mfilename);
      
      % display available catalogs
      for f=fieldnames(self.catalogs)'
        name = f{1};
        if ~isempty(self.catalogs.(name))
          num  = numel(self.catalogs.(name).RA);
          if isfield(self.catalogs.(name), 'Description')
            desc = self.catalogs.(name).Description;
          else desc = ''; end
          disp([ mfilename ': ' name ' with ' num2str(num) ' entries.' ]);
          disp([ '  ' desc ])
        end
      end
      
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
      if ~iscell(src), src = { src }; end
      
      img   = []; im = {};
      
      for index=1:numel(found)
      
        % read new image, or retrieve it
        [this_im, this_img] = imread_single(src{index}, self.images, flag);
        
        if isempty(this_img) continue; end  % invalid image
        
        % store if this is a new image
        if isnan(found(index))
          this_img.index = numel(self.images)+1;
          if isempty(self.images)
            self.images    = this_img;
          else 
            self.images(this_img.index) = this_img; 
          end
        end

        % add to return arguments
        if isempty(img), img        = this_img;
        else             img(end+1) = this_img; end
        im{end+1} = this_im;
      end % for
      
      if iscell(im) && numel(im) == 1
        im = im{1};
      end
      
    end % imread
    
    function img = load(self, source)
      % load: read an image, and store its information
      %
      %  img = load(self, filename or matrix or index)
      %    returns the image matrix and information
      
      % load images, but do not get matrix from files
      [~, img] = imread(self, source, 0);
      
    end % load
    
    function [found, src] = exist(self, source)
      % exist: check if an image is already loaded
      %
      %   found = exist(self, source)
      %     returns the index of any matching image, or nan (not found)
      
      % handle array of images
      found = []; src = {};
      if iscell(source) || (isnumeric(source) && ~isscalar(source) ...
        && numel(source) == max(size(source)))
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
          end
        end % for
      end
      
    end % exist
    
    function label(self, img, lab)
      % label: sort images into categories
      
      % we can sort intensity or sharpness
      % and show a listdlg to select, or select auto on threshold
      
    end % label
    
        
    function img = get.dark(self)
      % dark: compute the mean Dark frame, using labelled images
      sumDarks = 0;
      nbDarks  = 0;
      for index=1:mumel(self.images)
        this = self.images(index);
        if ~strcmp(img.type, 'dark')
          continue;
        end
        [img, im] = imread(self, index);  % read image (re-read if needed)
      end
    end
    
    % compute stuff ============================================================

    function im = flat(self, im)
      % flat: correct an image for vigneting and background
      if nargin < 2, return; end
      
      % correct for read-out/sensor noise (subtract)
      if ~isempty(self.dark)
        % must cast the dark to the image.
        if strcmp('uint8', class(im))
          dark = self.dark/255;
        else
          dark = self.dark;
        end
        im = im - dark;
      end
      
      % correct for flat field to compensate vigneting (divide)
      if ~isempty(self.flat)
        for l=1:size(im,3)
          layer = double(im(:,:,l));
          layer = layer./flat;
          im(:,:,l) = cast(layer, class(im));
        end
      end
      
    end
    
    
    function [points, img] = cpselect(self, img, im)
      % cpselect: Control Point automatic search
      %
      %   cpselect(self):     find control points for all images
      %   cpselect(self, img, im): find control points for single image, where
      %     arguments are those returned from: [im, img]=imread(self, 'image');    
      %
      %   Up to self.nbControlPoints are selected, ensuring that all of these are
      %   separated with at least self.toleranceTranslation*10 in % of image size.
      
      if nargin < 3
        % get the images (only the meta info). We shall read images 1-by-1
        [im, img] = imread(self, img, 0); % do not re-read image
      end
      points = [];
      
      if numel(img) == 1
        % search control points for a single image
        if isempty(img.points)
          [im, ~]   = imread(self, img);
          im        = rgb2gray(im);
          
          
          N         = self.nbControlPoints;
          tol_trans = self.toleranceTranslation;
          
          % search points=(x,y,m,sx,sy) and store
          img.points             = find_control_points(im, N, tol_trans*10);
          if numel(self.image) <= 1, self.images = img;
          else self.images(img.index) = img; end
        end
        points = img.points;
      else
        % search control points for an array of images
        for index=1:numel(img)
          this_points = cpselect(self, img(index));
          if isempty(points) 
            points = this_points;
          else
            points(end+1) = this_points;
          end
        end
      end
      
    end % cpselect
    
    function [ret_t, theta] = diff(self, img1, img2)
      % diff: compute the difference between img2 and img1 used as reference.
      %
      %   [ret_t, theta] = diff(self, img1, img2)
      %     returns the translation and rotation between img1 and img2.
      %   [ret_t, theta] = diff(self, img2)
      %     use the reference 'light' image as img1.
      
      ret_t=[]; theta=[];
      if nargin < 3
        if self.reference && self.reference < numel(self.images)
          img2 = img1;
          img1 = self.images(self.reference);
        else return; end
      elseif nargin < 2
        return
      end
      % get image and reference
      [~, img1] = imread(self, img1, 0);
      [~, img2] = imread(self, img2, 0);
      
      % here we need img1 and img2. Make sure we have control points
      cpselect(self, img1);
      cpselect(self, img2);
      
      % the images should have the same scaling
      tol_trans = self.toleranceTranslation;
      tol_rot   = self.toleranceRotation;
      [x1,y1,x2,y2, p1_orig,p2_orig, p1_axis,p2_axis] = ...
        analyse_dist_angles(img1.points, img2.points, tol_trans, tol_rot);
      
      % get the best similarity
      if ~p1_orig || ~p2_orig || ~p1_axis || ~p2_axis
        disp([ mfilename ': WARNING: not enough control points for ' img1.id ' ' img2.id ' (axis).' ])
        return
      end
  
      % identify the current points which are common to the reference, and match
      % distances AND rotations
      p1 = p1_orig; p2=p2_orig;
      d1  = sqrt( (x1-x1(p1)).^2 + (y1-y1(p1)).^2 );
      t1  = atan2( y1-y1(p1),       x1-x1(p1))*180/pi;
      d2  = sqrt( (x2-x2(p2)).^2 + (y2-y2(p2)).^2 );
      t2  = atan2( y2-y2(p2),       x2-x2(p2))*180/pi;
      [ok1,ok2] = find_similar2(d1,             d2,             tol_trans*img2.image_size(1), ...
                                t1-t1(p1_axis), t2-t2(p2_axis), tol_rot*pi/180);
                                
      if numel(ok1) <= 1 || numel(ok2) <= 1
        disp([ mfilename ': WARNING: not enough control points for ' img1.id ' ' img2.id ' (common).' ])
        return
      end
  
      % we make a check for wrong guesses
      delta = (t1(ok1)-t1(p1_axis)) - (t2(ok2)-t2(p2_axis));
      bad   = find(abs(delta) > tol_rot);
      ok1(bad) = [];
      ok2(bad) = [];
      if numel(ok1) <= 1 || numel(ok2) <= 1
        disp([ mfilename ': WARNING: not enough control points for ' img1.id ' ' img2.id ' (common2).' ])
        return
      end
  
      % compute the mean translation (x,y) and rotation angle wrt to reference image
      x1  = x1(ok1); y1=y1(ok1);
      x2  = x2(ok2); y2=y2(ok2);
      
      % theta2-theta1 is an estimate of the rotation angle
      % compute the affine transformatin
      [ret_R, ret_t] = rigid_transform_3D([x2 ; y2]', [x1 ; y1]');
      if isempty(ret_R)
        disp([ mfilename ': WARNING: invalid affine transformation. Skipping.'])
        return
      end
      % compute an estimate of the translation and rotation from the identified
      % control points orig and axis
      theta1 = t1(p1_axis); % t1(p1_orig) == 0
      theta2 = t2(p2_axis);
      theta = asind(ret_R(1,2));
      if abs(theta-(theta2-theta1)) > tol_rot/3
        disp([ mfilename ': WARNING: invalid affine rotation. Skipping.']);
        return
      end
  
      img2.angle        = theta;
      img2.translation  = ret_t;
      img2.rotation     = ret_R;
      img2.reference    = img1.index;
      self.images(img2.index) = img2;
    end % diff
    
    function stack(self, img)
      % stack: stack all images on the reference
      
      if nargin == 1
        img = 1:numel(self.images);
      end
      
      % get the images (only the meta info). We shall read images 1-by-1
      [~, img] = imread(self, img, 0); % do not re-read image
      
      [im,M] = imaffine(im, ret_R, ret_t);
      self.lightN=M0+M;
      clear M
      self.light = self.light + uint16(im); % add rotated/translated image on the reference
      
      for index=1:numel(img)
        this_img = img(index);
        if strcmp('light', this_img.type)
          % stack light images
          [im, ~] = imread(self, img); % we read the image
          [im,M] = imaffine(im, ret_R, ret_t);
          self.lightN=M0+M;
          clear M
          self.light = self.light + uint16(im); % add rotated/translated image on the reference
 
        end
      end
    end % stack
    
    % plotting stuff ===========================================================
    
    function h = plot(self, img)
      % plot: display loaded images
      %
      %   plot(self, img): display specified image. img can be filename, RGB image
      %     image structure or index.
      
      if nargin == 1
        img = 1:numel(self.images);
      end
      
      % get the images (only the meta info). We shall read images 1-by-1
      [~, img] = imread(self, img, 0); % do not re-read image

      % plot no more than ~20 images
      if numel(img) > 20
        ratio = ceil(numel(img)/20);
        img   = img(1:ratio:end);
      end
      
      % now build the grid
      n = ceil(sqrt(numel(img)));
      m = ceil(numel(img)/n);
      h = [];
      for index=1:numel(img)
        % use a grid to display images
        subplot(m,n,index);
        [im, ~] = imread(self, img); % we read the image and plot it
        if isempty(im), im=img.thumbnail; end
        h = [ h image(im) ];
        set(gca,'XTick',[], 'YTick',[]); title(img.id);

        % we add the control points, scaled to the thumbnail/image
        if ~isempty(img.points)
          hold on
          x = img.points.x*size(im,1)/img.image_size(1);
          y = img.points.y*size(im,2)/img.image_size(2);
          scatter(y,x,'o','g');
          hold off
        end
      end
      
    end % plot
    
    function info = about(self, img)
      % about: display image information
      info = {};
      if nargin == 1
        info{end+1} = 'Ma(e)stroStack for Matlab';
        info{end+1} = 'A Software to load, identify, align and stack astrophotography images.';
        info{end+1} = '(c) E. Farhi <https://github.com/farhi/matlab-mastrostack>';
        helpdlg(info, 'Ma(e)stroStack: About');
        return
      end
      [~, img] = imread(self, img, 0); % do not re-read image
      % print information
      f = {'FileName','Software', 'DateTime', 'WhiteBalance', 'ExposureProgram', ...
        'LightSource', 'ExposureMode', 'ExposureBiasValue', 'FNumber', ...
        'ExposureTime', 'DigitalCamera', 'XResolution', 'YResolution', ...
        'ISOSpeedRatings','Author', 'Comment', ...
        'Make','Model','Source'}; 
      f = sort(f);
      
      if numel(img) == 1
        % get single EXIF data
        exif = img.exif;
        for i=1:numel(f)
          if isfield(exif, 'DigitalCamera') && isfield(exif.DigitalCamera, f{i})
            exif.(f{i}) = exif.DigitalCamera.(f{i});
          end
          if isfield(exif, f{i})
            try; info{end+1} = [ f{i} ' = ' num2str(exif.(f{i})) ]; end
          end
        end
      else
        % an array of image EXIF data      
        for index=1:numel(img)
          info{end+1} = about(self, img(index));
        end
      end
      % display dialog when not getting output
      if nargout == 0
        helpdlg(info, img.id);
      end
    end % about
  
  end % methods
  
end % mastrostack

% ------------------------------------------------------------------------------













  
  


