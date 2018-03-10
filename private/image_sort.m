function sort()

  % auto sort into Light, Dark and Flat
    if isempty(type)
      if img.image_sum/prod(img.image_size) < self.thresholdDark
        img.type        = 'dark';
      elseif img.image_sum/prod(img.image_size) > self.thresholdFlat
        img.type        = 'flat';
      else
        img.type        = 'light';
      end
    else
      img.type = type;
    end
    % determine control points if needed
    if strcmp(img.type, 'light')
      if ~self.reference
        self.reference = img.index;
      end
      [~, img] = cpselect(self, img, im);  % find control points
    end
