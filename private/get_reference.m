function [nbtostack, self] = get_reference(self, index)
% get_reference: return a 'reference' image (when not set)
%   when working in 'extend' mode, the control points are compared and added by
%   performing a 'fast' stack with diff.
 
  if nargin < 2, index=[]; end
   
  nbtostack = 0; ref = self.reference;
  
  % count how many LIGHT images we have
  light       =  strcmp('light', { self.images.type }) | ...
                 strcmp('reference', { self.images.type }) | ...
                 cellfun(@isempty, { self.images.type });
  if ~any(light)
    disp([ mfilename ': there is no LIGHT image, only dark and flat images. Aborting.' ]);
    self.refrence = [];
    return
  else
    nbtostack = sum(light);
  end
  
  % when an index is given, set it from the index
  if ~isempty(index)
    if 1 <= index && index <= numel(light) && light(index)
      ref = self.images(index);
    else
      disp([ mfilename ': Invalid Reference index (not a Light image, or out of bounds).' ])
      disp([ '    Ignoring (reference unchanged).' ]);
      return
    end
  end

  % is the reference unset ?
  if isempty(ref) || ~isstruct(ref)

    % we get the first 'light' image which has the highest sharpness
    sharpness    = [ self.images.sharpness ];
    % only retain the 'light' images for the reference
    sharpness(light == 0) = -inf;
    [~,sharpest] = max(sharpness);
    ref = self.images(sharpest); % copy the reference image structure

    ref.type = 'reference';
  end
  
  disp([ mfilename ':   using reference as ' ref.id ' [' num2str(ref.index) ']' ])

  % was the reference changed ? restore the initial 'type' to 'light'
  if isstruct(self.reference) && isstruct(ref)
    if self.reference.index ~= ref.index
      disp([ mfilename ': image ' self.reference.id ' set as LIGHT (was REFERENCE)' ]);
    end
    if 1 <= self.reference.index && self.reference.index <= numel(light) && light(self.reference.index)
      self.images(self.reference.index).type = 'light';
    end
  end
  
  ref.type                     = 'reference';
  self.images(ref.index).type  = 'reference';
  self.reference               = ref;
  
  % compute initial control points for the reference, if needed
  if isempty(ref.points) || ~isstruct(ref.points) || ...
    ~isfield(ref.points, 'x') || numel(ref.points.x) < 2
    % compute control points
    disp([ mfilename ':   computing control points for ' ref.id ' [' num2str(ref.index) '] REF' ])
    [self.reference] = cpselect(self, self.reference);
  end
  

