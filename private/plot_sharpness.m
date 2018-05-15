function [h, x_light, y_light, x_skip, y_skip]=plot_sharpness(images, current)
  
  if nargin < 2, current = 0; end
  y = [ images.sharpness ] .* [ images.intensity ];
  y(~isfinite(y)) = nan;
  % we create data sets for Dark, Light, Flat and Skipped images
  x_dark=[]; x_light=[]; x_flat=[]; x_skip=[];
  y_dark=[]; y_light=[]; y_flat=[]; y_skip=[];
  for index=1:numel(images)
    typ = lower(images(index).type);
    if isempty(typ), typ = 'light'; end
    switch typ
    case 'dark'
      x_dark(end+1)=index;
      y_dark(end+1)=y(index);
    case 'flat'
      x_flat(end+1)=index;
      y_flat(end+1)=y(index);
    case 'skip'
      x_skip(end+1)=index;
      y_skip(end+1)=y(index);
    otherwise
      x_light(end+1)=index;
      y_light(end+1)=y(index);
    end
  end
  % plot
  h1=plot(x_light,y_light, 'bo');
  hold on
  h2=plot(x_dark, y_dark,  'ks');
  h3=plot(x_flat, y_flat,  'mv');
  h4=plot(x_skip, y_skip,  'gx');
  if current >=1 && current <= numel(images)
    h5=plot(current, y(current), 'c+'); set(h5, 'MarkerSize',10);
  else h5=[];
  end
  set(gca, 'YScale','log');
  % set-up the legend
  leg = {};
  if ~isempty(h1), leg{end+1} = 'Light'; end
  if ~isempty(h2), leg{end+1} = 'Dark'; end
  if ~isempty(h3), leg{end+1} = 'Flat'; end
  if ~isempty(h4), leg{end+1} = 'Skip'; end
  if ~isempty(h5), leg{end+1} = 'Current'; end
  hold off
  xlabel('Image index')
  ylabel('Sharpness (higher is better)');
  title('Sharpness: {\color{blue}Left-clik}:SKIP {\color{blue}right-click}:LIGHT {\color{blue}shift-click} OPEN');
  h = [ h1, h2, h3, h4, h5 ];
  legend(h, leg);
  xlim([ 0 numel(images)+1 ])
  
