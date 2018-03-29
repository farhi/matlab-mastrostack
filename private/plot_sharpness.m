function [h, x_light, y_light]=plot_sharpness(images)
  y = [ images.sharpness ]./[ images.width ];
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
  h=plot(x_light,y_light, 'bo', ...
         x_dark, y_dark,  'ks', ...
         x_flat, y_flat,  'mv', ...
         x_skip, y_skip,  'gx');
  xlabel('Image index')
  ylabel('Sharpness (higher is better)');
  title('Sharpness');
  legend(h, 'Light','Dark','Flat','Skip');
  xlim([ 0 numel(images)+1 ])
