function geoid = geoid_load_file(filename)
  geoid.file = filename;
  fid = fopen(geoid.file, 'r');
  if fid == -1, error(['Cannot open ' geoid.file]),  end
  done = 0;
  for i = 1:16
    if done == 3, break, end
    text = fgetl(fid);
    if ~ischar(text), continue, end
    text = strsplit(text, ' ');
    if length(text) == 3
      if strcmp(text{2}, 'Offset')
        geoid.offset = str2double(text{3});
        done = bitor(done, 1);
      elseif strcmp(text{2}, 'Scale')
        geoid.scale = str2double(text{3});
        done = bitor(done, 2);
      end
    elseif length(text) == 2
      break
    end
  end
  fclose(fid);
  if done ~=3, error('Cannot find scale and offset'), end
  geoid.im = imread(geoid.file);
  if ~isa(geoid.im, 'uint16'), error('Wrong image type'), end
  geoid.h = size(geoid.im, 1);
  geoid.w = size(geoid.im, 2);
  if ~(bitand(geoid.w, 1) == 0 && bitand(geoid.h, 1) == 1)
    error('Image width/height not even/odd')
  end
end
