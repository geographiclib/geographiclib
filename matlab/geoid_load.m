function geoid = geoid_load(name, dir)
  if nargin < 1
    file = geoid_file;
  elseif nargin < 2
    file = geoid_file(name);
  else
    file = geoid_file(name, dir);
  end
  geoid = geoid_load_file(file);
end
