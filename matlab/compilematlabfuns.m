function compilematlabfuns;
% compilematlabfuns.m
% Use mex to compile interface to GeographicLib
%
% Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
% the LGPL.  For more information, see http://geographiclib.sourceforge.net/
%
% $Id$
  funs = {'geoidheight', ...
    'utmupsforward', 'utmupsreverse', ...
    'mgrsforward', 'mgrsreverse'};
  if ispc,
    incdir='../include';
    libdir='../windows/Release';
    lib='GeographicLib';
  else
    incdir='/usr/local/include';
    libdir='/usr/local/lib';
    lib='Geographic';
  end
  for i=1,size(funs,2),
    fprintf('Compiling %s...', funs{i});
    mex( ['-I' incdir], ['-L' libdir], ['-l' lib], [funs{i} '.cpp'] );
    fprintf(' done.\n');
  end
end
