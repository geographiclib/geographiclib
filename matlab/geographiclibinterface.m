function geographiclibinterface;
% geographiclibinterface.m
% Use mex to compile interface to GeographicLib
%
% This has been tested with
%
%   Octave 3.2.2 and g++ 4.4.4 under Linux
%   Matlab 2007a and Visual Studio 2005 under Windows
%   Matlab 2008a and Visual Studio 2005 under Windows
%   Matlab 2008a and Visual Studio 2008 under Windows
%
% Note that the geoidheight causes Matlab to CRASH (in the last
% configuration above).  The crash happens on the second call.
%
% Copyright (c) Charles Karney (2010) <charles@karney.com> and licensed under
% the LGPL.  For more information, see http://geographiclib.sourceforge.net/
%
% $Id$
  funs = {'geodesicdirect', 'geodesicinverse', 'geodesicline', ...
	  'geoidheight', ...
	  'utmupsforward', 'utmupsreverse', ...
	  'mgrsforward', 'mgrsreverse'};
  if ispc,
    warning(['geoidheight (compiled with Visual Studio 2008) ',...
	     'causes Matlab (2008a) to CRASH']);
    incdir='../include';
    libdir='../windows/Release';
    lib='GeographicLib';
  else
    incdir='/usr/local/include';
    libdir='/usr/local/lib';
    lib='Geographic';
  end
  for i=1:size(funs,2),
    fprintf('Compiling %s...', funs{i});
    if ispc,
      mex( ['-I' incdir], ['-L' libdir], ['-l' lib], [funs{i} '.cpp'] );
    else
      mex( ['-I' incdir], ['-L' libdir], ['-l' lib], ['-Wl,-rpath=' libdir], ...
        [funs{i} '.cpp'] );
    end
    fprintf(' done.\n');
  end
end
