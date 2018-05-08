function initpkg(rm)
%INITPKG  Initialize the SBD-iPALM package
%
%  initpkg()  initializes the SBD-iPALM package by adding the appropriate
%  files to the MATLAB path.
%
%  initpkg(RM)  removes the same files from the MATLAB path if RM == TRUE.
%
%  Directories to be managed can be removed by changing the DIRS variable
%  at the beginning of the funciton.
%

% Directories + subdirectories to add to path
dirs = {
    '.',        false;
    'utils',    true;
};

warning('OFF', 'MATLAB:addpath:DirNotFound');
warning('OFF', 'MATLAB:rmpath:DirNotFound');

if nargin < 1 || isempty(rm);  rm = false;  end

pkgpath = mfilename('fullpath');
d = strsplit(pkgpath, {'\','/'});
pkgpath = pkgpath(1:end-numel(d{end}));

for d = 1:size(dirs,1)
    pathop([pkgpath dirs{d,1}], rm, dirs{d,2});
end

warning('ON', 'MATLAB:addpath:DirNotFound');
warning('ON', 'MATLAB:rmpath:DirNotFound');

end

function pathop(path, optype, recursive)
    if recursive
        path = genpath(path);
    end
    
    if optype
        rmpath(path);
    else
        addpath(path);
    end
end