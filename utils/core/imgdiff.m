classdef imgdiff < handle
properties
    stepsz;
end

properties (Access = private)
    imsz;
    dirs = 1:2;
    splits;
end

methods
function o = imgdiff(imsz, dirs)
    o.imsz = imsz;
    if nargin >= 2 && ~isempty(dirs);  o.dirs = dirs;  end
end

function v = funh(o, x)
    ndirs = numel(o.dirs);
    imsz = o.imsz;
    
    function [v, stepdec] = diff1(x)
        v = diff(x,1,1);  v = v(:);
        stepdec = 2;
    end
    
    function [v, stepdec] = diff2(x)        
        v = diff(x,1,2);  v = v(:);
        stepdec = 2;
    end
    
    % DIFF/DIV 3/4 -- IN PROGRESS
    function [v, stepdec] = diff3(x)
        v = conv([0; -1; zeros(imsz(1)-2,1); 1], x(:), 'full');
        v = v(1:(imsz(1)-1)*(imsz(2)-1))/sqrt(2);
        stepdec = sqrt(2);
    end
    
    function [v, stepdec] = diff4(x)
        [v, stepdec] = diff3(flipud(x));
    end
    
    dfuns = {@diff1 @diff2 @diff3 @diff4};
    v = [];  o.stepsz = 0;
    o.splits = zeros(ndirs+1,1);
    for i = 1:ndirs
        [tmp, stepdec] = dfuns{o.dirs(i)}(x);
        v = [v ; tmp];  %#ok<AGROW>
        o.splits(i+1) = numel(tmp);
        o.stepsz = o.stepsz + stepdec;
    end
    o.splits = cumsum(o.splits);
    o.stepsz = 1/o.stepsz;
end

function x = adjh(o, v)
    imsz = o.imsz;
    
    function x = div1(v)
        v = reshape(v, [(imsz(1)-1) imsz(2)]);
        x = [-v(1,:); -diff(v,1,1); v(end,:)];
    end
    
    function x = div2(v)
        v = reshape(v, [imsz(1) (imsz(2)-1)]);
        x = [-v(:,1) -diff(v,1,2) v(:,end)];
    end
    
    function x = div3(v)
        vpad = [zeros(imsz(1)-1,1); v];
        x = conv([1; zeros(imsz(1)-2,1); -1], vpad, 'full');
        x = [x; zeros(prod(imsz),1)];
        x = x(imsz(1)+(0:prod(imsz)-1));
        x = reshape(x, imsz);
    end
    
    function x = div4(v)
        x = flipud(div3(v));
    end
    
    dfuns = {@div1 @div2 @div3 @div4};
    x = zeros(o.imsz);
    for i = 1:numel(o.dirs)
        tmp = (o.splits(i)+1):o.splits(i+1);
        x = x + dfuns{o.dirs(i)}(v(tmp));
    end
end

end
end