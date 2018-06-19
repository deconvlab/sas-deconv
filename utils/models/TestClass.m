classdef TestClass < handle
properties
    x = 5;
end

methods
function out = abc(o, niter)
    for it = 1:niter
        o = o.step();
    end
    out = o.x;
end

function o = step(o)
    o.x = o.x + 1;
end
end
end