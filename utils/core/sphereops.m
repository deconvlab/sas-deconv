function [s] = sphereops()
s.inner = @inner;
s.proj = @proj;
s.proj2tan = @proj2tan;
s.e2rgrad = @proj2tan;
s.dist = @dist;
s.Exp = @Exp;
s.Log = @Log;
end

function [out] = inner(A,B)
out = dot(A(:), B(:));
end

function [out] = proj(in)
tmp = norm(in(:));
if tmp > 0
    out = in/tmp;
else
    out = in;
end
end

function [out] = proj2tan(A, in)
out = in - inner(A, in)*A;
end

function [out] = dist(A, B)
out = real(acos(inner(A, B)));
end

function [out] = Exp(A, del)
ndel = norm(del(:));
if ndel > 0
    out = cos(ndel)*A +  sin(ndel)*del/ndel;
else
    out = A;
end
out = proj(out);
end

function [out] = Log(A, B)
out = dist(A, B) * proj(proj2tan(A, B-A));
end