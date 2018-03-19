function [ob] = obops()
ob.inner = @inner;
ob.proj = @proj;
ob.proj2tan = @proj2tan;
ob.e2rgrad = @proj2tan;
ob.dist = @dist;
ob.Exp = @Exp;
ob.Log = @Log;
end

function [out] = inner(A,B)
out = dot(A(:), B(:));
end

function [out] = proj(in)
out = zeros(size(in));
for i = 1:size(in,3)
    tmp = in(:,:,i);
    if norm(tmp, 'fro') > 0
        out(:,:,i) = tmp/norm(tmp, 'fro');
    else
        out(:,:,i) = tmp;
    end
end
end

function [out] = proj2tan(A, in)
out = zeros(size(A));
for i = 1:size(in,3)
    out(:,:,i) = in - inner(A(:,:,i), in(:,:,i))*A(:,:,i);
end
end

function [out] = dist(A, B)
out = 0;
for i = 1:size(A,3)
   out = out + real(acos(inner(A(:,:,i), B(:,:,i))))^2;
end
out = sqrt(out);
end

function [out] = Exp(A, del)
out = zeros(size(A));
for i = 1:size(A,3)
    ndel = norm(del(:,:,i),'fro');
    if ndel > 0
        out(:,:,i) = cos(ndel)*A(:,:,i) +  sin(ndel)*del(:,:,i)/ndel;
    else
        out(:,:,i) = A(:,:,i);
    end
end
out = proj(out);
end

function [out] = Log(A, B)
out = proj(proj2tan(A, B-A));
for i = 1:size(A,3)
    out(:,:,i) = dist(A(:,:,i), B(:,:,i)) * out(:,:,i);
end
end