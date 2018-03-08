function [Achk] = rev(A, m)
if nargin < 2
    m = size(A);
end
Achk = zeros(m);
Achk(1:size(A,1), 1:size(A,2)) = A;
Achk = circshift(circshift(rot90(Achk,2), 1), 1, 2);
end