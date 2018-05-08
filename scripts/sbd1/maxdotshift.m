function [maxdot] = maxdotshift(aref, aest)
%MAXDOTSHIFT  Max abs. scaled dot product between shifts of two kernels.
%
%  Tests AEST against shift-truncations of AREF.
%

p = numel(aref);  aref = aref(:)/norm(aref(:));
p_ = numel(aest); aest = aest(:)/norm(aest(:)); 
assert(p <= p_, 'Ref. kernel must not be longer than estimate.');

tmp = zeros(2*p+2,1);
tmp(1:p) = aref;
aref = tmp;

maxdot = 0;
for tau = -ceil(p/2):ceil(p/2)
    tmp = circshift(aref, tau);
    tmp = tmp(1:p);
    tmp = xcorr(tmp/norm(tmp), aest);
    maxdot = max(maxdot, max(abs(tmp)));
end
end