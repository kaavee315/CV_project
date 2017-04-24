function [out] = epan(x,y)
%Depricated
normsq = x*x + y*y;
if normsq < 1
    out = 2*(1 - x*x - y*y) / 3.14;
else
    out = 0;
end

