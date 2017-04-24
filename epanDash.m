function [ out ] = epanDash( x,y )
%EPANDASH Summary of this function goes here
%   Depricated
% point = [x, y];
% disp(point);
    normsq = x*x + y*y;
    if normsq < 1
        out = 2 / 3.14;
    else
        out = 0;
    end

end

