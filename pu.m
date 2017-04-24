function [out] = pu(frame, y, h, binSize, numBins)
%PU Summary of this function goes here
%   Detailed explanation goes here
%     disp(y);
%     disp([size(frame, 1), size(frame,2)]);
    centrex = y(1);
    centrey = y(2);
    
    xmin = max(1,uint32(centrex-h));
    xmax = min(size(frame,2),uint32(centrex+h));
    ymin = max(1,uint32(centrey-h));
    ymax = min(size(frame,1),uint32(centrey+h));
    
    out = zeros(numBins,1);
%     disp([num2str(ymin),' ',num2str(ymax),' ',num2str(xmin), ' ',num2str(xmax)]);
    window = double(frame(ymin:ymax, xmin:xmax));
    window = floor(window/binSize) + 1;
    
    xrange = double(xmin:xmax);
    yrange = double(ymin:ymax);
    xrange = (xrange-centrex)/h;
    yrange = (yrange-centrey)/h;
    xrange = repmat(xrange,ymax-ymin+1,1);
    yrange = repmat(yrange.',1,xmax-xmin+1);
    epanKernel = (1 - (xrange.^2 + yrange.^2)) * (2/3.14);
    epanKernel(epanKernel <= 0) = 0;
    
    for i = 1:numBins
        out(i,1) = sum(epanKernel(window == i));
    end
    
%     for i = uint64(centrex - h):uint64(centrex + h)
%         for j = uint64(centrey - h):uint64(centrey + h)
%             if i > 0 && i < size(frame, 2) && j > 0 && j < size(frame, 1)
%                 out(b(frame(j,i), binSize),1) = out(b(frame(j,i), binSize),1) + epan((i-centrex)/h,(j-centrey)/h);
%                 count = count + 1;
%             end
%         end
%     end
    if sum(out) == 0
        error('quitting due to error in pu..');
    end
    out = out / sum(out);
end

