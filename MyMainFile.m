%% MyMainScript
tic;
clear;
% set up file path
folder = 'data';
% vid_file_path = strcat(folder,'/small.mp4');
vid_file_path = '../cv_proj_data/rahul.avi';

% find num frames
% inp = VideoReader(vid_file_path);
% num_frames = inp.NumberOfFrames;
% num_frames = 50;

%setup to read video data
inp = VideoReader(vid_file_path);
% inp.CurrentTime = 28;
vidHeight = inp.Height;
vidWidth = inp.Width;
% s = struct('cdata',zeros(inp_vid_obj.Height, inp_vid_obj.Width, 3));

%hyper parameters
binSize = 8;
numBins = 256/binSize;
h = 180;

%output video
outFileName = 'data/outvid_karan.avi';
ff = VideoWriter(outFileName);
ff.FrameRate = inp.FrameRate;
open(ff);

%read first frame and mark target
frame_prev = rgb2gray(readFrame(inp));
% rect = [140 160 60 150];
y0 = [710 800];
rect = [y0(1) - h, y0(2) - h, 2*h, 2*h];
frame_prev_out = rgb2gray(insertShape(frame_prev, 'Rectangle', rect, 'LineWidth', 8));
writeVideo(ff, frame_prev_out);
imshow(frame_prev_out);

%find qu of prev frame
qu = pu(frame_prev, y0, h, binSize, numBins);
% qu = zeros(numBins,1);
% for i = rect(1):rect(1) + rect(3)
%     for j = rect(2):rect(2) + rect(4)
%         qu(b(frame_prev(i,j), binSize),1) = qu(b(frame_prev(i,j), binSize),1) + epan((i-rect(1))/rect(3),(j-rect(2))/rect(4));
%     end
% end
% qu = qu / sum(qu);
% y0 = [uint8(rect(1)+(rect(3)/2)), uint8(rect(2)+(rect(4)/2))];

% for i = 2:100
i = 1;
while 1
    if hasFrame(inp)
        frame_curr = rgb2gray(readFrame(inp));
        while 1
%             disp(y0);
            puy = pu(frame_curr, y0, h, binSize, numBins);

            BC = sqrt(puy).' * sqrt(qu);
            weights_pool = sqrt(qu ./ puy);
            weights_pool(~isfinite(weights_pool)) = 0;
            
            %calculating y1
            y1 = [0 0];
            denominator = 0; 
            
            xmin = max(1,uint32(y0(1)-h));
            xmax = min(size(frame_curr,2),uint32(y0(1)+h));
            ymin = max(1,uint32(y0(2)-h));
            ymax = min(size(frame_curr,1),uint32(y0(2)+h));
            
            xrange = double(xmin:xmax);
            yrange = double(ymin:ymax);
            xrangeNorm = (xrange-y0(1))/h;
            yrangeNorm = (yrange-y0(2))/h;
            xrangeNorm = repmat(xrangeNorm,ymax-ymin+1,1);
            yrangeNorm = repmat(yrangeNorm.',1,xmax-xmin+1);
            epanDashKernel = (xrangeNorm.^2 + yrangeNorm.^2);
            epanDashKernel(epanDashKernel < 1) = 2/3.14;
            epanDashKernel(epanDashKernel >= 1) = 0;
            
            weights = double(frame_curr(ymin:ymax, xmin:xmax));
            weights = floor(weights/binSize) + 1;
            
            for j = 1:numBins
                weights(weights == j) = weights_pool(j);
            end
            
            xrange = repmat(xrange,ymax-ymin+1,1);
            yrange = repmat(yrange.',1,xmax-xmin+1);
            
            denominatorMatrix = weights .* epanDashKernel;
            y1_1 = xrange .* denominatorMatrix;
            y1_2 = yrange .* denominatorMatrix;
            
            y1(1) = sum(y1_1(:));
            y1(2) = sum(y1_2(:));
            
            y1 = y1 / sum(denominatorMatrix(:));
            
%             for j = uint64(y0(1) - h):uint64(y0(1) + h)
%                 for k = uint64(y0(2) - h):uint64(y0(2)+ h)
%                     if j > 0 && j < size(frame_curr, 2) && k > 0 && k < size(frame_curr, 1)
%                         w = weights_pool(b(frame_curr(k,j), binSize));
% %                         disp([j,k]);
%                         y1(1) = y1(1) + j * w * epanDash((y0(1)-j)/h, (y0(2)-k)/h);
%                         y1(2) = y1(2) + k * w * epanDash((y0(1)-j)/h, (y0(2)-k)/h);
% %                         disp(num2str(epanDash((y0(1)-j)/h, (y0(2)-k)/h)));
%                         denominator = denominator + w * epanDash((y0(1)-j)/h, (y0(2)-k)/h);
%                     end
%                 end
%             end
%             y1 = y1 / denominator;

            puy = pu(frame_curr, y1, h, binSize, numBins);
            BCnew = sqrt(puy).' * sqrt(qu);
%             disp(['this is BC ', num2str(BC)]);
            while(BCnew < BC)
%                 disp(['hello ', num2str(BCnew)]);
                y1 = (y0 + y1)/2;
                puy = pu(frame_curr, y1, h, binSize, numBins);
                BCnew = sqrt(puy).' * sqrt(qu);
%                 if norm(y1 - y0) < 1
%                     break;
%                 end
            end

            if norm(y1 - y0) < 1
                break;
            else
                y0 = y1; 
            end
        end
        y0 = y1;
%         disp(y0);
        qu = pu(frame_curr, y0, h, binSize, numBins);
        
        %mark frame and add to video file
        rect = [y0(1) - h, y0(2) - h, 2*h, 2*h];
        frame_curr_out = rgb2gray(insertShape(frame_curr, 'Rectangle', rect, 'LineWidth', 8));
        writeVideo(ff, frame_curr_out);
        disp(['writing frame ', num2str(i)]);
        i = i+1;
    else
        break;
    end
    
end

close(ff);
%read video
% k = 1;
% while hasFrame(inp_vid_obj)
%     s(k).cdata = rgb2gray(readFrame(inp_vid_obj));
%     if k == 50
%         break
%     end
%     k = k + 1;
% end

%setup output video data structure
% out_vid = zeros(vidHeight, vidWidth, 50);

% insert circle
% for i = 1:num_frames
%     out_vid(:,:,i) = rgb2gray(insertShape(s(i).cdata, 'circle' , [50 50 30], 'LineWidth', 5));
% end

% writevideo('data/newvid.avi', out_vid/max(out_vid(:)), inp.FrameRate);

toc;