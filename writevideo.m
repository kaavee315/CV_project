function writevideo (ff, img)
%depricated


% open(ff);
% disp(size(img));
% figure, imshow(img);
writeVideo(ff, img);
% close(ff);

% [H,W,NF] = size(vid);
% ff = VideoWriter(filename);
% ff.FrameRate = framerate;
% 
% open (ff);
% for i=1:NF
%     writeVideo(ff,vid(:,:,i));
% end
% close (ff);