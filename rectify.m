function [rectimg] = rectify(img, R, xdim, ydim, KK, fc)

% FUNCTION: rectify, rectifies an image given camera parameters
% PARAMETERS:
%               - img: an image 
%               - R: rotation matrix (3x3)
%               - xdim, ydim: dimensions of the image
%               - KK: camera matrix with intrinsic parameters
%               - fc = focal length
% PRODUCES: 
%               - rectimg: the rectified image
%               

% Create new image plane
rectimg = zeros(ydim,xdim);
% Generate rectification matrix
M = KK*inv(R)*inv(KK);

for i=1:ydim
    for j=1:xdim
        point = M*[i j 1]';
        point = point.*point(3);
        point(point<1)=1;
        
        % Correct cordinates if out of bounds
        if point(1)>ydim
            point(1) = ydim-1;
        end
         if point(2)>xdim
            point(2) = xdim-1;
        end
            
        right_down = ceil([point(1) point(2)]);
        left_up = floor([point(1) point(2)]);
        left_down = [ceil(point(1)) floor(point(2))];
        right_up = [floor(point(1)) ceil(point(2))];
        
        
        % Interpolate brightness from four nearest neighbors
        neighbors = [img(right_down(1), right_down(2)) img(left_up(1), left_up(2))...
            img(left_down(1), left_down(2)) img(right_up(1), right_up(2))];
        
        interp = mean(neighbors);
        
        rectimg(i,j) = interp;
    end
end