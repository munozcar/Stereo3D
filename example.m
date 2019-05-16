%% Load stereo images 
imgl = imread('./Test_images/left_rectified.bmp');
imgr = imread('./Test_images/right_rectified.bmp');

%imgl = im2double(rgb2gray(imgl));
%imgr = im2double(rgb2gray(imgr));

imgl = im2double(imgl);
imgr = im2double(imgr);

dims = size(imgl);
xdim = dims(2);
ydim = dims(1);

% imgl = padarray(imgl,[600,50],'both');
% imgr = padarray(imgr,[600,50],'both');
%  
% dims = size(imgl);
% xdim = dims(2);
% ydim = dims(1);

addpath ~weinman/courses/CSC262/toolbox

%% Rectify Images (skip if images already rectified)
% Load camera parameters
load Calibration/Calib_Results_stereo.mat
R = R;
T = -T;

e1 = T/norm(T);
e2 = [-T(2) T(1) 0]/(sqrt(T(1)^2+T(2)^2));
e3 = cross(e1,e2);
Rrect = [e1'; e2; e3];

Rl = Rrect;
Rr = R*Rrect;

imgl = rectify(imgl,Rl,xdim,ydim,KK_left, fc_left(1));
imgr = rectify(imgr,Rr,xdim,ydim,KK_right, fc_right(1));
%%

figure(1);
imshow(imgl)
title('Left Image');
pause(0.5);

figure(2);
imshow(imgr);
title('Right Image');
pause(0.5);

%% Find correspondences

[xx, yy] = feat_matcher(imgl, imgr);

figure(1);
hold on;
plot(xx(:, 1), yy(:, 1), 'rx', 'markersize', 10, 'linewidth',2); 
figure(2);
hold on;
plot(xx(:, 2), yy(:, 2), 'rx', 'markersize', 10, 'linewidth',2); 

%% Estimate fundamental matrix and epipoles
[F, el, er] = epipolar(xx,yy);
N = size(xx,1);

% Plot correspondances and epipolar lines in the right image
for i = 1:N
    p_L = [xx(i, 1) yy(i, 1) 1]; 
    p_R = [xx(i, 2) yy(i, 2) 1]; 
    u = F*p_L';
    ur = p_R*F;
    x_R = [1 xdim]; 
    
    y_R1 = -(u(1)*x_R(1) + u(3))/u(2); 
    y_R2 = -(u(1)*x_R(2) + u(3))/u(2); 
    
    yr_R1 = (ur(1)*x_R(1) + ur(3))/ur(2); 
    yr_R2 = (ur(1)*x_R(2) + ur(3))/ur(2); 
    
    figure(2); 
    hold on; 
    plot([x_R(1) x_R(2)], [y_R1 y_R2], 'linewidth',0.5, 'Color', [0,1,0]);
    
    figure(1); 
    hold on; 
    plot([x_R(1) x_R(2)], [yr_R1 yr_R2], 'linewidth',0.5, 'Color', [0,1,0]);
end  

%% Reconstruction
load Calibration/Calib_Results_stereo.mat

load ./Calibration/Calib_Results_left.mat
Rl = Rc_1;
Tl = Tc_1;
load ./Calibration/Calib_Results_right.mat
Rr = Rc_1;
Tr = Tc_1;

Pl = KK_left*cat(2, Rl, Tl);
Pr = KK_right*cat(2, Rr, Tr);

ptsl = [xx(:,1) yy(:,1)]';
ptsr = [xx(:,2) yy(:,2)]';


set3D = reconstructor(Pl, Pr, ptsl, ptsr);

%% Plot figure
X = set3D(1, :)';
Y = set3D(2, :)';
Z = set3D(3, :)';
Tri = delaunay(X, Y);

figure;
trisurf(Tri, X, Y, Z, 'EdgeColor', 'none', 'FaceColor', 'red');
lighting phong;
view([90 90]);
camlight headlight;

figure;
scatter3(X,Y,Z);
view([90 90]);

