% Epipolar Geometry
% Getting Correspondences

% Load stereo images 
img_left = im2double(rgb2gray(imread('/home/munozcar/Desktop/moto1.png')));
img_right = im2double(rgb2gray(imread('/home/munozcar/Desktop/moto2.png'))); 

figure(1);
imshow(img_left);
title('Left Image');
pause(0.5);

figure(2);
imshow(img_right);
title('Right Image');
pause(0.5);


features_left = kpdet(img_left);
features_right = kpdet(img_right);

descriptors_left = kpfeat(img_left, features_left);
descriptors_right = kpfeat(img_right, features_right);

[xx, yy] = kpmatch(features_left, features_right, descriptors_left, descriptors_right);

F = kpfund(xx,yy);

% Visualizing Epipoles
% Now that we have an estimate of the fundamental matrix, we can find the
% epipolar lines for each correspondence. 

img_width = size(img_left, 2); 

% Choose some point on left image
p_L = [100 50 1]; 
% Find right image epipolar line coefficient for this point
u = F*p_L';

x_R = [1 img_width]; 
y_R1 = -(u(1)*x_R(1) + u(3))/u(2); 
y_R2 = -(u(1)*x_R(2) + u(3))/u(2); 
figure(3);
imshow(img_left);
hold on; 
plot(p_L(1), p_L(2), 'rs', 'markersize', 10, 'linewidth',3); 
title("Left Image"); 
pause(0.5);
figure(4); 
imshow(img_right);
hold on; 
plot([x_R(1) x_R(2)], [y_R1 y_R2], 'r', 'linewidth',2);
title("Right Image"); 

% Choose some point on right image
p_R = [520 250 1]; 
% Find left image epipolar line coefficient for this point
u = p_R*F;

x_R = [1 img_width]; 
y_R1 = -(u(1)*x_R(1) + u(3))/u(2); 
y_R2 = -(u(1)*x_R(2) + u(3))/u(2); 
figure(4);
hold on; 
plot(p_R(1), p_R(2), 'gs', 'markersize', 10, 'linewidth',3); 
title("Left Image"); 
pause(0.5);
figure(3); 
hold on; 
plot([x_R(1) x_R(2)], [y_R1 y_R2], 'g', 'linewidth',2);
title("Right Image"); 
% Finally, repeat the above algorithm for all 20 correspondences and
% plot the epipolar lines on the right image.

% Plot correspondances and epipolar lines for all in the right image
for i = 1:N
    p_L = [xx(i, 1) yy(i, 1) 1]; 
    u = F*p_L';
    x_R = [1 img_width]; 
    y_R1 = -(u(1)*x_R(1) + u(3))/u(2); 
    y_R2 = -(u(1)*x_R(2) + u(3))/u(2); 
    figure(1); 
    hold on; 
    plot(xx(i, 1), yy(i, 1), 'rs', 'markersize', 10, 'linewidth',3); 
    figure(2); 
    hold on; 
    plot(xx(i, 2), yy(i, 2), 'rs', 'markersize', 10, 'linewidth',3); 
    plot([x_R(1) x_R(2)], [y_R1 y_R2], 'linewidth',2);
    hold on; 
    plot([x_R(1) x_R(2)], [y_R1 y_R2], 'linewidth',2);
end

%% Acknowledgements
% The stereo images were acquired by Jerod Weinman in his office and are
% released into the public domain. This algortihm was completed following
% instruction provided by Jerod Weinman, available at
% https://www.cs.grinnell.edu/~weinman/courses/CSC262/2019S/labs/epipolar-geometry.html.
