% Epipolar Geometry

% Getting Correspondences

% Load and show stereo images 
img_left = imread("/home/weinman/courses/CSC262/images/left.jpg");
img_right = imread("/home/weinman/courses/CSC262/images/right.jpg"); 

figure(1);
imshow(img_left);
title('Left Image');
pause(0.5);

figure(2);
imshow(img_right);
title('Right Image');
pause(0.5);

N = 20;
% Check if points already selected, is so, load coordinates
% Make sure corr.mat file is in same directory as epipolar.m
if isfile("corr.mat")
    load 'corr.mat';
    for i = 1 : N
        figure(1); 
        hold on; 
        plot(xx(i,1), yy(i,1)); 
        figure(2); 
        hold on; 
        plot(xx(i,2), yy(i,2));
    end
else
    % Loop allows user to select N correspondance points 
    xx = zeros(N,2);
    yy = zeros(N,2);

    for i = 1 : N
       figure(1); 
       [xx(i,1), yy(i,1)] = ginput(1); 
       hold on; 
       plot(xx(i,1), yy(i,1)); 
       figure(2); 
       [xx(i,2), yy(i,2)] = ginput(1);
    end
    save('corr.mat', 'xx', 'yy'); 
end



% Constructing the System Matrix
% Next, construct the system matrix S, a 20x9 matrix whose rows are given
% by [xLxR, xLyR, xL, yLxR, yLyR, yL, xR, yR, 1], where xL,xR and yL,yR are
% the x and y coordinates, respectively, of each correspondence point of
% the images. 
corr_matrix = zeros(N,9);
for i=1:N
    corr_matrix(i,:) = [xx(i, 1)*xx(i, 2) xx(i, 1)* yy(i, 2) ...
        xx(i, 1) yy(i, 1)*xx(i, 2) yy(i, 1)*yy(i, 2) ...
        yy(i, 1) xx(i, 2) yy(i, 2) 1]; 
end

% Estimating the Fundamental matrix
% Apply singular value decomposition to the system matrix $S =
% UWV^t$ in order to estimate the fundamental matrix. For an initial
% estimate, extract the column of V corresponding to the smallest
% singular value. The rank of this initial estimate of the fundamental
% matrix is 3. Then, apply singular value decomposition to this initial
% estimate. The extremely small magnitude of the third singular value
% relative to the other two is desirable, since this allows us to set it
% equal to zero in order to enforce a singularity and reduce the rank.

[U, S, V] = svd(corr_matrix); 
smallest_col_V = V(:, 8); 
F = reshape(smallest_col_V, [3, 3]); 
rank(F); % This matrix has rank 3

[U_f, S_f, V_f] = svd(F);
S_f(3,3) = 0; 
F = U_f*S_f*V_f'; 
rank(F); % Now, after SVD, it has rank 2

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
