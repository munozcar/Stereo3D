%% Epipolar Geometry
%
%% Overview
% We aim to analyze the epipolar geometry describing a pair of images of a
% bookshelf pictured from different reference points. Our goal is to obtain
% corresponding points in both images with which we construct a system
% matrix. Then, we perform singular value decomposition on this matrix in
% order to estimate the fundamental matrix which allows for the full
% reconstruction of the epipolar geometry. We hope that this exercise will
% allow us to increase our knowledge of stereopsis and develop a greater
% intuition for the linear algebra concepts involved in this topic.


%% A. Getting Correspondences
% To begin with, we manually select 20 points of correspondence on two
% images of a bookshelf, shown below. These 20 points are selected such
% that as few as posible are coplanar and we try to span a relatively large
% region.

% Load and show images 
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



%% B. Constructing the System Matrix
% Next, we construct the system matrix S, a 20x9 matrix whose rows are given
% by [xLxR, xLyR, xL, yLxR, yLyR, yL, xR, yR, 1], where xL,xR and yL,yR are
% the x and y coordinates, respectively, of each correspondence point of
% the images. 
corr_matrix = zeros(N,9);
for i=1:N
    corr_matrix(i,:) = [xx(i, 1)*xx(i, 2) xx(i, 1)* yy(i, 2) ...
        xx(i, 1) yy(i, 1)*xx(i, 2) yy(i, 1)*yy(i, 2) ...
        yy(i, 1) xx(i, 2) yy(i, 2) 1]; 
end

%% C. Estimating the Fundamental matrix
% We further apply singular value decomposition to the system matrix $S =
% UWV^t$ in order to estimate the fundamental matrix. For an initial
% estimate, we extract the column of V corresponding to the smallest
% singular value. The rank of this initial estimate of the fundamental
% matrix is 3. Then, we apply singular value decomposition to this initial
% estimate, which yields the singular values: 0.8435, 0.5372, and
% 5.8225e-06. The extremely small magnitude of the third singular value
% relative to the other two is desirable, since this allows us to set it
% equal to zero in order to enforce a singularity and reduce the rank of
% the fundamental matrix. It's low value further lets us know that
% noise and numerical errors by inaccurate correspondances do not play a
% significant role. The first two singular values have magnitudes of
% similar size, as expected, and allow us to reconstruct the
% fundamental matrix. We find that the reconstruction of the adjusted
% fundamental matrix has rank = 2.
% 


[U, S, V] = svd(corr_matrix); 
smallest_col_V = V(:, 8); 
F = reshape(smallest_col_V, [3, 3]); 
rank(F); % This matrix has rank 3

[U_f, S_f, V_f] = svd(F);
S_f(3,3) = 0; 
F = U_f*S_f*V_f'; 
rank(F); % Now, after SVD, it has rank 2

%% D. Visualizing Epipoles
% Now that we have an estimate of the fundamental matrix, we can find the
% epipolar lines for each correspondence. For a point in the left image
% $p_l = (100,50)$, we find the epipolar line coefficients as $u_R=Fp_L^t$.
% The epipolar line in the other image is described by $p_R^tu_R=0$. We do
% the inverse for the point $p_l = (520,250)$ on the right image. The
% figures below show each point and their corresponding epipolar line in
% the other image. In both cases, observe that the epipolar line in one
% images passes through its respective correspondence point in the other
% image. Different colors are provided for convenience.

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
%%
% Finally, we repeat the above algorithm for all 20 correspondences and
% plot the epipolar lines on the right image. The correspondences and lines
% are displayed below. 

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

%%
% Note that all lines go through (or very close) to their respective points
% of correspondence. The epipolar lines further intersect at almost the
% same point: the epipole, which we (manually) determine is at
% approximately x = 276.6737, y = 193.02892. The epipole tells us where, on
% the right image plane, we can find the location at which the camera was
% when it obtained the left image. The perspective we observe in the left
% image further agrees with this result. There is, however, no way to know
% if these pictures were taken at the same time by a pair of cameras or if
% one of the images was taken after the other by moving the camera: There
% is no temporal information in any of the variables determined through
% this exercise.
% 

%% Conclusion
% We succesfully determined the fundamental matrix and found epipolar lines
% for a set of correspondences between two images. We were able to describe
% the epipolar geometry of the given setup and found the epipole with great
% accuracy. We are satisfied with the precision of our values. The results
% from this laboratory build upon our previous work on feature detection
% and matching, and will soon allow us to implement even more complex
% algorithms involving 3-D reconstruction.

%% Acknowledgements
% The stereo images were acquired by Jerod Weinman in his office and are
% released into the public domain. This lab was completed by following the
% instructions provided by Jerod Weinman, available at
% https://www.cs.grinnell.edu/~weinman/courses/CSC262/2019S/labs/epipolar-geometry.html.
% Mathematical concepts for stereo taken from Trucco & Verri, Int. Tech. for
% 3-D Comp. Vision, 1998, Chapter 7. 