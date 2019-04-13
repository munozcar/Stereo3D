% Epipolar Geometry
% Getting Correspondences

% Load stereo images 
img_left = im2double(rgb2gray(imread('/home/weinman/courses/CSC262/images/left.jpg')));
img_right = im2double(rgb2gray(imread('/home/weinman/courses/CSC262/images/right.jpg'))); 

features_left = kpdet(img_left);
features_right = kpdet(img_right);

descriptors_left = kpfeat(img_left, features_left);
descriptors_right = kpfeat(img_right, features_right);

[rows1, cols1] = find(features_left);
[rows2, cols2] = find(features_right);

feat_matrix = zeros(size(cols1, 1), 2);

% Calculating translations for each feature
for i=1:size(feat_matrix, 1)
    key = descriptors_left(i,:);
    diff_key = bsxfun(@minus, key, descriptors_right);
    
    % Checking if the feature is NaN 
    if isnan(key(1)) 
        feat_matrix(i,1) = NaN;
        feat_matrix(i,2) = NaN;
        continue;
    end
   % Finding euclidean distances 
    euclidean_dist_loop = sum((diff_key).^2,2);

    % Finding the closest match
    sorted_euclidean_dist_loop = sort(euclidean_dist_loop);
    
    %If the feature doesn't match well enough, set it to nan
    if ((sorted_euclidean_dist_loop(1)/sorted_euclidean_dist_loop(2)) > 0.5)
        feat_matrix(i,1) = NaN;
        feat_matrix(i,2) = NaN;
        continue;
    end
    min_ind = find(euclidean_dist_loop == sorted_euclidean_dist_loop(1));
    Tc = cols2(min_ind) - cols1(i) ;
    Tr = rows2(min_ind) - rows1(i) ;
    
    feat_matrix(i,1) = Tr;
    feat_matrix(i,2) = Tc;
end

%%

rowSortedTrans = sort(feat_matrix(:,1));
colSortedTrans = sort(feat_matrix(:,2));

%Find first NaN
firstNaN = find(isnan(rowSortedTrans));
firstNaN = firstNaN(1);

Tr_median = rowSortedTrans(floor((firstNaN-1)/2));
Tc_median = colSortedTrans(floor((firstNaN-1)/2));

translationMedian = [Tr_median Tc_median];

stitchedImageMedian = stitch(img_right, img_left, translationMedian);
figure;
imshow(stitchedImageMedian, []);
title('Stitched image with median');



%%

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
