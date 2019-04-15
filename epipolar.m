% Epipolar Geometry
% Getting Correspondences

% Load stereo images 
img_left = im2double(rgb2gray(imread('/home/munozcar/Desktop/left.jpg')));
img_right = im2double(rgb2gray(imread('/home/munozcar/Desktop/right.jpg'))); 

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

[rows1, cols1] = find(features_left);
[rows2, cols2] = find(features_right);

feat_matrix = zeros(size(cols1, 1), 2);
xx = zeros(size(cols1, 1), 3);
yy = zeros(size(cols1, 1), 3);

%% 
%

% Calculating translations for each feature

% Set initial threshold
t = 0.05;
matches = 0;

while matches < 8
matches = 0;
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
    if ((sorted_euclidean_dist_loop(1)/sorted_euclidean_dist_loop(2)) > t)
        feat_matrix(i,1) = NaN;
        feat_matrix(i,2) = NaN;
        xx(i,1) = NaN;
        xx(i,2) = NaN;
        yy(i,1) = NaN;
        yy(i,2) = NaN;
        yy(i,3) = NaN;
        xx(i,3) = NaN;
        continue;
    end
    mindist = sorted_euclidean_dist_loop(1);

    min_ind = find(euclidean_dist_loop == mindist);
    Tc = cols2(min_ind) - cols1(i) ;
    Tr = rows2(min_ind) - rows1(i) ;
    
    if (ismember(cols1(i), xx(:,1)))
        corrind = find(xx(:,1) == cols1(i));
        if (xx(corrind,3) < mindist)
            xx(i,1) = NaN;
            xx(i,2) = NaN;
            yy(i,1) = NaN;
            yy(i,2) = NaN;
            yy(i,3) = NaN;
            xx(i,3) = NaN;
            feat_matrix(i,1) = NaN;
            feat_matrix(i,2) = NaN;
        else % overwrite best match for given feature
            xx(i,1) = cols1(i);
            xx(i,2) = cols2(min_ind);
            yy(i,1) = rows1(i);
            yy(i,2) = rows2(min_ind);
            yy(i,3) = mindist;
            xx(i,3) = mindist;
            feat_matrix(i,1) = Tr;
            feat_matrix(i,2) = Tc;
            xx(corrind,1) = NaN;
            xx(corrind,2) = NaN;
            yy(corrind,1) = NaN;
            yy(corrind,2) = NaN;
            yy(corrind,3) = NaN;
            xx(corrind,3) = NaN;
            feat_matrix(corrind,1) = NaN;
            feat_matrix(corrind,2) = NaN;
        end
    else
        xx(i,1) = cols1(i);
        xx(i,2) = cols2(min_ind);
        yy(i,1) = rows1(i);
        yy(i,2) = rows2(min_ind);
        yy(i,3) = mindist;
        xx(i,3) = mindist;
        feat_matrix(i,1) = Tr;
        feat_matrix(i,2) = Tc;
        
        matches = matches+1;
    end
    
end
t=t+0.02;
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
figure(8);
imshow(stitchedImageMedian, []);
title('Stitched image with median');
pause(0.5);
%%
nonzeroxx = [xx(xx(:,1) ~= 0,1), xx(xx(:,1) ~= 0,2)];
notnanxx = [nonzeroxx(~isnan(nonzeroxx(:,1)),1), nonzeroxx(~isnan(nonzeroxx(:,1)),2)] ;
xx = notnanxx;

nonzeroyy = [yy(yy(:,1) ~= 0,1), yy(yy(:,1) ~= 0,2)];
notnanyy = [nonzeroyy(~isnan(nonzeroyy(:,1)),1), nonzeroyy(~isnan(nonzeroyy(:,1)),2)] ;
yy = notnanyy;

N = size(xx,1);
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
