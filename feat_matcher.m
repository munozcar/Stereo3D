function [xx, yy] = feat_matcher(imgl, imgr)
% FUNCTION: feat_matcher, a fast feature matcher.
% PARAMETERS:
%               - imgl: left stereo image 
%               - imgr: right stereo image
% PRODUCES: 
%               - xx: matched K features in the left image. (Kx2
%               array)
%               - yy: matched K features in the right image.
%               (Kx2 array)
% PURPOSE:      Detects key features in a pair of stereo images by
%               computing the harmonic mean. Matches said features using a
%               nearest neighbor kd-tree. Finally, estimates the
%               fundamental matrix of the system.

%% KEY FEATURE DETECTION

matches = 0;

disp("Finding features...");

if size(imgl,3)>1
    imgl = im2double(rgb2gray(imgl));
    imgr = im2double(rgb2gray(imgr));
end

% Initial detection threshold
threshold = 0.000005;

% Create Gaussian filters 
gauss = gkern(1); gauss1 = gkern(1, 1); gauss_blur = gkern(1.5^2);

% Compute x partial derivative of imgl
Ixl = conv2(gauss, gauss1, imgl, 'same');
% Compute y partial derivative of imgl
Iyl = conv2(gauss1, gauss, imgl, 'same');

% Compute x partial derivative of imgr
Ixr = conv2(gauss, gauss1, imgr, 'same');
% Compute y partial derivative of imgr
Iyr = conv2(gauss1, gauss, imgr, 'same');

% Calculate entries in matrix A (left)
Ix2l = Ixl.^2; Iy2l = Iyl.^2; IxIyl = Ixl.*Iyl;

% Calculate entries in matrix A (right)
Ix2r = Ixr.^2; Iy2r = Iyr.^2; IxIyr = Ixr.*Iyr;

% create blurred versions of the entries in A for both images
Ix2l_blur = conv2(gauss_blur, gauss_blur, Ix2l, 'same');
Iy2l_blur = conv2(gauss_blur, gauss_blur, Iy2l, 'same');
IxIyl_blur = conv2(gauss_blur, gauss_blur, IxIyl, 'same');

Ix2r_blur = conv2(gauss_blur, gauss_blur, Ix2r, 'same');
Iy2r_blur = conv2(gauss_blur, gauss_blur, Iy2r, 'same');
IxIyr_blur = conv2(gauss_blur, gauss_blur, IxIyr, 'same');


% calculate "image" with each entry containing the determinant of A
imgl_det = (Ix2l_blur.*Iy2l_blur-IxIyl_blur.^2);
imgr_det = (Ix2r_blur.*Iy2r_blur-IxIyr_blur.^2);

% calculate image with each entry containing the trace of A
imgl_tr = Ix2l_blur + Iy2l_blur;
imgr_tr = Ix2r_blur + Iy2r_blur;

% calculate image containing the ratio of the det and trace
det_imgl = imgl_det ./ imgl_tr; Ml = maxima(det_imgl);
det_imgr = imgr_det ./ imgr_tr; Mr = maxima(det_imgr);

% Use only local maxima and determine the detections
det_imgl = det_imgl.*Ml; det_imgl = det_imgl>threshold;
det_imgr = det_imgr.*Mr; det_imgr = det_imgr>threshold;

%% DESCRIPTORS
disp("Building simple feature descriptors...");

% Get indices of features on both images
[rowsl, colsl] = find(det_imgl); [rowsr, colsr] = find(det_imgr);

% blur image with gaussian kernel in preparation for downsampling
GG = gkern(5^2, 0);
blurredl = conv2(GG, GG, imgl, 'same');
blurredr = conv2(GG, GG, imgr, 'same');

% downsample images
downsampledl = blurredl(1:5:end, 1:5:end);
downsampledr = blurredr(1:5:end, 1:5:end);

% Pre-allocate space for tensor of 7x7 descriptors
SSDl = zeros(size(rowsl, 1), 49);
SSDr = zeros(size(rowsr, 1), 49);

% Use minimum number of features
minfeat = min(size(rowsl,1), size(rowsr,1));

% Generate left descriptors
for k = 1:size(rowsl,1)
    row = floor(rowsl(k)/5); column = floor(colsl(k)/5);
    % Find edges of patch
    lcol = column-3; rcol = column+3; drow = row-3; urow = row+3;
    % verify that the patch does not overhang the downsampled image
    if lcol < 1 || drow < 1 || urow > size(downsampledl,1) || ...
            rcol > size(downsampledl, 2)
        SSDl(k,:) = NaN;
        continue;
    end
    % Select patch from image
    patch = downsampledl(drow:urow,lcol:rcol);
    % Normalize for bias
    bias_norm_patch = patch - mean(patch(:));
    % Normalize bias patch for gain
    gain_norm_patch = bias_norm_patch ./ std(bias_norm_patch(:));
    patch = gain_norm_patch;
    % Add patch to output
    SSDl(k,:) = patch(:);
end

% Generate left descriptors
for k = 1:size(rowsr,1)
    row = floor(rowsr(k)/5); column = floor(colsr(k)/5);
    % Find edges of patch
    lcol = column-3; rcol = column+3; drow = row-3; urow = row+3;
    % verify that the patch does not overhang the downsampled image
    if lcol < 1 || drow < 1 || urow > size(downsampledr,1) || ...
            rcol > size(downsampledr, 2)
        SSDr(k,:) = NaN;
        continue;
    end
    % Select patch from image
    patch = downsampledr(drow:urow,lcol:rcol);
    % Normalize for bias
    bias_norm_patch = patch - mean(patch(:));
    % Normalize bias patch for gain
    gain_norm_patch = bias_norm_patch ./ std(bias_norm_patch(:));
    patch = gain_norm_patch;
    % Add patch to output
    SSDr(k,:) = patch(:);
end

rowsl(any(isnan(SSDl),2),:)=[];
colsl(any(isnan(SSDl),2),:)=[];
SSDl(any(isnan(SSDl),2),:)=[];

rowsr(any(isnan(SSDr),2),:)=[];
colsr(any(isnan(SSDr),2),:)=[];
SSDr(any(isnan(SSDr),2),:)=[];

minfeat = min(size(rowsl,1), size(rowsr,1));

%% MATCH FEATURES

disp("Computing nearest neighbor tree...");

% Compute nearest neighbors once 
[IDX,D] = knnsearch(SSDl, SSDr, 'K',2,'Distance','euclidean');
IDX = IDX(1:minfeat,1);
D = D(1:minfeat,:);

% Mathched features
xxr = colsr(1:minfeat); xxl = colsl(IDX); 
yyr = rowsr(1:minfeat); yyl = rowsl(IDX);

for i=1:minfeat
    ratio = D(i,1)/D(i,2);
    if ratio>0.35
        xxr(i) = NaN; xxl(i) = NaN; yyr(i) = NaN; yyl(i) = NaN;
    else
        continue;
    end
end

% Remove nans
xxr = xxr(~isnan(xxr)); xxl = xxl(~isnan(xxl)); 
yyr = yyr(~isnan(yyr)); yyl = yyl(~isnan(yyl));

matches = size(xxr,1);

if matches<8
    disp("Not enough matches found. Try lowering the threshold.");
    xx = NaN;
    yy = NaN;
else
    xx = [xxl xxr];
    yy = [yyl yyr];
    disp("Matches found. Done.");
end

end