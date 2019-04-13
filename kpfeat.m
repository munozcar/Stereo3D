function feat_descriptor_matrix = kpfeat(img,features)
%	KPFEAT obtains normalized 8x8 patches from $img at $features locations.
%   	KPFEAT takes in two parameters:
%   	img, an image
%     	features, a logical matix the same size as $img with true values
%                 denoting feature locations
%   	KPFEAT returns:
%     	feat_descriptor_matrix, an count(features) by 64 matrix of mops
%     	patches
%
%   For every feature location described in $features, KPFEAT extracts the
%   surronding 8x8 area from a copy of $img downsampled by a factor of 5,
%   and then preforms gain and bias normalization on the result. These 8x8
%   matrices are then returned in $feat_descriptor_matrix at the index
%   corresponding to its count (left to right, top to bottom) in the
%   features matrix. Any features whos patches overhang the edge will have
%   NaN as an entry instead of a brightness array.

    % Get indices of features
    [feat_x, feat_y] = find(features);
    
    % pre-allocate return value
    feat_descriptor_matrix = zeros(size(feat_x, 1), 64);

    % blur image with gaussian kernel in preparation for downsampling
    GG = gkern(5^2, 0);
    blurred = conv2(GG, GG, img, 'same');
    
    % downsample image
    downsampled = blurred(1:5:end, 1:5:end);

    for k = 1:size(feat_x, 1)
        % Find feature "center" in downsampled image
        row = floor(feat_x(k) / 5);
        column = ceil(feat_y(k) / 5);
        
        % Find edges of 8x8 patch
        ul_x = row - 3;
        ul_y = column - 4;

        br_x = row + 4;
        br_y = column + 3;

        % verify that the patch does not overhang the downsampled image
        if ul_x < 1 || ul_y < 1 || br_x > size(downsampled,1) || ...
                br_y > size(downsampled, 2)
            feat_descriptor_matrix(k,:) = NaN;
            continue;
        end

        % Select patch from downsampled image
        patch = downsampled(ul_x:br_x, ul_y:br_y);
        
        % Normalize for bias
        bias_norm_patch = patch - mean(patch(:));
        
        % Normalize bias patch for gain
        gain_norm_patch = bias_norm_patch ./ std(bias_norm_patch(:));
        
        % Add normalized patch to output
        feat_descriptor_matrix(k,:) = gain_norm_patch(:);
    end
end

