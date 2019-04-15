function [xx, yy] = kpmatch(feat1, feat2, mops1, mops2)

%% FUNCTION kpmatch.m
% PURPOSE: kpmatch will match, based on the minimum eucledian distance of
% their corresponding MOPS descriptors, the features of a pair of stereo
% images.
% PARAMETERS: feat1, feat2: NxM matrices with the same size as the images
% to compare, with entries either 1 (for featues) or 0 (absence of
% features).
%             mops1, mops2: 8x8 MOPS descriptors for each of the
%             aforementioned features.
% RETURNS: xx,yy, where xx and yy are 3 by K matrices, with K being the
% number of matches. First and second entries in xx and yy are x and y
% coordinates, respectively, of the matches in each image. Third entry is
% minimum euclidean distance in both cases. 

[rows1, cols1] = find(feat1);
[rows2, cols2] = find(feat2);

xx = zeros(size(cols1, 1), 3);
yy = zeros(size(cols1, 1), 3);


% Calculating translations for each feature

% Set initial threshold
t = 0.2;
matches = 0;

while matches < 8
matches = 0;
for i=1:size(cols1, 1)
    key = mops1(i,:);
    diff_key = bsxfun(@minus, key, mops2);
    
    % Checking if the feature is NaN 
    if isnan(key(1)) 
        continue;
    end
   % Finding euclidean distances 
    euclidean_dist_loop = sum((diff_key).^2,2);
    % Finding the closest match
    sorted_euclidean_dist_loop = sort(euclidean_dist_loop);
    
    %If the feature doesn't match well enough, set it to nan
    if ((sorted_euclidean_dist_loop(1)/sorted_euclidean_dist_loop(2)) > t)
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
    
    % Easiest to not allow repeated columns at all to avoid double matching
    if (ismember(cols1(i), xx(:,1))) || (ismember(cols2(i), xx(:,2)))
        if (ismember(cols1(i), xx(:,1)))
            corrind = find(xx(:,1) == cols1(i));
        else
            corrind = find(xx(:,2) == cols2(i));
        end
        if (xx(corrind,3) < mindist)
            xx(i,1) = NaN;
            xx(i,2) = NaN;
            yy(i,1) = NaN;
            yy(i,2) = NaN;
            yy(i,3) = NaN;
            xx(i,3) = NaN;
        else % overwrite best match for given feature
            xx(i,1) = cols1(i);
            xx(i,2) = cols2(min_ind);
            yy(i,1) = rows1(i);
            yy(i,2) = rows2(min_ind);
            yy(i,3) = mindist;
            xx(i,3) = mindist;
            xx(corrind,1) = NaN;
            xx(corrind,2) = NaN;
            yy(corrind,1) = NaN;
            yy(corrind,2) = NaN;
            yy(corrind,3) = NaN;
            xx(corrind,3) = NaN;
        end
    else
        xx(i,1) = cols1(i);
        xx(i,2) = cols2(min_ind);
        yy(i,1) = rows1(i);
        yy(i,2) = rows2(min_ind);
        yy(i,3) = mindist;
        xx(i,3) = mindist;
        
        matches = matches+1;
    end
    
end
t=t+0.02;
end

%%
nonzeroxx = [xx(xx(:,1) ~= 0,1), xx(xx(:,1) ~= 0,2)];
notnanxx = [nonzeroxx(~isnan(nonzeroxx(:,1)),1), nonzeroxx(~isnan(nonzeroxx(:,1)),2)] ;
xx = notnanxx;

nonzeroyy = [yy(yy(:,1) ~= 0,1), yy(yy(:,1) ~= 0,2)];
notnanyy = [nonzeroyy(~isnan(nonzeroyy(:,1)),1), nonzeroyy(~isnan(nonzeroyy(:,1)),2)] ;
yy = notnanyy;

end