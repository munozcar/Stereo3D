function [xx, yy, invalid] = kpmatch(feat1, feat2, mops1, mops2)

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
t = 0.05;
matches = 0;
minfeat = 1:min(size(cols1,1), size(cols2,1));
% shuffle order of features
minfeat = randperm(length(minfeat));
features = min(size(cols1,1), size(cols2,1));
disp('Number of features found is');
disp(features);
iter = 0;
invalid = 0;
while matches < 20
matches = 0;
for i=minfeat
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
    mindistratio = sorted_euclidean_dist_loop(1)/sorted_euclidean_dist_loop(2);
    mindist = sorted_euclidean_dist_loop(1);
    min_ind = find(euclidean_dist_loop == mindist,1);
    
    % Easiest to not allow repeated columns at all to avoid double matching
    if (ismember(cols1(i), xx(:,1))) || (ismember(cols2(min_ind), xx(:,2)))
        if (ismember(cols1(i), xx(:,1)))
            corrind = find(xx(:,1) == cols1(i));
        else
            corrind = find(xx(:,2) == cols2(min_ind));
        end
        if (xx(corrind,3) < mindistratio)
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
            yy(i,3) = mindistratio;
            xx(i,3) = mindistratio;
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
        yy(i,3) = mindistratio;
        xx(i,3) = mindistratio;
        
        matches = matches+1;
    end
    
end
t=t+0.005;
minfeat = randperm(length(minfeat));
if t > 0.5
    disp('Not enough matches.');
%     t = 0.05;
%     matches = 0;
%     iter = iter +1;
    invalid = 1;
    break
end
% if iter > 10
%     invalid = 1;
%     disp('killing kpmatch');
%     break
% end
end

%%
nonzeroxx = [xx(xx(:,1) ~= 0,1), xx(xx(:,1) ~= 0,2), xx(xx(:,1) ~= 0,3)];
notnanxx = [nonzeroxx(~isnan(nonzeroxx(:,1)),1), nonzeroxx(~isnan(nonzeroxx(:,1)),2), nonzeroxx(~isnan(nonzeroxx(:,1)),3)] ;
xx = notnanxx;

nonzeroyy = [yy(yy(:,1) ~= 0,1), yy(yy(:,1) ~= 0,2), yy(yy(:,1) ~= 0,3)];
notnanyy = [nonzeroyy(~isnan(nonzeroyy(:,1)),1), nonzeroyy(~isnan(nonzeroyy(:,1)),2), nonzeroyy(~isnan(nonzeroyy(:,1)),3)] ;
yy = notnanyy;


end