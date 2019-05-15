function [F,el,er] = epipolar(xx,yy)

% FUNCTION: epipolar
% PARAMETERS:
%               - xx: matched features' x coordinates in left and right
%               images. (Nx2 matrix, where N is the number of matches)
%               - yy: matched features' y coordinates in left and right
%               images.
% PRODUCES: 
%               - F: An estimate of the fundamental matrix
%               - el: left epipole
%               - er: right epipole
% PURPOSE:      Determines the epipolar geometry of a stereo scene, given
%               enough (>8) non-coplanar matching features. To avoid
%               numerical innacuracy, the coordinates are isotropically
%               normalized.
% Last edited: May 15, Carlos Munoz


N = size(xx,1); % number of matches

% Construct the system matrix S, a matrix whose rows are given
% by [xLxR, xLyR, xL, yLxR, yLyR, yL, xR, yR, 1], where xL,xR and yL,yR are
% the x and y coordinates, respectively, of each correspondence point of
% the images. 

% Normalize coordinates with isotropic scaling ( based on
% https://www.ecse.rpi.edu/~qji/CV/8point.pdf, by Philip Lamoureux)
% This extra process gets rid of numerical instabilities and helps reduce
% the effects of outliers.

xbar_l = mean(xx(:,1));
ybar_l = mean(yy(:,1));
xbar_r = mean(xx(:,2));
ybar_r = mean(yy(:,2));

d_l = 1/sqrt(2) * mean(sqrt((xx(:,1)-xbar_l).^2  + (yy(:,1)-ybar_l).^2 )); 
d_r = 1/sqrt(2) * mean(sqrt((xx(:,2)-xbar_r).^2  + (yy(:,2)-ybar_r).^2 )); 

IsoT_l = [1/d_l 0 0; 0 1/d_l 0; -xbar_l/d_l -ybar_l/d_l 1]';
IsoT_r = [1/d_r 0 0; 0 1/d_r 0; -xbar_r/d_r -ybar_r/d_r 1]';

for i=1:N
    new_l = IsoT_l*[xx(i,1) yy(i,1) 1]';
    new_r = IsoT_r*[xx(i,2) yy(i,2) 1]';

    xx(i,1) = new_l(1); yy(i,1) = new_l(2);
    xx(i,2) = new_r(1); yy(i,2) = new_r(2);
end

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

[~, S, V] = svd(corr_matrix);   
smallest_col_V = V(:,end); 
F = reshape(smallest_col_V, [3, 3]); 

[U_f, S_f, V_f] = svd(F);
S_f(end,end)=0;
F_estimate = U_f*S_f*V_f'; 

F = IsoT_r'*F_estimate*IsoT_l;
% Epipoles
[U, D, V] = svd(F);
el = V(:,3);
er = U(:,3);

end