function [F] = kpfund(xx,yy)

N = size(xx,1); % number of matches

% Construct the system matrix S, a 20x9 matrix whose rows are given
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
end