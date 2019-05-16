function [set3D] = reconstructor(Pl, Pr, ptl, ptr)

% FUNCTION: reconstructor.
% PARAMETERS:
%               - Pl: left camera matrix (3x4)
%               - Pr: right camera matrix (3x4)
%               - ptl: left image projections (2xn)
%               - ptr: right image projections (2xn)
% PRODUCES: 
%               -set3D: set of Euclidean reconstructed 3D points (3xn)
%
% PURPOSE:       function to reconstruct a set of 3D points from
%                their 2-view image projections and the corresponding
%                camera matrices. Utilizes the mid-point method.
%%
  % determine the number n of correspondents
  n = size(ptl, 2);
  
  % initialize the values
  A = zeros(3, 2);
  b = zeros(3, 1);
  set3D = zeros(3, n);
  
  % extract 3x3 rotation matrix Rl from Pl
  Rl = Pl(1:3, 1:3);
  Tl = Pl(1:3, 4);
  % invert Rl
  Rl=inv(Rl);
  % calculate cl from Rl and Pl
  cl=-Rl*Pl(:, 4);
  
  % extract 3x3 rotation matrix Mr from Pr
  Rr=Pr(1:3, 1:3);
  Tr = Pr(1:3, 4);
  % invert matrix Rr
  Rr=inv(Rr);
  % calculate cr from Mr and Pr
  cr=-Rr*Pr(:, 4);
  
  % calculate the rotational matrix
  %R = Rr*Rl';
  % calculate the translational vector
  %T = Tl - R'*Tr;
  
  % calculate the 3D point of each pair of correspondent.
  for i=1:n
    A=[Rl*[ptl(:,i);1], -Rr*[ptr(:,i);1]];
    % solve for scale factor b
    b=cr-cl;
    % Orthogonal-triangular decomposition of A
    [Q, R]=qr(A); 
    % solve for scale factor a
    a = R\(Q'*b); 
    set3D(:,i)=(cl+a(1)*A(:,1) + cr-a(2)*A(:,2))*0.5;

  end
end



  


  
  
  
