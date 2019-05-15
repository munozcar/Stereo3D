% 3D Streo Reconstructor 

% function to reconstruct a set of 3D points from their 2-view image
% projections and the corresponding camera matrices. Using the mid-point
% method to determine the 3D point as the closest midpoint of the two rays
% this method is fast but not optimal with respect to the reprojection error
%
% inputs: 
% P0, P1 - 3x4 projection matrices
% pts0, pts1 - 2xn corresponding image projections
%
% output: 
% pts3D - 3xn set of Euclidean reconstructed 3D points

function [pts3D] = reconstructor(Pl, Pr, ptsl, ptsr)
  
  % determine the number of points
  npts=size(ptsl, 2);
  
  
  A=zeros(3, 2);
  b=zeros(3, 1);
  pts3D=zeros(3, npts);
  
  Ml=Pl(1:3, 1:3);
  %inverse matrix M0
  Ml=inv(Ml);
  cl=-Ml*Pl(:, 4);
  Mr=Pr(1:3, 1:3);
  Mr=inv(Mr);
  cr=-Mr*Pr(:, 4);
  for i=1:npts
    A=[Ml*[ptsl(:,i);1], -Mr*[ptsr(:,i);1]];
    b=cr-cl;
    [Q, R]=qr(A); a=R\(Q'*b); % LS with QR
    pts3D(:,i)=(cl+a(1)*A(:,1) + cr-a(2)*A(:,2))*0.5;

  end
end

  


  
  
  
