function [rectified_l, rectified_r] = kprect(matches_l, matches_r)

% Translation and rotation, focal length
T = [1, 1, 1];
R = ones(2,2);
f = 1;
% New basis
e1 = T./norm(T);
e2 = 1/(sqrt((T(1))^2+(T(2))^2)) .* [-T(2) T(1) 0]';
e3 = cross(e1, e2);

R_rect = [e1' e2' e3']';

Rl = R_rect;
Rr = R*R_rect;

transf_l = Rl*matches_l';
transf_r = Rr*matches_r';

% Rectified points
rectified_l = f*transf_l.*(1/transf_l(3));
rectified_r = f*transf_r.*(1/transf_r(3));



