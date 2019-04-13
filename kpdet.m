function [detection_image] = kpdet(img)
% create gaussian
gauss = gkern(1);

% create gaussian first derivative
gauss1 = gkern(1, 1);

% create big gaussian for blur
gauss_blur = gkern(1.5^2);


% x partial derivative
Ix = conv2(gauss, gauss1, img, 'same');

% y partial derivative
Iy = conv2(gauss1, gauss, img, 'same');


% calculate entries in matrix A
Ix2 = Ix.^2;
Iy2 = Iy.^2;
IxIy = Ix.*Iy;

% create blurred versions of the entries in A
Ix2_blur = conv2(gauss_blur, gauss_blur, Ix2, 'same');
Iy2_blur = conv2(gauss_blur, gauss_blur, Iy2, 'same');
IxIy_blur = conv2(gauss_blur, gauss_blur, IxIy, 'same');

% calculate "image" with each entry containing the determinant of A
img_det = (Ix2_blur.*Iy2_blur-IxIy_blur.^2);

% calculate image with each entry containing the trace of A
img_tr = Ix2_blur + Iy2_blur;

% calculate image containing the ratio of the det and trace
detection_image = img_det ./ img_tr;

threshold = 0.003;
M = maxima(detection_image);
detection_image = M > threshold;

end

