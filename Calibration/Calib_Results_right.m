% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 2051.811082571509814 ; 2053.369763625307314 ];

%-- Principal point:
cc = [ 577.173548522354963 ; 363.547672040730333 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.497754324136366 ; 6.165537267085685 ; -0.005281012073691 ; -0.002195998285051 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 141.603459789803367 ; 139.051248002706842 ];

%-- Principal point uncertainty:
cc_error = [ 63.851610561944298 ; 24.801684867937631 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.227458233038958 ; 3.607000252000328 ; 0.006269679556496 ; 0.012779294765903 ; 0.000000000000000 ];

%-- Image size:
nx = 1280;
ny = 960;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 3;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ -2.217150e+00 ; -2.168293e+00 ; 5.998956e-02 ];
Tc_1  = [ -1.698253e+03 ; -8.673446e+02 ; 6.995745e+03 ];
omc_error_1 = [ 2.283656e-02 ; 2.132900e-02 ; 4.423476e-02 ];
Tc_error_1  = [ 2.184773e+02 ; 8.609081e+01 ; 5.001945e+02 ];

%-- Image #2:
omc_2 = [ 2.037772e+00 ; 1.997903e+00 ; 3.229353e-01 ];
Tc_2  = [ -1.723080e+03 ; -8.564533e+02 ; 7.411220e+03 ];
omc_error_2 = [ 1.943995e-02 ; 1.927157e-02 ; 3.677722e-02 ];
Tc_error_2  = [ 2.317001e+02 ; 9.137735e+01 ; 5.242798e+02 ];

%-- Image #3:
omc_3 = [ -2.078418e+00 ; -2.031317e+00 ; 3.335291e-01 ];
Tc_3  = [ -2.291782e+03 ; -8.238975e+02 ; 9.258466e+03 ];
omc_error_3 = [ 2.343430e-02 ; 2.128296e-02 ; 3.710078e-02 ];
Tc_error_3  = [ 2.887537e+02 ; 1.143245e+02 ; 6.307511e+02 ];

