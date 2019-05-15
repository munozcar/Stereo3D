% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 2236.733091852901453 ; 2229.056290370880561 ];

%-- Principal point:
cc = [ 537.287007625079809 ; 353.175369582589440 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.190455394849409 ; 1.285125643207117 ; -0.002898695621085 ; -0.012833126630603 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 163.426516313381171 ; 160.386357815467790 ];

%-- Principal point uncertainty:
cc_error = [ 123.688548685543338 ; 67.067017127869136 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.273868012745139 ; 1.165150953423627 ; 0.006017007535644 ; 0.026076273034279 ; 0.000000000000000 ];

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
omc_1 = [ -2.202449e+00 ; -2.126153e+00 ; 1.281704e-01 ];
Tc_1  = [ 1.132847e+02 ; -8.056755e+02 ; 7.576703e+03 ];
omc_error_1 = [ 3.197009e-02 ; 3.282943e-02 ; 6.304309e-02 ];
Tc_error_1  = [ 4.185347e+02 ; 2.288691e+02 ; 5.517053e+02 ];

%-- Image #2:
omc_2 = [ 2.085269e+00 ; 2.031061e+00 ; 2.596183e-01 ];
Tc_2  = [ 1.373291e+02 ; -7.937626e+02 ; 8.083494e+03 ];
omc_error_2 = [ 3.809150e-02 ; 3.240807e-02 ; 6.105397e-02 ];
Tc_error_2  = [ 4.475214e+02 ; 2.431545e+02 ; 5.906822e+02 ];

%-- Image #3:
omc_3 = [ -2.036480e+00 ; -1.966479e+00 ; 4.381994e-01 ];
Tc_3  = [ -2.812210e+02 ; -7.413428e+02 ; 1.012695e+04 ];
omc_error_3 = [ 4.074580e-02 ; 3.615046e-02 ; 5.862692e-02 ];
Tc_error_3  = [ 5.601279e+02 ; 3.038369e+02 ; 6.933257e+02 ];

