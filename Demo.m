% REDLOG 2016 IEEE TIP
% Please cite our paper: 
% Reduced-Reference Quality Assessment Based on the Entropy of 
% DWT Coefficients of Locally Weighted Gradient Magnitudes
% S. Alireza Golestaneh: sgolest1@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Copyright Notice:
%  Copyright (c) 2016-2017 Arizona Board of Regents. 
%  All Rights Reserved.
%  Contact: Alireza Golestaneh (sgolest1@asu.edu) and Lina Karam (karam@asu.edu)
%  Image, Video, and Usabilty (IVU) Lab, http://ivulab.asu.edu
%  Arizona State University
%  This copyright statement may not be removed from any file containing it 
% or from modifications to these files.
%  This copyright notice must also be included in any file or product 
%  that is derived from the source files. 
%  
%  Redistribution and use of this code in source and binary forms, 
%  with or without modification, are permitted provided that the 
%  following conditions are met:  
%  - Redistribution's of source code must retain the above copyright 
%    notice, this list of conditions and the following disclaimer. 
%  - Redistribution's in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution. 
%  - The Image, Video, and Usability Laboratory (IVU Lab, 
%    http://ivulab.asu.edu) is acknowledged in any publication that 
%    reports research results using this code, copies of this code, or 
%    modifications of this code.  
%  
% The code and our papers are to be cited in the bibliography as:
% 
% Alireza Golestaneh and Lina Karam. "Reduced-Reference Quality Assessment Based on
% the Entropy of DWT Coefficients of Locally
% Weighted Gradient Magnitudes 2016-TIP.
% 
% 
% DISCLAIMER:
% This software is provided by the copyright holders and contributors 
% "as is" and any express or implied warranties, including, but not 
% limited to, the implied warranties of merchantability and fitness for 
% a particular purpose are disclaimed. In no event shall the Arizona 
% Board of Regents, Arizona State University, IVU Lab members, authors, or 
% contributors be liable for any direct, indirect, incidental, special,
% exemplary, or consequential damages (including, but not limited to, 
% procurement of substitute goods or services; loss of use, data, or 
% profits; or business interruption) however caused and on any theory 
% of liability, whether in contract, strict liability, or tort 
% (including negligence or otherwise) arising in any way out of the use 
% of this software, even if advised of the possibility of such damage. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
path(pathdef)

Org_Img=imread('I01.BMP');
Dist_Img=imread('I01_01_1.bmp');
figure,imshow(Org_Img,[])
figure,imshow(Dist_Img,[])
ShowImagMap=0; % Change to zero if you dont want to see the normalized maps

[REDLOG_Quality_Score Org_Img_Map Dist_Img_Map] = REDLOG_Index(Org_Img,Dist_Img,ShowImagMap);
REDLOG_Quality_Score
disp('Please note that you should receive REDLOG_Quality_Score=5.2526 for Demo, otherwise it means there is something wrong in your settings!')