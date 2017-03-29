% Producing the results in the paper for different databases:


% REDLOG 2016 IEEE TIP
% Please make sure to reference the new version of  our paper:
% Reduced-Reference Quality Assessment Based on the Entropy of
% DWT Coefficients of Locally Weighted Gradient Magnitudes
% S. Alireza Golestaneh: sgolest1@asu.edu

% This code was tested on MATLAB2015a-b and 2016a

%=====================================================================


% Copyright Notice:
% Copyright (c) 2016-2017 Arizona Board of Regents.
% All Rights Reserved.
% Contact: Alireza Golestaneh(sgolest1@asu.edu) and Lina Karam (karam@asu.edu)
% Image, Video, and Usabilty (IVU) Lab, ivulab.asu.edu
% Arizona State University
% This copyright statement may not be removed from this file or from
% modifications to this file.
% This copyright notice must also be included in any file or product
% that is derived from this source file.
%
% Redistribution and use of this code in source and binary forms,
% with or without modification, are permitted provided that the
% following conditions are met:
% - Redistribution's of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% - Redistribution's in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% - The Image, Video, and Usability Laboratory (IVU Lab,
% http://ivulab.asu.edu) is acknowledged in any publication that
% reports research results using this code, copies of this code, or
% modifications of this code.
% The code and our papers are to be cited in the bibliography as:
%
%Alireza Golestaneh and Lina Karam. "Reduced-Reference Quality Assessment Based on
%the Entropy of DWT Coefficients of Locally
%Weighted Gradient Magnitudes 2016-TIP.
%
% DISCLAIMER:
% This software is provided by the copyright holders and contributors
% "as is" and any express or implied warranties, including, but not
% limited to, the implied warranties of merchantability and fitness for
% a particular purpose are disclaimed. In no event shall the Arizona
% Board of Regents, Arizona State University, IVU Lab members, or
% contributors be liable for any direct, indirect, incidental, special,
% exemplary, or consequential damages (including, but not limited to,
% procurement of substitute goods or services; loss of use, data, or
% profits; or business interruption) however caused and on any theory
% of liability, whether in contract, strict liability, or tort
% (including negligence or otherwise) arising in any way out of the use
% of this software, even if advised of the possibility of such damage.
%=====================================================================clc
clc
close all
clear all
path(pathdef)


%% % LIVE
clear
load('LIVE.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores,DMOS,zeros(size(REDLOG_Scores)), 0);
disp('LIVE')
[CC SROCC rmse]

%% % CSIQ
clear
load('CSIQ.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores,DMOS,zeros(size(REDLOG_Scores)), 0);
disp('CSIQ')
[CC SROCC rmse]


%% % TID2008
clear
load('TID2008.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores, MOS,zeros(size(REDLOG_Scores)), 1);
disp('TID2008')
[CC SROCC rmse]

%% % TID 2013
clear
load('TID2013.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores, MOS,zeros(size(REDLOG_Scores)), 1);
disp('TID2013')
[CC SROCC rmse]

%% % Toyoma
clear
load('Toyoma.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores, MOS,zeros(size(REDLOG_Scores)), 1);
disp('Toyoma')
[CC SROCC rmse]

%% % IVC
clear
load('IVC.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores, MOS,zeros(size(REDLOG_Scores)), 1);
disp('IVC')
[CC SROCC rmse]

%% % QualTEX
clear
load('QualTEX.mat')
[x_n, CC, SROCC, OR, OD, rmse, residual] =...
    perf_eval_( REDLOG_Scores, MOS,zeros(size(REDLOG_Scores)), 1);
disp('QualTEX')
[CC SROCC rmse]