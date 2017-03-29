% REDLOG 2016 IEEE TIP
% Please cite our paper in any published work if you use this code.
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

function [REDLOG_Score Img_Grad_Nomalized_O Img_Grad_Nomalized_D]=REDLOG_Index(Org_Img,Dist_Img,showimage )
% Pooling
Num_or = 4;
Num_scales = 6;
[REDLOG_Features_O Img_Grad_Nomalized_O]=REDLOG_Features_Computing(Org_Img,showimage,Num_or,Num_scales);
  [REDLOG_Features_D Img_Grad_Nomalized_D]=REDLOG_Features_Computing(Dist_Img,showimage,Num_or,Num_scales);
  
Pooled = [(sum(abs(REDLOG_Features_O(1:Num_scales)-REDLOG_Features_D(1:Num_scales)).^2)) REDLOG_Features_D(Num_scales+1) ];
REDLOG_Score=log(((Pooled(:,1).*(Pooled(:,2)+1)*50)+1));
end

function [REDLOG_Features Img_Grad_Nomalized]=REDLOG_Features_Computing(Input_Img,showimage,Num_or,Num_scales)
 

Size=size(Input_Img);
if length(Size)==3; Input_Img=rgb2gray(Input_Img);  end

Input_Img=double(Input_Img);

% Making CSF
csf = make_csf( Size(1), Size(2), 32 )';
%CSF of Image
Img_CSF = double(real( ifft2( ifftshift( fftshift( fft2( Input_Img ) ).* csf ) ) ));

BL=8;
sigm= 0.5;
Gauss_window  = fspecial('gaussian',BL,sigm);
Img_CSF=filter2(Gauss_window,Img_CSF );

% Gradient
[Img_Grad,Img_GDir] = imgradient_MATLAB( Img_CSF);
% Normalizing Gradient
A= (sum(cat(3,(Img_Grad), Img_CSF).^2,3))./2;

% Normaliztion
window1 = fspecial('gaussian',BL,sigm);
window1 = window1/sum(window1(:));
Alpha =  sqrt(filter2(window1,A,'same'));
FAC=Alpha<256/10;FAC=256/10.*FAC;
Alpha=Alpha+FAC;
Img_Grad_Nomalized = ((Img_Grad)./Alpha);
if showimage==1
    figure,imshow(Img_Grad_Nomalized,[])
end


% Wavelet Coef
Img_Input_Wavelet=  Img_Grad_Nomalized;
[Img_Pyr Img_Pind] = buildSFpyr(Img_Input_Wavelet,Num_scales,Num_or-1);
[Img_subbandDNT Img_subbandGSM OrgImg_Cov_Mat OrgImg_size_band] =...
    DWTCOEF(Img_Pyr,Img_Pind,Num_scales,Num_or,1,1,3,3,100);
clc
ScaleNum=0;
for ii = 1:Num_or: length(Img_subbandGSM)
    OrtNumber=0;
    ScaleNum=ScaleNum+1;
    for iii=ii:ii+Num_or -1
        OrtNumber=OrtNumber+1;
        Img_Coeff = Img_subbandGSM{iii};
        Param1(1,OrtNumber) = log(1+entropy(uint8((Img_Coeff))));
        Param2(1,OrtNumber) =  (mean2(abs(Img_Coeff(:))));
    end
    SumFeatur1(1,ScaleNum)=sum(abs(Param1));
    SumFeatur2(1,ScaleNum)=sum(abs(Param2));
    
end
%  RR features
REDLOG_Features=[];
REDLOG_Features=[SumFeatur1 sum(SumFeatur2)];
end

function [res] = make_csf(x, y, nfreq)
[xplane,yplane]=meshgrid(-x/2+0.5:x/2-0.5, -y/2+0.5:y/2-0.5);	% generate mesh
plane=(xplane+1i*yplane)/y*2*nfreq;
radfreq=abs(plane);				% radial frequency

w=0.7;
s=(1-w)/2*cos(4*angle(plane))+(1+w)/2;
radfreq=radfreq./s;

csf = 2.6*(0.0192+0.114*radfreq).*exp(-(0.114*radfreq).^1.1);
f=find( radfreq < 7.8909 ); csf(f)=0.9809+zeros(size(f));
res = csf;
end


% Following functions are from Eero Simoncelli
% [PYR, INDICES] = buildSFpyrLevs(LODFT, LOGRAD, XRCOS, YRCOS, ANGLE, HEIGHT, NBANDS)
% Recursive function for constructing levels of a steerable pyramid.  This
% is called by buildSFpyr, and is not usually called directly.
% Eero Simoncelli, 5/97.
function [pyr,pind] = buildSFpyrLevs(lodft,log_rad,Xrcos,Yrcos,angle,ht,nbands);

if (ht <= 0)
    lo0 = ifft2(ifftshift(lodft));
    pyr = real(lo0(:));
    pind = size(lo0);
else
    bands = zeros(prod(size(lodft)), nbands);
    bind = zeros(nbands,2);
    %  log_rad = log_rad + 1;
    Xrcos = Xrcos - log2(2);  % shift origin of lut by 1 octave.
    lutsize = 1024;
    Xcosn = pi*[-(2*lutsize+1):(lutsize+1)]/lutsize;  % [-2*pi:pi]
    order = nbands-1;
    %% divide by sqrt(sum_(n=0)^(N-1)  cos(pi*n/N)^(2(N-1)) )
    %% Thanks to Patrick Teo for writing this out :)
    const = (2^(2*order))*(factorial(order)^2)/(nbands*factorial(2*order));
    Ycosn = sqrt(const) * (cos(Xcosn)).^order;
    himask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
    
    for b = 1:nbands
        anglemask = pointOp(angle, Ycosn, Xcosn(1)+pi*(b-1)/nbands, Xcosn(2)-Xcosn(1));
        banddft = ((-sqrt(-1))^(nbands-1)) .* lodft .* anglemask .* himask;
        band = ifft2(ifftshift(banddft));
        
        bands(:,b) = real(band(:));
        bind(b,:)  = size(band);
    end
    dims = size(lodft);
    ctr = ceil((dims+0.5)/2);
    lodims = ceil((dims-0.5)/2);
    loctr = ceil((lodims+0.5)/2);
    lostart = ctr-loctr+1;
    loend = lostart+lodims-1;
    log_rad = log_rad(lostart(1):loend(1),lostart(2):loend(2));
    angle = angle(lostart(1):loend(1),lostart(2):loend(2));
    lodft = lodft(lostart(1):loend(1),lostart(2):loend(2));
    YIrcos = abs(sqrt(1.0 - Yrcos.^2));
    lomask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
    lodft = lomask .* lodft;
    [npyr,nind] = buildSFpyrLevs(lodft, log_rad, Xrcos, Yrcos, angle, ht-1, nbands);
    pyr = [bands(:); npyr];
    pind = [bind; nind];
end
end

% [PYR, INDICES, STEERMTX, HARMONICS] = buildSFpyr(IM, HEIGHT, ORDER, TWIDTH)
% Construct a steerable pyramid on matrix IM, in the Fourier domain.
% This is similar to buildSpyr, except that:
%    + Reconstruction is exact (within floating point errors)
%    + It can produce any number of orientation bands.
%    - Typically slower, especially for non-power-of-two sizes.
%    - Boundary-handling is circular.
% HEIGHT (optional) specifies the number of pyramid levels to build. Default
% is maxPyrHt(size(IM),size(FILT));
% The squared radial functions tile the Fourier plane, with a raised-cosine
% falloff.  Angular functions are cos(theta-k\pi/(K+1))^K, where K is
% the ORDER (one less than the number of orientation bands, default= 3).
% TWIDTH is the width of the transition region of the radial lowpass
% function, in octaves (default = 1, which gives a raised cosine for
% the bandpass filters).
% PYR is a vector containing the N pyramid subbands, ordered from fine
% to coarse.  INDICES is an Nx2 matrix containing the sizes of
% each subband.  This is compatible with the MatLab Wavelet toolbox.
% See the function STEER for a description of STEERMTX and HARMONICS.
% Eero Simoncelli, 5/97.
% See http://www.cns.nyu.edu/~eero/STEERPYR/ for more
% information about the Steerable Pyramid image decomposition.

function [pyr,pind,steermtx,harmonics] = buildSFpyr(im, ht, order, twidth)
%-----------------------------------------------------------------
%% DEFAULTS:
max_ht = floor(log2(min(size(im)))) - 2;

if (exist('ht') ~= 1)
    ht = max_ht;
else
    if (ht > max_ht)
        error(sprintf('Cannot build pyramid higher than %d levels.',max_ht));
    end
end

if (exist('order') ~= 1)
    order = 3;
elseif ((order > 15)  | (order < 0))
    fprintf(1,'Warning: ORDER must be an integer in the range [0,15]. Truncating.\n');
    order = min(max(order,0),15);
else
    order = round(order);
end
nbands = order+1;

if (exist('twidth') ~= 1)
    twidth = 1;
elseif (twidth <= 0)
    fprintf(1,'Warning: TWIDTH must be positive.  Setting to 1.\n');
    twidth = 1;
end

%-----------------------------------------------------------------
%% Steering stuff:

if (mod((nbands),2) == 0)
    harmonics = [0:(nbands/2)-1]'*2 + 1;
else
    harmonics = [0:(nbands-1)/2]'*2;
end

steermtx = steer2HarmMtx(harmonics, pi*[0:nbands-1]/nbands, 'even');

%-----------------------------------------------------------------

dims = size(im);
ctr = ceil((dims+0.5)/2);

[xramp,yramp] = meshgrid( ([1:dims(2)]-ctr(2))./(dims(2)/2), ...
    ([1:dims(1)]-ctr(1))./(dims(1)/2) );
angle = atan2(yramp,xramp);
log_rad = sqrt(xramp.^2 + yramp.^2);
log_rad(ctr(1),ctr(2)) =  log_rad(ctr(1),ctr(2)-1);
log_rad  = log2(log_rad);

%% Radial transition function (a raised cosine in log-frequency):
[Xrcos,Yrcos] = rcosFn(twidth,(-twidth/2),[0 1]);
Yrcos = sqrt(Yrcos);

YIrcos = sqrt(1.0 - Yrcos.^2);
lo0mask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
imdft = fftshift(fft2(im));
lo0dft =  imdft .* lo0mask;

[pyr,pind] = buildSFpyrLevs(lo0dft, log_rad, Xrcos, Yrcos, angle, ht, nbands);

hi0mask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
hi0dft =  imdft .* hi0mask;
hi0 = ifft2(ifftshift(hi0dft));

pyr = [real(hi0(:)) ; pyr];
pind = [size(hi0); pind];
end


function [subbandDNT subbandGSM Cov_Mat size_band] = DWTCOEF(pyro,pind,Nsc,Nor,parent,neighbor,blSzX,blSzY,nbins)
pyro = real(pyro);
Nband = size(pind,1)-1;
band = [1 3 6 8 9 11];
zRange = 15;

p = 1;
for scale=1:Nsc
    for orien=1:Nor
        nband = (scale-1)*Nor+orien+1; % except the ll
        aux = pyrBand(pyro, pind, nband);
        [Nsy,Nsx] = size(aux);
        prnt = parent & (nband < Nband-Nor);   % has the subband a parent?
        BL = zeros(size(aux,1),size(aux,2),1 + prnt);
        BL(:,:,1) = aux;
        if prnt,
            auxp = pyrBand(pyro, pind, nband+Nor);
            auxp = real(imresize(auxp,2)); %
            %  fprintf('parent band and size is %d %d %d \n',nband+Nor,Nsy,Nsx)
            BL(:,:,2) = auxp(1:Nsy,1:Nsx);
        end
        y=BL;
        [nv,nh,nb] = size(y);
        block = [blSzX blSzY];
        %
        nblv = nv-block(1)+1;	% Discard the outer coefficients
        nblh = nh-block(2)+1;   % for the reference (centrral) coefficients (to avoid boundary effects)
        nexp = nblv*nblh;			% number of coefficients considered
        Ly = (block(1)-1)/2;		% block(1) and block(2) must be odd!
        Lx = (block(2)-1)/2;
        Cov_Mat{p} = 0;
        o_c = aux(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
        o_c = (o_c(:));
        o_c = o_c - mean(o_c);
        subbandGSM{p} = o_c;
        subbandDNT=0;
        size_band(p,:) = 0;
        p = p+1;
    end
end
end

% MTX = steer2HarmMtx(HARMONICS, ANGLES, REL_PHASES)
% Compute a steering matrix (maps a directional basis set onto the
% angular Fourier harmonics).  HARMONICS is a vector specifying the
% angular harmonics contained in the steerable basis/filters.  ANGLES
% (optional) is a vector specifying the angular position of each filter.
% REL_PHASES (optional, default = 'even') specifies whether the harmonics
% are cosine or sine phase aligned about those positions.
% The result matrix is suitable for passing to the function STEER.
% Eero Simoncelli, 7/96.
function mtx = steer2HarmMtx(harmonics, angles, evenorodd)

%%=================================================================
%%% Optional Parameters:

if (exist('evenorodd') ~= 1)
    evenorodd = 'even';
end

% Make HARMONICS a row vector
harmonics = harmonics(:)';

numh = 2*size(harmonics,2) - any(harmonics == 0);

if (exist('angles') ~= 1)
    angles = pi * [0:numh-1]'/numh;
else
    angles = angles(:);
end

%%=================================================================

if isstr(evenorodd)
    if strcmp(evenorodd,'even')
        evenorodd = 0;
    elseif strcmp(evenorodd,'odd')
        evenorodd = 1;
    else
        error('EVEN_OR_ODD should be the string  EVEN or ODD');
    end
end

%% Compute inverse matrix, which maps Fourier components onto
%% steerable basis.
imtx = zeros(size(angles,1),numh);
col = 1;
for h=harmonics
    args = h*angles;
    if (h == 0)
        imtx(:,col) = ones(size(angles));
        col = col+1;
    elseif evenorodd
        imtx(:,col) = sin(args);
        imtx(:,col+1) = -cos(args);
        col = col+2;
    else
        imtx(:,col) = cos(args);
        imtx(:,col+1) = sin(args);
        col = col+2;
    end
end

r = rank(imtx);
if (( r ~= numh ) & ( r ~= size(angles,1) ))
    fprintf(2,'WARNING: matrix is not full rank');
end

mtx = pinv(imtx);

end


% [RES] = shift(MTX, OFFSET)
% Circular shift 2D matrix samples by OFFSET (a [Y,X] 2-vector),
% such that  RES(POS) = MTX(POS-OFFSET).

function res = shift(mtx, offset)

dims = size(mtx);

offset = mod(-offset,dims);

res = [ mtx(offset(1)+1:dims(1), offset(2)+1:dims(2)),  ...
    mtx(offset(1)+1:dims(1), 1:offset(2));        ...
    mtx(1:offset(1), offset(2)+1:dims(2)),          ...
    mtx(1:offset(1), 1:offset(2)) ];
end

% RES = pyrBandIndices(INDICES, BAND_NUM)
%
% Return indices for accessing a subband from a pyramid
% (gaussian, laplacian, QMF/wavelet, steerable).
% Eero Simoncelli, 6/96.
function indices =  pyrBandIndices(pind,band)

if ((band > size(pind,1)) | (band < 1))
    error(sprintf('BAND_NUM must be between 1 and number of pyramid bands (%d).', ...
        size(pind,1)));
end

if (size(pind,2) ~= 2)
    error('INDICES must be an Nx2 matrix indicating the size of the pyramid subbands');
end

ind = 1;
for l=1:band-1
    ind = ind + prod(pind(l,:));
end

indices = ind:ind+prod(pind(band,:))-1;
end

% RES = pyrBand(PYR, INDICES, BAND_NUM)
% Access a subband from a pyramid (gaussian, laplacian, QMF/wavelet,
% or steerable).  Subbands are numbered consecutively, from finest
% (highest spatial frequency) to coarsest (lowest spatial frequency).
% Eero Simoncelli, 6/96.

function res =  pyrBand(pyr, pind, band)

res = reshape( pyr(pyrBandIndices(pind,band)), pind(band,1), pind(band,2) );
end

% RES = pointOp(IM, LUT, ORIGIN, INCREMENT, WARNINGS)
%
% Apply a point operation, specified by lookup table LUT, to image IM.
% LUT must be a row or column vector, and is assumed to contain
% (equi-spaced) samples of the function.  ORIGIN specifies the
% abscissa associated with the first sample, and INCREMENT specifies the
% spacing between samples.  Between-sample values are estimated via
% linear interpolation.  If WARNINGS is non-zero, the function prints
% a warning whenever the lookup table is extrapolated.
% This function is much faster than MatLab's interp1, and allows
% extrapolation beyond the lookup table domain.  The drawbacks are
% that the lookup table must be equi-spaced, and the interpolation is
% linear.
% Eero Simoncelli, 8/96.

function res = pointOp(im, lut, origin, increment, warnings)

%% NOTE: THIS CODE IS NOT ACTUALLY USED! (MEX FILE IS CALLED INSTEAD)

% fprintf(1,'WARNING: You should compile the MEX version of "pointOp.c",\n         found in the MEX subdirectory of matlabPyrTools, and put it in your matlab path.  It is MUCH faster.\n');

X = origin + increment*[0:size(lut(:),1)-1];
Y = lut(:);

res = reshape(interp1(X, Y, im(:), 'linear', 'extrap'),size(im));
end

% RES = innerProd(MTX)
% Compute (MTX' * MTX) efficiently (i.e., without copying the matrix)
function res = innerProd(mtx)

% fprintf(1,['WARNING: You should compile the MEX version of' ...
% 	   ' "innerProd.c",\n         found in the MEX subdirectory' ...
% 	   ' of matlabPyrTools, and put it in your matlab path.' ...
% 	   ' It is MUCH faster and requires less memory.\n']);
res = mtx' * mtx;
end

% [X, Y] = rcosFn(WIDTH, POSITION, VALUES)
% Return a lookup table (suitable for use by INTERP1)
% containing a "raised cosine" soft threshold function:
%    Y =  VALUES(1) + (VALUES(2)-VALUES(1)) *
%              cos^2( PI/2 * (X - POSITION + WIDTH)/WIDTH )
% WIDTH is the width of the region over which the transition occurs
% (default = 1). POSITION is the location of the center of the
% threshold (default = 0).  VALUES (default = [0,1]) specifies the
% values to the left and right of the transition.
% Eero Simoncelli, 7/96.

function [X, Y] = rcosFn(width,position,values)
%------------------------------------------------------------
% OPTIONAL ARGS:
if (exist('width') ~= 1)
    width = 1;
end
if (exist('position') ~= 1)
    position = 0;
end

if (exist('values') ~= 1)
    values = [0,1];
end
%------------------------------------------------------------
sz = 256;  %% arbitrary!

X    = pi * [-sz-1:1] / (2*sz);

Y = values(1) + (values(2)-values(1)) * cos(X).^2;

%    Make sure end values are repeated, for extrapolation...
Y(1) = Y(2);
Y(sz+3) = Y(sz+2);

X = position + (2*width/pi) * (X + pi/4);
end

function [Gmag, Gdir] = imgradient_MATLAB(varargin) 
narginchk(1,2);

[I, Gx, Gy, method] = parse_inputs(varargin{:});

% Compute directional gradients
if (isempty(I))     
    % Gx, Gy are given as inputs
    if ~isfloat(Gx)
        Gx = double(Gx);
    end
    if ~isfloat(Gy)
        Gy = double(Gy);
    end
    
else   
    % If Gx, Gy not given, compute them. For all others except Roberts
    % method, use IMGRADIENTXY to compute Gx and Gy. 
    if (strcmpi(method,'roberts'))        
        if ~isfloat(I)
            I = double(I);
        end
        Gx = imfilter(I,[1 0; 0 -1],'replicate');         
        Gy = imfilter(I,[0 1; -1 0],'replicate'); 
        
    else        
        [Gx, Gy] = imgradientxy(I,method);
        
    end
end

% Compute gradient magnitude
Gmag = hypot(Gx,Gy);

% Compute gradient direction
if (nargout > 1)
    if (strcmpi(method,'roberts'))
        
        Gdir = zeros(size(Gx));
        
        % For pixels with zero gradient (both Gx and Gy zero), Gdir is set
        % to 0. Compute direction only for pixels with non-zero gradient.
        xyNonZero = ~((Gx == 0) & (Gy == 0)); 
        Gdir(xyNonZero) = atan2(Gy(xyNonZero),-Gx(xyNonZero)) - (pi/4);
        Gdir(Gdir < -pi) = Gdir(Gdir < -pi) + 2*pi; % To account for the discontinuity at +-pi.
        
        Gdir = Gdir*180/pi; % Radians to degrees
        
    else
        
        Gdir = atan2(-Gy,Gx)*180/pi; % Radians to degrees
    end    
end

end
%======================================================================
function [I, Gx, Gy, method] = parse_inputs(varargin)

methodstrings = {'sobel','prewitt','roberts','centraldifference', ...
            'intermediatedifference'};
I = []; 
Gx = []; 
Gy = [];
method = 'sobel'; % Default method

if (nargin == 1)
    I = varargin{1};
    validateattributes(I,{'numeric','logical'},{'2d','nonsparse','real'}, ...
                       mfilename,'I',1);
        
else % (nargin == 2)
    if ischar(varargin{2})
        I = varargin{1};
        validateattributes(I,{'numeric','logical'},{'2d','nonsparse', ...
                           'real'},mfilename,'I',1);
        method = validatestring(varargin{2}, methodstrings, ...
            mfilename, 'Method', 2);
    else
        Gx = varargin{1};
        Gy = varargin{2}; 
        validateattributes(Gx,{'numeric','logical'},{'2d','nonsparse', ...
                           'real'},mfilename,'Gx',1);
        validateattributes(Gy,{'numeric','logical'},{'2d','nonsparse', ...
                           'real'},mfilename,'Gy',2);
        if (~isequal(size(Gx),size(Gy)))
            error(message('images:validate:unequalSizeMatrices','Gx','Gy'));
        end
    end
         
end

end
%----------------------------------------------------------------------
