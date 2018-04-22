function blur = blur_metric(original, m)
% original : entry image

% The idea is from "The Blur Effect: Perception and Estimation with a New No-Reference Perceptual Blur Metric"
% Crt-Roffet F., Dolmiere T., Ladret P., Nicolas M. - GRENOBLE - 2007
% In SPIE proceedings - SPIE Electronic Imaging Symposium Conf Human Vision and Electronic Imaging, tats-Unis d'Amrique (2007)

% Written by DO Quoc Bao, PhD Student in L2TI Laboratory, Paris 13 University, France
% Email: doquocbao@gmail.com, do.quocbao@l2ti.univ-paris13.fr
% Last changed: 21-09-2008

%%%%%%%%%%%%%%%%%%Note: the output (blur) is in [0,1]; 0 means sharp, 1 means blur%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the author when you use this code. All remarks are welcome. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% https://fr.mathworks.com/matlabcentral/fileexchange/24676-image-blur-metric
%  Copyright (c) 2009, Do Quoc Bao 
%  All rights reserved.

%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:

%  * Redistributions of source code must retain the above copyright 
%  notice, this list of conditions and the following disclaimer. 
%  * Redistributions in binary form must reproduce the above copyright 
%  notice, this list of conditions and the following disclaimer in 
%  the documentation and/or other materials provided with the distribution

%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

if isempty(original), blur = inf; return; end
I = double(rgb2gray(original));
[y x] = size(I);

if nargin < 2, m=5; end

B_Ver = I;
B_Hor = I;
for index=2:m
  % blur the input image in vertical direction
  B_Ver(1:(end-m+1), :) = B_Ver(1:(end-m+1), :) + I(m:end, :);
  % blur the input image in horizontal direction
  B_Hor(:, 1:(end-m+1)) = B_Hor(:, 1:(end-m+1)) + I(:, m:end);
end
B_Ver = B_Ver/m;
B_Hor = B_Hor/m;

D_F_Ver = abs(I(:,1:x-1) - I(:,2:x));%variation of the input image (vertical direction)
D_F_Hor = abs(I(1:y-1,:) - I(2:y,:));%variation of the input image (horizontal direction)

D_B_Ver = abs(B_Ver(:,1:x-1)-B_Ver(:,2:x));%variation of the blured image (vertical direction)
D_B_Hor = abs(B_Hor(1:y-1,:)-B_Hor(2:y,:));%variation of the blured image (horizontal direction)

T_Ver = D_F_Ver - D_B_Ver;%difference between two vertical variations of 2 image (input and blured)
T_Hor = D_F_Hor - D_B_Hor;%difference between two horizontal variations of 2 image (input and blured)

V_Ver = max(0,T_Ver);
V_Hor = max(0,T_Hor);

S_D_Ver = sum(sum(D_F_Ver(2:y-1,2:x-1)));
S_D_Hor = sum(sum(D_F_Hor(2:y-1,2:x-1)));

S_V_Ver = sum(sum(V_Ver(2:y-1,2:x-1)));
S_V_Hor = sum(sum(V_Hor(2:y-1,2:x-1)));

blur_F_Ver = (S_D_Ver-S_V_Ver)/S_D_Ver;
blur_F_Hor = (S_D_Hor-S_V_Hor)/S_D_Hor;

blur = max(blur_F_Ver,blur_F_Hor);


