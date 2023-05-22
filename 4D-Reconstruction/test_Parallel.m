% Instruction: This code serves as the 4D reconstruction algorithm for
% zebrafish heart. The idea of this method originates from the paper: X.
% Zhang et al., "4D Light-Sheet Imaging and Interactive Analysis of 
% Cardiac Contractility in Zebrafish Larvae".It combines parallel computing
% with retrospective synchronization for the four-dimensional (4D) 
% reconstruction of quasi-periodic cardiac dynamics.
%
%
% Disclaimer: This code was initiated by Xinyuan Zhang on 
% 7/5/2022. For any information, please reach out to Xinyuan at 
% xinyuan.zhang@utdallas.edu or Yichen at yichen.ding@utdallas.edu.
%
% Copyright (c) 2023, Xinyuan Zhang, Yichen Ding
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the
%       distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



clear;

%% Set before running the code
% Due to limiation of file size, the example data has 30 2D image slices. 
% The full dataset including 130 slices is available from the corresponding 
% author upon request.

% set data folder 
baseDir = 'example\cmlc2-nucGFP-7dpf';
dirNames = split(baseDir,"\"); dataName = char(dirNames(end));
% set total number of 2D image slices
numOfSlice = 30; 
% set number of images in each 2D image slice
numOfImage = 300;
% set exposure time
h_T = 5;   
% set number of periods you want to reconstruct
numOfPeriod = 2; 
% check any 2D image slice,find frame # of first and fourth systole
systolicPoint_1st = 40; 
systolicPoint_4st = 280; 
%-------------------------------------------------------------------------%
subDir = ['output\' num2str(numOfSlice) 'slices' num2str(numOfImage) 'images_' dataName]; disp(subDir);
outputDir = [baseDir '\' subDir]; disp(numOfSlice);
%open parallel pool
parfor i=1
end
t0 = tic();
t_period = zeros(1,numOfSlice); % timer
t_relativeshift = zeros(1,numOfSlice*2-2); % timer

%% Get period (with parallel computation)
t_start1 = toc(t0);
periodTh1 = (systolicPoint_4st-systolicPoint_1st)/3*h_T*0.85 - 30;
periodTh2 = (systolicPoint_4st-systolicPoint_1st)/3*h_T*1.15 + 30;
t_p_candidate = zeros(1,numOfSlice);
% *can't use global variables in parfor loop
parfor i = 1:numOfSlice
    t_p_candidate(i) = getPeriod_wrapper([baseDir '\' int2str(i)], h_T, periodTh1, periodTh2,numOfImage); %get the estimated heartbeat period in each 2D image slice
    t_period(i) = toc(t0)-t_start1;
    disp(i);
end
t_p = sum(t_p_candidate)/length(t_p_candidate); % get the average estimated heartbeat period of all slices
disp(t_p); 
t_period_all = t_period(numOfSlice); %timer
disp("Got Period");disp(t_period_all);

%% Get relative shift (with parallel computation)
t_start2 = toc(t0);
rserial = 1; %timer
maxSliceConsidered=2;
numOfSliceMinus1 = numOfSlice-1;
Q = -500.*ones(numOfSlice, numOfSlice); % creat a 54*54 matrix with all -500 value

parfor i = 1:numOfSlice-1
    for j = 1:numOfSliceMinus1
        if( j < i || j > i+2)
            continue;
        end
        if i == j % diagonals are 0
            Q(i,j) = 0;            
        elseif j > i % anti-symmetric matrix
            Q(i,j) = getRelativeShift_wrapper( baseDir, i, j, t_p, h_T, numOfImage );       
%             t_relativeshift(rserial) = toc(t0)-t_start2; rserial = rserial+1;            
        end        
    end
    fprintf('%d is complete\n',i);
    disp(toc(t0)-t_start2);
end

for i = 1:numOfSlice-1
    for j = 1:numOfSlice-1
        if i < j
            Q(j,i) = -Q(i,j);
        end
    end
end
t_relativeshift(rserial) = toc(t0)-t_start2;
t_relativeshift_all = t_relativeshift(rserial); %timer
disp("Got Relative Shift");disp(t_relativeshift_all);

%% Get absolute shift
t_start3 = toc(t0);
% maxSliceConsidered=2;
maxTermInS = (numOfSlice)*(numOfSlice-1)/2 - (numOfSlice-maxSliceConsidered)*(numOfSlice-maxSliceConsidered-1)/2;
A = zeros(maxTermInS, numOfSlice);
s = zeros(maxTermInS,1);
A(1,1) = 1;
W = zeros(maxTermInS,maxTermInS);
W(1,1) = 1;
k = 2;  
weight = 0.8;
for i = 1:maxSliceConsidered
   for j = 1:numOfSlice
      if( j+i > numOfSlice )
          continue;
      end
      s(k) = Q(j,j+i); 
      A(k,j) = 1;
      A(k,j+i) = -1;
      W(k,k) = weight;
      k = k+1;
   end
   weight = weight * 0.8;
   if i==maxSliceConsidered-1
       weight = 0;
   end
end
t = -inv(A'*W'*W*A) * A'*W'*W*s; % t is the absolute shift calculated
t(1) = 0;
t = mod(t,t_p/h_T);
t = floor(t);
t_absoluteshift_all = toc(t0)-t_start3; %timer
disp("Got Abosolute Shift");disp(t_absoluteshift_all);

%% Final aligning and resample
t_start4 = toc(t0);
numOfImage = round(numOfPeriod * t_p/h_T);
regulizeData_wrapper(baseDir, outputDir, t_p, t, h_T, numOfSlice, numOfPeriod, numOfImage); %resample data and save bySlice output
t_finalize = toc(t0)-t_start4; %timer
disp("Got Resampled");disp(t_finalize);
output_wrapper(outputDir, numOfSlice, numOfImage); %save byStaet output as 3D tiff files
rmdir(append(outputDir,'\bySlice'), 's') %*delete bySlice output

%% Record used time
t_align_all = t_relativeshift_all + t_absoluteshift_all;
t_all = t_period_all + t_align_all; 
timerParallel(outputDir,t_period_all,t_align_all,t_all,t_period,t_relativeshift,t_relativeshift_all,t_absoluteshift_all);

%% Record average period(t_p); potial period(t_p_candidate); registration result(t)
dataRecorder(outputDir,t_p,t_p_candidate,t);
