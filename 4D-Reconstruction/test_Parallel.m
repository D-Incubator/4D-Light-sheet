clear;

%% Set before running the code
% Due to limiation of file size, the example data has 30 2D image slices. 
% The full dataset including 130 slices is available from the corresponding 
% author upon request.

% set data folder 
baseDir = 'C:\Users\Austin\Desktop\DLoad\4D-example-data\4D-example-data\cmlc2-nucGFP-7dpf';
dirNames = split(baseDir,"\"); 
dataName = char(dirNames(end));
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
outputDir = [baseDir '\' subDir]; %disp(numOfSlice);

% get height and width of images
imageList = dir([baseDir '\' int2str(1) '\*.tif']);
tempIMG = imread([baseDir '\' int2str(1) '\' imageList(1).name]); % read in one image
[imgHeight, imgWidth] = size(tempIMG);
clear tempIMG
clear imageList

%open parallel pool
parfor i=1
end
t0 = tic();
t_period = zeros(1,numOfSlice); % timer
t_relativeshift = zeros(1,numOfSlice*2-2); % timer

%% Get period (with parallel computation)
disp('Getting period...')

t_start1 = toc(t0);
periodTh1 = (systolicPoint_4st-systolicPoint_1st)/3*h_T*0.85 - 30;
periodTh2 = (systolicPoint_4st-systolicPoint_1st)/3*h_T*1.15 + 30;
t_p_candidate = zeros(1,numOfSlice);

% *can't use global variables in parfor loop
parfor i = 1:numOfSlice
    imageList = dir([baseDir '\' int2str(i) '\*.tif']); % obtain image list    

    images= uint16.empty( imgHeight , imgWidth , 0 ); % create matrix to store all images

    for j = 1:numOfImage
        images(:,:,j) = imread([baseDir '\' int2str(i) '\' imageList(j).name]); % read all images in to the matrix
    end
    
    t_p_candidate(i) = getPeriod(periodTh1,periodTh2,images,h_T,numOfImage); % get period from input data

    t_period(i) = toc(t0)-t_start1;
    disp(['====> Iteration ' int2str(i) ' in period loop is complete...']);
end
t_p = sum(t_p_candidate)/length(t_p_candidate); % get the average estimated heartbeat period of all slices
%disp(t_p); 
t_period_all = max(t_period); %timer
t_period_all2 = toc(t0) - t_start1;
disp(['Getting period took ' num2str(t_period_all2) ' seconds overall...']);

%% Get relative shift (with parallel computation)
disp('Getting relative shift...')

t_start2 = toc(t0);
rserial = 1; %timer
maxSliceConsidered=2;
numOfSliceMinus1 = numOfSlice-1;
Q = -500.*ones(numOfSlice, numOfSlice); % creat a 54*54 matrix with all -500 value
startIndex = floor(t_p/h_T)+1;

parfor i = 1:numOfSlice-1
    % first get images1, instead of getting it every time you get images2
    imageList = dir([baseDir '\' int2str(i) '\*.tif']); % obtain image list    
    images1= int16.empty( imgHeight , imgWidth , 0 ); % create matrix to store all images
        
    for k = 1:numOfImage 
        images1(:,:,k) = imread([baseDir '\' int2str(i) '\' imageList(k).name]); % read all images in to the matrix
    end

    images1 = gpuArray(images1);

    for j = 1:numOfSliceMinus1
        if( j < i || j > i+2)
            continue;
        end
        if i == j % diagonals are 0
            Q(i,j) = 0;            
        elseif j > i % anti-symmetric matrix
            imageList = dir([baseDir '\' int2str(j) '\*.tif']); % obtain image list    
            images2= int16.empty( imgHeight , imgWidth , 0 ); % create matrix to store all images
            
            for k = 1:numOfImage 
                images2(:,:,k) = imread([baseDir '\' int2str(j) '\' imageList(k).name]); % read all images in to the matrix
            end

            images2 = gpuArray(images2)
            
            cost = inf;
            for s = -floor((startIndex-1)/2):floor(startIndex/2)   
                energy = 0;
                for ind = startIndex:3*startIndex-1
                    if s+ind <= numOfImage
                        buffer = images1(:,:,ind) - images2(:,:,s+ind); %subtraction
                        buffer = buffer.^2;                       
                        energy = energy + sum(sum(buffer));
                        % if energy is already too high, no need to keep
                        % calculating
                        if energy >= cost
                            break
                        end
                    end
                end
                if energy < cost           
                    Q(i,j) = s;
                    cost = energy;
                end
            end

        end        
    end
    disp(['====> Iteration ' int2str(i) ' in relative shift loop is complete...']);
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
disp(['Getting relative shift took ' num2str(t_relativeshift_all) ' seconds overall...']);

%% Get absolute shift
disp('Getting absolute shift...')

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
disp(['Getting absolute shift took ' num2str(t_absoluteshift_all) ' seconds overall...']);

%% Final aligning and resample
disp('Resampling and writing output...')

% create the necessary directories
makeOutputDirs(outputDir);    
byStateDir = [outputDir '\byState'];

t_start4 = toc(t0);
numOfImage = round(numOfPeriod * t_p/h_T);

% resample and write output by slice
for i = 1:numOfSlice
    imageList = dir([baseDir '\' int2str(i) '\*.tif']);
    images= zeros( imgHeight , imgWidth , floor(numOfPeriod*t_p/h_T)+2 , 'uint16' );       
    count = 1;
    % align all slice
    for j = t(i)+1 : t(i)+1+floor(numOfPeriod*t_p/h_T)
        images(:,:,count) = imread([baseDir '\' int2str(i) '\' imageList(j).name]);            
        count = count + 1;
    end     

    % save resample data(temporarily)         
    images_resampled = getResample(images, t_p , numOfPeriod , numOfImage, h_T );
    for j = 1:numOfImage
        combine_name = [byStateDir '\' 'state_' num2str(j,'%04d') '.tif'];

        if (i == 1)
            imwrite(uint16(images_resampled(:,:,j)),combine_name);            
        else   
            imwrite(uint16(images_resampled(:,:,j)),combine_name, 'WriteMode', 'append');
        end
    end         
end

t_resample = toc(t0)-t_start4; %timer
disp(['Resampling and writing took ' num2str(t_resample) ' seconds overall...']);

%% Record used time
t_align_all = t_relativeshift_all + t_absoluteshift_all;
t_all = t_period_all + t_align_all + t_resample; 
timerParallel(outputDir,t_period_all,t_align_all, t_resample,t_all,t_period,t_relativeshift,t_relativeshift_all,t_absoluteshift_all);
disp(['Total time elapsed: ' num2str(toc(t0)) ' seconds...'])

%% Record average period(t_p); potial period(t_p_candidate); registration result(t)
dataRecorder(outputDir,t_p,t_p_candidate,t);
