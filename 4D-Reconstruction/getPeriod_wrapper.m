function [T_p, images_resampled] = getPeriod_wrapper( path , h_T, th1, th2, numOfImage)

%%  read images
    imageList = dir([path '\*.tif']); % obtain image list    
    % numOfImage = length(imageList); % obtain num of images    
    
    tempIMG = imread([path '\' imageList(1).name]); % read in one image
    [height, width] = size(tempIMG); % get the height and length of the images     
        
    images= zeros( height , width , numOfImage ); % create matrix to store all images
    for i = 1:numOfImage
        
        images(:,:,i) = imread([path '\' imageList(i).name]); % read all images in to the matrix
        images = gpuArray(images);
    end

%% get period
    T_p = getPeriod(th1,th2,images,h_T); % get period from input data
    
    

    