function [Q_i_j] = getRelativeShift_wrapper( baseDir, index1, index2, t_p, h_T,numOfImage)

        
%%  read images of slice i
    imageList = dir([baseDir '\' int2str(index1) '\*.tif']); % obtain image list    
%     numOfImage1 = length(imageList); % obtain num of images
    numOfImage1 = numOfImage;
    
    tempIMG = imread([baseDir '\' int2str(index1) '\' imageList(1).name]); % read in one image
    [height, width] = size(tempIMG); % get the height and length of the images
    
    startIndex = floor(t_p/h_T)+1;
    
    %%%%%%%%%%%%%%%%###########################!!!!!!!!!!!!!
%     startIndex = floor(numOfImage/3);
    
    images1= zeros( height , width , numOfImage1 ); % create matrix to store all images
        
    for i = 1:numOfImage1
        
        images1(:,:,i) = imread([baseDir '\' int2str(index1) '\' imageList(i).name]); % read all images in to the matrix
        images1 = gpuArray(images1);
    end
    
 %%  read images of slice j
    imageList = dir([baseDir '\' int2str(index2) '\*.tif']); % obtain image list    
%     numOfImage2 = length(imageList); % obtain num of images   
    numOfImage2 = numOfImage;


    images2= zeros( height , width , numOfImage2 ); % create matrix to store all images
    
    for i = 1:numOfImage2
        
        images2(:,:,i) = imread([baseDir '\' int2str(index2) '\' imageList(i).name]); % read all images in to the matrix
        images2 = gpuArray(images2);
    end
    
%%  calculate relative shift between slice i and j
    Q_i_j = getRelativeShift(startIndex,images1,images2);