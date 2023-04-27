function output_wrapper(outputDir, numOfSlice, numOfImage)  
    for i = 1:numOfImage      
        % combine image to one tiff and save
        combine_name = [byStateDir '\' 'state_' num2str(i,'%04d') '.tif'];                     
         % save reconstructed 3D images
        for j = 1:numOfSlice
            imagedata1 = imread([outputDir '\bySlice\' int2str(j) '\' int2str(i) '.tif']);
            if (j == 1)
                imwrite(imagedata1,combine_name);            
            else   
                imwrite(imagedata1,combine_name, 'WriteMode', 'append');
            end
        end          
    end