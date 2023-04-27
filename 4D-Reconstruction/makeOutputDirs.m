function makeOutputDirs(outputDir)
    if exist(outputDir,'dir') == 7     
        disp('output folder exist, result will be overwritten, press any key to continue, ctrl+c to abort');
        pause;
        rmdir(outputDir, 's');
    end
    mkdir(outputDir); 

    byStateDir = [outputDir '\byState'];
    if exist(byStateDir,'dir') == 7
        rmdir(byStateDir, 's');
    end
    mkdir(byStateDir);  