% Period Determination

% Parameters Desciption
%   x: horizontal pixel index
%   y: vertical pixel index
%   z: slice index
%   h_T: imaging sampling step
%   T_p: candidate Period.


% Input
%   images: all the images considered
%   h_T: imaging sampling step


% Output
%   T_p: period


function [T_p] = getPeriod(th1,th2,images,h_T,numOfImages)

%% optimization objective function
    T_p = fminbnd(@(T_p)getPeriodMatchEnergy(T_p,images,h_T,numOfImages), th1 , th2 ); % matlab build in non linear optimization function
    
    
    
