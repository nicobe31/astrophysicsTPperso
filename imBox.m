function [subImage] = imBox(name)
%%
% [subregion] = imBox(name)
%
% This function determine a sub-image around the spot light.
% 
% Argument: -name: string containing the name or the path to the image
%            name= 'im/IMG_0001.CR2'
% Output: Matrix corresponding to the sub-image
%%

close all

displayProgressiveBox = false;

[xx,map,alpha] = imread(name);
X = xx(:,:,1);
Xsize = size(X);

% noise, standard deviation and maximum of the image
noiseBI = mean2(X);
sBI = std2(X);
[maxX,iX] = max(X,[],2);
[maxInt,iMaxY] = max(maxX,[],1);
iMaxX = iX(iMaxY);

% threshold for the spot detection and factor for the box enlargement
threshold = double(noiseBI + maxInt/3);
thresholdMean = double(noiseBI + sBI*2.5);
factorInd = 1.7;
limiteWhile = 20;

%% determine box around the spot
indh1 = find(X(iMaxY,1:iMaxX)>threshold,1,'first');
indh2 = find(X(iMaxY,iMaxX:end)>threshold,1,'last') + iMaxX;
indv1 = find(X(1:iMaxY,iMaxX)>threshold,1,'first');
indv2 = find(X(iMaxY:end,iMaxX)>threshold,1,'last') + iMaxY;

if displayProgressiveBox
    squareX = [indh2,indh1,indh1,indh2,indh2];
    squareY = [indv1,indv1,indv2,indv2,indv1];
    figBoxConstru = figure('Name','BoxConstruction');
    hold on
    image(X);
    plot(squareX,squareY,'-o', 'MarkerEdgeColor','r','MarkerSize',5)
    axis ij 
% hold off
end


% optimization
% North side = X(indv1,indh1:indh2);
% West side = X(indv1:indv2,indh1);
% South side = X(indv2,indh1:indh2);
% East side = X(indv1:indv2,indh2);
test = 4;
iter = 0;
while test > 0 && iter<limiteWhile
    test = 4;
    
    % Norht side:
    meanN = mean(X(indv1,indh1:indh2));
    if meanN >= thresholdMean
        [~,indN] = max(X(indv1,indh1:indh2));
        indN = indN + indh1 -1;
        indv1Prev = indv1;
        indv1 = find(X(1:indv1,indN)<threshold,1,'last');
        if indv1 <= indv1Prev -1
            test = test - 1;
        end
    else 
        test = test - 1;
    end

    
    % West side
    meanW = mean(X(indv1:indv2,indh1));
    if meanW >= thresholdMean
        [~,indW] = max(X(indv1:indv2,indh1));
        indW = indW + indv1-1;
        indh1Prev = indh1;
        indh1 = find(X(indW,1:indh1)<threshold,1,'last');
        if indh1 <= indh1Prev -1
            test = test -1;
        end
    else 
        test = test - 1;
    end
    
    % South side
    meanS = mean(X(indv2,indh1:indh2));
    if meanS >= thresholdMean
        [~,indS] = max(X(indv2,indh1:indh2));
        indS = indS + indh1-1;
        indv2Prev = indv2;
        indv2 = find(X(indv2:end,indS)<threshold,1,'first') + indv2;
        if indv2 <= indv2Prev +1
            test = test -1;
        end
    else 
        test = test - 1;
    end
    
    % East side
    meanE = mean(X(indv1:indv2,indh2));
    if meanE >= thresholdMean 
        [~,indE] = max(X(indv1:indv2,indh2));
        indE = indE + indv1-1;
        indh2Prev = indh2;
        indh2 = find(X(indE,indh2:end)<threshold,1,'first') + indh2;
        if indh2 <= indh2Prev +1
            test = test -1;
        end
    else 
        test = test - 1;
    end
    
    if displayProgressiveBox
        squareX = [indh2,indh1,indh1,indh2,indh2];
        squareY = [indv1,indv1,indv2,indv2,indv1];
        figure(1)
        hold on
        plot(squareX,squareY,'-o', 'MarkerEdgeColor','r','MarkerSize',5)
        axis ij 
        pause  
    end
    iter = iter+1;
end

% enlargement
indh1 = indh1 + floor((indh2-indh1)/2 - factorInd*(indh2-indh1)/2);
indh2 = indh1 + floor((indh2-indh1)/2 + factorInd*(indh2-indh1)/2);
indv1 = indv1 + floor((indv2-indv1)/2 - factorInd*(indv2-indv1)/2);
indv2 = indv1 + floor((indv2-indv1)/2 + factorInd*(indv2-indv1)/2);


%% plot
squareX = [indh2,indh1,indh1,indh2,indh2];
squareY = [indv1,indv1,indv2,indv2,indv1];
figBox = figure('Name','ImageBox');
hold on
image(X);
plot(iMaxX,iMaxY,'or', 'MarkerSize', 5);
plot(squareX,squareY,'-o', 'MarkerEdgeColor','c','MarkerSize',7)
% axis ij
legend('Max','Box');
hold off

%% 

promp = sprintf('Is this box ok? y/n  ');
ok = input(promp,'s');
while ~strcmpi(ok,'y')
    
    disp('Use the mouse to define up-left and bottom-right corners')
    
    figure(figBox)
    % Use ginput to select corner points of a rectangular
    % region by pointing and clicking the mouse twice
    p = ginput(2); 

    % Get the x and y corner coordinates as integers
    indh1 = min(floor(p(1)), floor(p(2))); %xmin
    indv1 = min(floor(p(3)), floor(p(4))); %ymin
    indh2 = max(ceil(p(1)), ceil(p(2)));   %xmax
    indv2 = max(ceil(p(3)), ceil(p(4)));   %ymax

    % Index into the original image to create the new image
    MM = X(indv1:indv2, indh1:indh2,:);

    % Display the subsetted image with appropriate axis ratio
    if ~exist('figOwnBox','var')
        figOwnBox = figure('Name','OwnBoxConstruction');
    end
    figure(figOwnBox)
    MM2 = imread(name,'PixelRegion',{[1 indv2], [1 indh2]});
    image(MM2);
    
    ok = input(promp,'s');
    
    %% 
    % Mathworks, _Read, Write, and Query Image Files_, _Subsetting a Graphics Image (Cropping)_,
    % <https://nl.mathworks.com/help/matlab/creating_plots/reading-writing-and-querying-graphics-image-files.html#f2-19897>,
    % visited on (11/12/2016)
    %
end

%% Subregion definition

% subregion = {[indv1,indv2], [indh1,indh2]};
subImage = X(indv1:indv2, indh1:indh2,:);

end


% [mInt, linInd] = max(X(:));
% indX = ceil(linInd/Xsize(1));
% indY = mod(linInd,Xsize(1));

