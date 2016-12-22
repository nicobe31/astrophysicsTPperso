clear all
close all

name = '40m20-8sec-5_6ouv-ISO800.CR2';

[subImage] = imBox(name);

image(subImage)

sizeImage = size(subImage);

threshold = max(subImage(:))*0.9;

[pks,locs] = findpeaks(double(subImage(:)),'MinPeakHeight',threshold);
    
yy = mod(locs,sizeImage(1));
xx = ceil(locs/sizeImage(1));

figure
hold on
surf(double(subImage));
plot3(xx,yy,pks,'or')
hold off