function A = angle(subImage)
%%
% A = angle(name,subImage)
%
% This function compute the angle in degree of the fringes with respect to
% the horizontal.
% 
% Argument: subImage: the sub-image obtain with imBox (not the limite of
%                     the box but the matrix).
% Output: A: the angle of the fringes with respect to the horizon [°]
%%

close all
num = 47;
name = num2str(num,'im/IMG_%.4d.CR2');

[subImage] = imBox(name);

image(subImage)

sizeImage = size(subImage);

threshold = max(subImage(:))*0.8;

[pks1,locs1] = findpeaks(double(subImage(:)),'MinPeakHeight',threshold);
ssubImage = subImage';
[pks2,locs2] = findpeaks(double(ssubImage(:)),'MinPeakHeight',threshold);
    
y1 = mod(locs1,sizeImage(1));
x1 = ceil(locs1/sizeImage(1));
x2 = mod(locs2,sizeImage(2));
y2 = ceil(locs2/sizeImage(2));
x = [x1;x2];
y = [y1;y2];
pks = [pks1; pks2];
nPks = length(x);

figure
hold on
surf(double(subImage));
plot3(x1,y1,pks1,'or')
plot3(x2,y2,pks2,'og')
plot3(x,y,pks,'*k')
hold off

figConstru = figure;
hold on
surf(double(subImage));
plot3(x,y,pks,'*k')

% dx = zeros(nPks,nPks-1);
% dy = zeros(nPks,nPks-1);
% 
% for i=2:nPks
%     xx = [x(i:end);x(1:i-1)];
%     yy = [y(i:end);y(1:i-1)];
%     
%     dx(:,i-1) = x-xx;
%     dy(:,i-1) = y-yy;
% end
% 
% ds = dx.^2 + dy.^ 2;
radius = 3^2;

pksInd = [1:nPks]';
pksdoneInd = 1;
pksdone = [];
grpNbre = 1;

while length(pksdone)<nPks 
    grp=0;
    grpInd = 1;
    pksList = pksInd;
    pksList(pksdone) = [];
    indSearch = pksList(1);
    pksdone(pksdoneInd) = indSearch;
    grp(grpInd) = indSearch;
    grpInd = grpInd+1;
    pksdoneInd = pksdoneInd+1;
    
    j=1;
    while ~isempty(j)
        dx = x - x(indSearch);
        dy = y - y(indSearch);
        ds = dx.^2 + dy.^2;
        ds([pksdone;indSearch]) = [];
        pksSearch = pksInd;
        pksSearch([pksdone;indSearch]) = [];

        [j] = find(ds<radius);

        nbrePt = length(j);
        grp(grpInd:grpInd+nbrePt-1,1) = pksSearch(j);

        figure(figConstru);
        plot3(x(indSearch),y(indSearch),pks(indSearch),'gx')
        plot3(x(grp),y(grp),pks(grp),'or');
    %     pause

    %     [~, imax] = max(ds(indSearch,pksInd(j)));
    %     indSearch = pksInd(imax);
        pksdone(pksdoneInd:pksdoneInd+nbrePt-1,1) = pksSearch(j);
        indSearch = grp(end,1);
        grpInd = grpInd+nbrePt;
        pksdoneInd = pksdoneInd+nbrePt;

    %     pksInd(j) = [];  
    end
    Grp{grpNbre} = grp;
    grpNbre = grpNbre +1;
end

sym = ['or','og','ok','ob','oc','om'];
figure
hold on
surf(double(subImage));
for i = 1:grpNbre-1
    plot3(x(Grp{i}),y(Grp{i}),pks(Grp{i}),sym(mod(i-1,6)+1));
end
hold off

ind = 1;
for i = 1:grpNbre-1
    if length(Grp{i})>5
        mdx(ind) = sum(diff(x(Grp{i})));
        mdy(ind) = sum(diff(y(Grp{i})));
        ind = ind+1;
    end
end
mA = tand(mdy./mdx);
A = mean(mA);

lineLength = 2;
% lineX = [floor(sizeImage(2)/2 - lineLength/tand(A)),floor(sizeImage(2)/2 + lineLength/tand(A))];
% lineY = [floor(sizeImage(1)/2 - lineLength*tand(A)),floor(sizeImage(1)/2 + lineLength*tand(A))];
lineX = [floor(sizeImage(2)/2 - lineLength*mean(mdx)),floor(sizeImage(2)/2 + lineLength*mean(mdx))];
lineY = [floor(sizeImage(1)/2 - lineLength*mean(mdy)),floor(sizeImage(1)/2 + lineLength*mean(mdy))];

figure
hold on
surf(double(subImage));
plot3(lineX,lineY,max(pks)*ones(2,1),'o-r');
hold off

end

% figure
% hold on
% surf(double(subImage));
% plot3(x(i),y(j),300*ones(length(i),1),'or')
% hold off

