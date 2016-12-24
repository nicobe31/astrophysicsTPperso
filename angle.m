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

sizeImage = size(subImage);

displayConstruction = false;

%% detection of local peaks
threshold = max(subImage(:))*0.8;

[pks1,locs1] = findpeaks(double(subImage(:)),'MinPeakHeight',threshold);
ssubImage = subImage';
[pks2,locs2] = findpeaks(double(ssubImage(:)),'MinPeakHeight',threshold);

% transform lin indice in matrix indices
y1 = mod(locs1,sizeImage(1));
x1 = ceil(locs1/sizeImage(1));
x2 = mod(locs2,sizeImage(2));
y2 = ceil(locs2/sizeImage(2));
x = [x1;x2];
y = [y1;y2];
pks = [pks1; pks2];
nPks = length(x);

if displayConstruction
    figure
    hold on
    surf(double(subImage));
    plot3(x,y,pks,'*k')
    hold off
end

%% grouping of the fringes

if displayConstruction
    figConstru = figure;
    hold on
    surf(double(subImage));
    plot3(x,y,pks,'*k')
end

    % radius^2 used to follow the fringe 
    radius = 3^2; % [px^2]

    
pksInd = [1:nPks]';
pksdoneInd = 1;
pksdone = [];
grpNbre = 1;

while length(pksdone)<nPks % loop on the fringes
    
    % store the index of for 1 fringe
    grp=0;
    grpInd = 1;
    
    % starting pt of the fringes
    pksList = pksInd;
    pksList(pksdone) = [];
    indSearch = pksList(1);
    
    %store the pts already ine fringes
    pksdone(pksdoneInd) = indSearch;
    grp(grpInd) = indSearch;
    
    % update
    grpInd = grpInd+1;
    pksdoneInd = pksdoneInd+1;
    
    j=1;
    while ~isempty(j) % follow the fringe
        % distance to other points
        dx = x - x(indSearch);
        dy = y - y(indSearch);
        ds = dx.^2 + dy.^2;
        ds([pksdone;indSearch]) = []; %remove pt already in fringes and current pt
        
        % indices corresponding to the distances
        pksSearch = pksInd;
        pksSearch([pksdone;indSearch]) = [];
        
        % search the close pt
        [j] = find(ds<radius);
        
        % save the point 
        nbrePt = length(j);
        grp(grpInd:grpInd+nbrePt-1,1) = pksSearch(j);
        
        if displayConstruction
            figure(figConstru);
            plot3(x(indSearch),y(indSearch),pks(indSearch),'gx')
            plot3(x(grp),y(grp),pks(grp),'or');
            legend('Centeral pt','Pt in fringes')
            pause
        end

        % update
        pksdone(pksdoneInd:pksdoneInd+nbrePt-1,1) = pksSearch(j);
        indSearch = grp(end,1);
        grpInd = grpInd+nbrePt;
        pksdoneInd = pksdoneInd+nbrePt;
    end
    
    iCheck = 1;
    while iCheck < grpInd  % check after missing point
        % distance to other points
        dx = x - x(grp(iCheck));
        dy = y - y(grp(iCheck));
        ds = dx.^2 + dy.^2;
        ds(pksdone) = []; %remove pt already in fringes and current pt
        
        % indices corresponding to the distances
        pksSearch = pksInd;
        pksSearch(pksdone) = [];
        
        % search the close pt
        [j] = find(ds<radius);
        
        % save the point 
        nbrePt = length(j);
        grp(grpInd:grpInd+nbrePt-1,1) = pksSearch(j);
        
        if displayConstruction
            figure(figConstru);
            plot3(x(grp(iCheck)),y(grp(iCheck)),pks(grp(iCheck)),'mx')
            plot3(x(grp),y(grp),pks(grp),'oc');
            legend('Centeral pt','Pt in fringes')
            pause
        end

        % update
        pksdone(pksdoneInd:pksdoneInd+nbrePt-1,1) = pksSearch(j);
        grpInd = grpInd+nbrePt;
        pksdoneInd = pksdoneInd+nbrePt;
     
        iCheck = iCheck+1;
    end
    
    % save the followed fringe in a variable
    Grp{grpNbre} = grp;
    grpNbre = grpNbre +1;
end

s = sprintf('%i fringes',grpNbre-1);
disp(s);

% display of the grouped fringes
sym = ['or','og','om','oc','ok','ob'];
figure
hold on
surf(double(subImage));
for i = 1:grpNbre-1
    plot3(x(Grp{i}),y(Grp{i}),pks(Grp{i}),sym(mod(i-1,6)+1));
end
hold off


%% Angle computation
ind = 1;
for i = 1:grpNbre-1
    if length(Grp{i})>5 % don't use small fringe to compute the angle
        mdx(ind) = sum(diff(x(Grp{i})));
        mdy(ind) = sum(diff(y(Grp{i})));
        ind = ind+1;
    end
end
mA = tand(mdy./mdx);
A = mean(mA);

%% Display angle
lineLength = 2;
% lineX = [floor(sizeImage(2)/2 - lineLength/tand(A)),floor(sizeImage(2)/2 + lineLength/tand(A))];
% lineY = [floor(sizeImage(1)/2 - lineLength*tand(A)),floor(sizeImage(1)/2 + lineLength*tand(A))];
lineX = [floor(sizeImage(2)/2 - lineLength*mean(mdx)),floor(sizeImage(2)/2 + lineLength*mean(mdx))];
lineY = [floor(sizeImage(1)/2 - lineLength*mean(mdy)),floor(sizeImage(1)/2 + lineLength*mean(mdy))];

figLine = figure;
hold on
surf(double(subImage));
plot3(lineX,lineY,max(pks)*ones(2,1),'o-r');
hold off

%% Ask promp

promp = sprintf('Is the angle ok? y/n  ');
ok = input(promp,'s');
while ~strcmpi(ok,'y')
    
    disp('Use the mouse to define the begining pt and ending pt of a line aligned with the fringes')
    
    figure(figLine);
    [xp,yp] = ginput(2); 
    A = atand((yp(2)-yp(1))/(xp(2)-xp(1)));
    
    lineX = [floor(sizeImage(2)/2 - lineLength*(xp(2)-xp(1))),floor(sizeImage(2)/2 + lineLength*(xp(2)-xp(1)))];
    lineY = [floor(sizeImage(1)/2 - lineLength*(yp(2)-yp(1))),floor(sizeImage(1)/2 + lineLength*(yp(2)-yp(1)))];
    figure(figLine);
    hold on
    plot3(lineX,lineY,max(pks)*ones(2,1),'o-m');
    
    ok = input(promp,'s');
end

end

