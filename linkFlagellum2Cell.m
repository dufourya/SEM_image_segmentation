function flg = linkFlagellum2Cell(flg,coord,cellId)

%%
flg.cellId = unique([flg.cellId, cellId]);

% idx = sum(flg.endpointsCoordinates-coord==0,2)==2;

endpoints = bwmorph(flg.maskSkeleton,'endpoints');
coord = coord - flg.boundingBox(1:2);

D = bwdistgeodesic(flg.maskSkeleton,coord(2),coord(1),'quasi-euclidean');
D = round(D*8)/8;

endpoints = D==max(D(endpoints)) & endpoints;
[ex, ey] = find(endpoints);
ex = [coord(1) ex(1)];
ey = [coord(2) ey(1)];

flg.longestEndpointsCoordinates = [ex, ey];

D1 = bwdistgeodesic(flg.maskSkeleton, ey(1), ex(1), 'quasi-euclidean');
D2 = bwdistgeodesic(flg.maskSkeleton, ey(end), ex(end), 'quasi-euclidean');

D = D1 + D2;
D = round(D * 8) / 8;

D(isnan(D)) = inf;
paths = imregionalmin(D);

solution_path = bwmorph(paths, 'thin', inf);
flg.longestPath = min(D(solution_path));
flg.meanIntensityLongest = mean(flg.image(solution_path));
flg.maskLongestPath = solution_path;

%     imshow(labeloverlay(uint8(imflg),flgmask0+flgmask+solution_path));
%     hold on;
%     plot(ey, ex,'w*');
%     plot(by, bx,'wo');
%     hold off;
%     pause;

end
