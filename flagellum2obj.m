function flg = flagellum2obj(imflg,flgmask0)

flg = struct(...
    'file','',...
    'id',nan, ...
    'boundingBox',[],...
    'image',imflg,...
    'mask',flgmask0,...
    'maskSkeleton',[],...
    'meanIntensity',nan,...
    'modeBackgroundIntensity',nan,...
    'sumMask',[],...
    'numberEndPoints',nan,...
    'numberPixelBranches',[],...
    'branchingCoordinates',[],...
    'endpointsCoordinates',[],...
    'longestPath',nan,...
    'longestEndpointsCoordinates',[],...
    'meanIntensityLongest',[],...
    'maskLongestPath',[],...
    'sinuosity',[],...
    'cellId',[]...
    );

flg.meanIntensity = mean(imflg(flgmask0));
flg.modeBackgroundIntensity = mode(round(imflg(~flgmask0)));
flg.sumMask = sum(flgmask0(:));

flgmask = bwmorph(flgmask0,'thin', inf);
flg.maskSkeleton = flgmask;

endpoints = bwmorph(flgmask,'endpoints');
[ex, ey] = find(endpoints);
flg.endpointsCoordinates = [ex, ey];
flg.numberEndPoints = numel(ex);

if flg.numberEndPoints>1
    branchpoints = bwmorph(flgmask,'branchpoints');
    [bx, by] = find(branchpoints);
    flg.branchingCoordinates = [bx, by];
    %
    sections = flgmask;
    sections(imdilate(branchpoints,strel('disk',1,4))) = 0;
    
    branches = bwselect(sections,ey,ex,8);
    loops = bwmorph(flgmask-branches,'spur');
    branches = flgmask-loops;
    
    branches = bwlabel(branches);
    flg.numberPixelBranches = histcounts(bwlabel(branches),(0:max(branches(:)))+0.5);
    %
    D = zeros(size(flgmask));
    
    for k = 1:numel(ex)
        D = max(D,bwdistgeodesic(flgmask,ey(k),ex(k),'quasi-euclidean'));
    end
    D = round(D*8)/8;
    
    endpoints = D==max(D(endpoints)) & endpoints;
    [ex, ey] = find(endpoints,2);

    flg.longestEndpointsCoordinates = [ex, ey];
    
    D1 = bwdistgeodesic(flgmask, ey(1), ex(1), 'quasi-euclidean');
    D2 = bwdistgeodesic(flgmask, ey(end), ex(end), 'quasi-euclidean');
    
    D = D1 + D2;
    D = round(D * 8) / 8;
    
    D(isnan(D)) = Inf;
    paths = imregionalmin(D);
    
    solution_path = bwmorph(paths, 'thin', inf);
    flg.longestPath = min(D(solution_path));
    flg.meanIntensityLongest = mean(imflg(solution_path));
    flg.maskLongestPath = solution_path;
    flg.sinuosity = flg.longestPath / sqrt(diff(ex)^2+diff(ey)^2);
    
%     imshow(labeloverlay(uint8(imflg),flgmask0+flgmask+solution_path));
%     hold on;
%     plot(ey, ex,'w*');
%     plot(by, bx,'wo');
%     hold off;
%     pause;
    
end
end






