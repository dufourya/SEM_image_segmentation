function runSEMSingle(folder)

folder_short = strsplit(folder,filesep);
folder_short = folder_short{end};
f = dir(strcat(folder,filesep,'*.tif'));
marker = cell(numel(f),1);
bar = cell(numel(f),1);
mag = cell(numel(f),1);
imsize = cell(numel(f),1);
for k = 1:numel(f)
    % im = imread(f(k).name);
    fid = fopen(fullfile(folder,strrep(f(k).name,'.tif','.txt')));
    meta = textscan(fid,'%s %[^\n\r]');
    fclose(fid);
    marker{k} = meta{2}{strcmp(meta{1},'$$SM_MICRON_MARKER')};
    bar{k} = meta{2}{strcmp(meta{1},'$$SM_MICRON_BAR')};
    mag{k} = meta{2}{strcmp(meta{1},'$CM_MAG')};
    imsize{k} = meta{2}{strcmp(meta{1},'$CM_FULL_SIZE')};
end
%%
% unique(marker)
% unique(bar)
% unique(mag)
% unique(imsize)
%%
fprintf('Reading images...');
ind = (strcmp(mag,'1600'));
f = f(ind);
im = zeros(1920,2560,numel(f));
for k =1:numel(f)
    tmp = imread(fullfile(folder,f(k).name));
    ims = sscanf(imsize{k},'%d');
    im(:,:,k) = tmp(1:ims(2),1:ims(1));
end
fprintf('done.\n');
%%
fprintf('Analyzing images...000/000');

flgAll = cell(size(im,3),1);
cellAll = cell(size(im,3),1);
maskCells = zeros(size(im));
maskFlagella = zeros(size(im));
n = size(im,3);

for k = 1:n
    
    fprintf('\b\b\b\b\b\b\b%0.3d/%0.3d',k,n);
    imk = im(:,:,k);
    maskCell = SEMcellMask(imk,4,0.05,400);
    maskFlg = SEMflagellaMask(imk,bwmorph(maskCell,'thicken',1),1,0.05,100);
    
    %     cell_props = regionprops(maskCell,'EquivDiameter',...
    %         'MinoraxisLength','Eccentricity','EulerNumber');
    
    %     good_cells = find([cell_props.EulerNumber] == 1 & ...
    %         [cell_props.EquivDiameter] > 25 & [cell_props.EquivDiameter] < 55 & ...
    %         [cell_props.MinorAxisLength] > 15 & [cell_props.MinorAxisLength] < 40 & ...
    %         [cell_props.Eccentricity] > 0.9);
    
    %     maskCell(~ismember(maskCell(:),good_cells)) = 0;
    %     imshow(labeloverlay(uint8(imk),maskCell),[],'InitialMagnification',100);
    
    %     flg_props = regionprops(maskFlg, ...
    %         'MajoraxisLength');
    
    %     good_flg = find([flg_props.MajorAxisLength] > 10);
    
    %     maskFlg(~ismember(maskFlg(:),good_flg)) = 0;
    %     imshow(labeloverlay(uint8(imk),maskFlg),[],'InitialMagnification',100);
    
    maskCells(:,:,k) = maskCell;
    maskFlagella(:,:,k) = maskFlg;
    
    %     pause;
    
    flagella = [];
    idx = unique(maskFlg(:))';
    idx = idx(idx>0);
    for j = idx
        flgmask = maskFlg == j;
        y = find(sum(flgmask,1));
        x = find(sum(flgmask,2));
        bb = [max(1,min(x)-10) max(1,min(y)-10) min(max(x)+10,size(im(:,:,k),1)) min(max(y)+10,size(im(:,:,k),2))];
        imflg = im(bb(1):bb(3),bb(2):bb(4),k);
        flgmask = flgmask(bb(1):bb(3),bb(2):bb(4));
        flgd = flagellum2obj(imflg,flgmask);
        flgd.id = j;
        flgd.file = f(k).name;
        flgd.dir = folder_short;
        flgd.boundingBox = bb;
        flagella = [flagella flgd];
    end
    flgAll{k} = flagella;
    
    cells = [];
    idx = unique(maskCell(:))';
    idx = idx(idx>0);
    for j = idx
        cellmask = maskCell == j;
        y = find(sum(cellmask,1));
        x = find(sum(cellmask,2));
        bb = [max(1,min(x)-10) max(1,min(y)-10) min(max(x)+10,size(im(:,:,k),1)) min(max(y)+10,size(im(:,:,k),2))];
        imcell = im(bb(1):bb(3),bb(2):bb(4),k);
        cellmask = cellmask(bb(1):bb(3),bb(2):bb(4));
        celld = cell2mesh(imcell,cellmask);
        celld.id = j;
        celld.file = f(k).name;
        celld.dir = folder_short;
        celld.boundingBox = bb;
        cells = [cells celld];
    end
    cellAll{k} = cells;
    
end

fprintf('...done.\n');

%%
fprintf('Writing masked images...');
for k = 1:size(im,3)
    imk = im(:,:,k);
    imk = imk-min(imk(:));
    imk = imk/max(imk(:));
    imc = labeloverlay(imk,(maskCells(:,:,k)>0)+3*(maskFlagella(:,:,k)>0));
    %     imshow(imc);
    imwrite(imc,fullfile(folder,strrep(f(k).name,'.tif','.png')));
end
fprintf('done.\n');
%%
% fprintf('Filtering flagella and cells...');
% 
% for k = 1:numel(cellAll)
%     
%     flagella = flgAll{k};
%     ind = [flagella.longestPath] > 10;
%     flgAll{k} = flagella(ind);
%     
%     cells = cellAll{k};
%     ind = [cells.meanIntensity] - [cells.modeBackgroundIntensity] >0 & sqrt([cells.Area]) > 25;
%     cellAll{k} = cells(ind);
% end
% 
% fprintf('done.\n');

%%
fprintf('Connecting flagella to cells...');

for k = 1:numel(cellAll)
    flagella = flgAll{k};
    flgcoord = [];
    for j = 1:numel(flagella)
        flagella(j).id = j;
        coord = flagella(j).endpointsCoordinates;
        if ~isempty(coord)
            bb = flagella(j).boundingBox;
            coord(:,1) = coord(:,1) + bb(1);
            coord(:,2) = coord(:,2) + bb(2);
            coord = cat(2,coord,flagella(j).id*ones(size(coord,1),1));
            flgcoord = [flgcoord; coord];
        end
    end
    
    cells = cellAll{k};
    polcoord = [];
    for j = 1:numel(cells)
        cells(j).id = j;
        coord = cells(j).meshPoleCoordinates;
        if ~isempty(coord)
            bb = cells(j).boundingBox;
            coord(:,1) = coord(:,1) + bb(1);
            coord(:,2) = coord(:,2) + bb(2);
            coord = cat(2,coord,cells(j).id*ones(size(coord,1),1));
            polcoord = [polcoord; coord];
        end
    end
    
    [indCell,dCell] = dsearchn(polcoord(:,1:2),flgcoord(:,1:2));
    [indFlg,~] = dsearchn(flgcoord(:,1:2),polcoord(:,1:2));
    
    ii = find(indFlg(indCell)' == 1:numel(indCell));
    ii = ii(dCell(ii)<10);
    %     jj = find(indCell(indFlg)' == 1:numel(indFlg));
    jj = indCell(ii);
    
    %     jj = jj(dFlg(jj)<10);
    %     ii = ii(dCell(ii)<10);
    
    for j = 1:numel(ii)
        idx = flgcoord(ii(j),3);
        flagella(idx) = linkFlagellum2Cell(flagella(idx), flgcoord(ii(j),1:2), polcoord(jj(j),3));
        cells(polcoord(jj(j),3)).flagellumId = unique([cells(polcoord(jj(j),3)).flagellumId idx]);
        
    end
    
    flgAll{k}= flagella;
    cellAll{k} = cells;
    
    %     indCells = [flagella.cellId];
    %     indFlg = [cells.flagellumId];
    
    %     imshow(labeloverlay(uint8(im(:,:,k)),...
    %         (4*(maskCells(:,:,k)>0) + 3*(maskFlagella(:,:,k)>0) + 2*(ismember(maskCells(:,:,k),indCells)) + 4*(ismember(maskFlagella(:,:,k),indFlg)))));
    
    %     hold on; scatter(flgcoord(ii,2), flgcoord(ii,1), 'go');scatter(polcoord(jj,2), polcoord(jj,1), 'g*');
    %     pause;
end

flgAll = [flgAll{:}];
cellAll = [cellAll{:}];
fprintf('done.\n');

%%
fprintf('Saving analysis output...');

save(fullfile(folder,strcat(folder_short,'.SEM.mat')), 'flgAll', 'cellAll');
reduced_flgAll = struct2table(rmfield(flgAll,{'boundingBox','image','mask','maskSkeleton','branchingCoordinates','endpointsCoordinates','longestEndpointsCoordinates','maskLongestPath','numberPixelBranches'}));
reduced_cellAll = struct2table(rmfield(cellAll,{'boundingBox','image','mask','mesh','meshRibsLength','meshMedialAxis','meshPoleCoordinates','meshCurvature'}));
writetable(reduced_flgAll,fullfile(folder,strcat(folder_short,'_summary_flagella.SEM.txt')));
writetable(reduced_cellAll,fullfile(folder,strcat(folder_short,'_summary_cells.SEM.txt')));

fprintf('done.\n');
%%
fprintf('Writing filtered masked images...');

flgMaskLong = zeros(size(im));
% flgCellLong = zeros(size(im));

for k = 1:size(im,3)
    indFlg = find(contains({flgAll.file},f(k).name));
    indCell = find(contains({cellAll.file},f(k).name));
    
    tmp = flgMaskLong(:,:,k);
    
    for j = indFlg
        %         if sum(flgAll(j).maskLongestPath(:))/sum(flgAll(j).maskSkeleton(:)) >0.1
        [x,y] = find(imdilate(flgAll(j).maskLongestPath,strel('disk',2)));
        x = x + flgAll(j).boundingBox(1)-1;
        y = y + flgAll(j).boundingBox(2)-1;
        ind = sub2ind(size(tmp),x,y);
        if ~isempty(flgAll(j).cellId)
            tmp(ind) = 1;
        else
            tmp(ind) = 2;
        end
        %         end
    end
    
    for j = indCell
        %         if cellAll(j).meshLength < 200 && cellAll(j).meshRibsLengthMax < 25
        [x,y] = find(cellAll(j).mask);
        x = x + cellAll(j).boundingBox(1)-1;
        y = y + cellAll(j).boundingBox(2)-1;
        ind = sub2ind(size(tmp),x,y);
        if ~isempty(cellAll(j).flagellumId)
            tmp(ind) = 3;
        else
            tmp(ind) = 4;
        end
        %         end
    end
    flgMaskLong(:,:,k) = tmp;
end

c = lines(4);
for k = 1:size(im,3)
    imwrite(labeloverlay(uint8(im(:,:,k)),flgMaskLong(:,:,k),'Colormap',c),fullfile(folder,strrep(f(k).name,'.tif','_filtered.png')));
end

fprintf('done.\n');
end
