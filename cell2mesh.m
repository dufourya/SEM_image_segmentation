function celld = cell2mesh(imcell,cellmask)

celld = struct(...
    'file','',...
    'id',nan,...
    'boundingBox',[],...
    'image',imcell,...
    'mask',cellmask,...
    'meanIntensity',nan,...
    'stdIntensity',nan,...
    'modeBackgroundIntensity',nan,...
    'mesh',[],...
    'meshArea',nan,...
    'meshSolidity',nan,...
    'meshRibsLength',[],...
    'meshRibsLengthMax',nan,...
    'meshRibsLengthMean',nan,...
    'meshRibsLengthStd',nan,...
    'meshMedialAxis',[],...
    'meshLength',nan,...
    'meshPoleCoordinates',[],...
    'meshSinuosity',nan,...
    'meshCurvature',[],...
    'meshCurvatureMean',nan,...
    'meshCurvatureMax',nan,...
    'meshCurvatureInflection',nan,...
    'meshVolume',nan,...
    'flagellumId',[]...
    );

cellmask = bwmorph(cellmask,'thicken',2);
imc = contourc(imgaussfilt(imcell.*cellmask,3),1);
imc = imc(:,2:end);

[xq, yq] = meshgrid(1:size(imcell,2),1:size(imcell,1));
in = inpolygon(xq-0.5,yq-0.5,imc(1,:),imc(2,:));

celld.meanIntensity = mean(imcell(in));
celld.stdIntensity = std(imcell(in));
celld.modeBackgroundIntensity = mode(round(imcell(~in)));

celld.meshArea = polyarea(imc(1,:),imc(2,:));
ch = convhull(imc');
celld.meshSolidity = celld.meshArea/polyarea(imc(1,ch),imc(2,ch));

[skel, mesh] = model2mesh(imc',0.5,0.001,100);
celld.mesh = mesh;

if mesh
    celld.meshRibsLength = sqrt((mesh(:,3)-mesh(:,1)).^2 + (mesh(:,4)-mesh(:,2)).^2);
    celld.meshRibsLengthMax = max(celld.meshRibsLength);
    celld.meshRibsLengthMean = nanmean(celld.meshRibsLength);
    celld.meshRibsLengthStd = nanstd(celld.meshRibsLength);
    
    [xpol, ypol, ipol] = polyxpoly(skel(:,1),skel(:,2),imc(1,:),imc(2,:));
    
    if numel(xpol)==2
        
        [ipol, iipol] = sort(ipol(:,1));
        skel = [xpol(iipol(1)), ypol(iipol(1)); skel((ipol(1)+1):ipol(2),:); xpol(iipol(2)), ypol(iipol(2))];
        
        celld.meshMedialAxis = skel;
        celld.meshPoleCoordinates = [ypol, xpol];
        celld.meshLength = sum(sqrt(sum(diff(skel).^2,2)));
        
        celld.meshSinuosity = celld.meshLength / sqrt(sum(diff([xpol, ypol]).^2,2));
        
        va = diff(skel);
        vb = va(2:end,:);
        va = va(1:end-1,:);
        ang = sign(atan2d(va(:,1).*vb(:,2)-va(:,2).*vb(:,1),va(:,1).*vb(:,1)+va(:,2).*vb(:,2)));
        
        skelRadius = zeros(size(skel,1)-2,1);
        
        for r = 1:numel(skelRadius)
            [~,~,skelRadius(r),~] = circfit(skel(r:r+2,1),skel(r:r+2,2));
        end
        
        celld.meshCurvature = [NaN; ang .* 1./skelRadius; NaN];
        celld.meshCurvatureMax = max(abs(celld.meshCurvature));
        celld.meshCurvatureMean = nanmean(abs(celld.meshCurvature));
        celld.meshCurvatureInflection = sum(ang(2:end) ~= ang(1:end-1));
        
        celld.meshVolume = celld.meshLength * mean(celld.meshRibsLength/2)^2 * pi;
        
    end
    
%     imshow(labeloverlay(uint8(imcell),cellmask),[]);
%     hold on;
%     plot(imc(1,:),imc(2,:),'g-');
%     plot(skel(:,1),skel(:,2),'g-');
%     plot(xpol,ypol,'ro');
%     hold off;
end
end