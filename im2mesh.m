imc = contourc(impg,[1 1] * max([-100; impg(impg<0)]));
imc = imc(:,2:end);

celld.cellMeshArea = polyarea(imc(1,:),imc(2,:));
ch = convhull(imc');
celld.cellMeshSolidity = celld.cellMeshArea/polyarea(imc(1,ch),imc(2,ch));

[skel, cellmesh] = model2mesh(imc',0.5,10^-(CONST.superSeggerOpti.MAGIC_RADIUS-1),celld.cellLength(2));
[xpol, ypol, ipol] = polyxpoly(skel(:,1),skel(:,2),imc(1,:),imc(2,:));

[ipol, iipol] = sort(ipol(:,1));

skel = [xpol(iipol(1)), ypol(iipol(1)); skel((ipol(1)+1):ipol(2),:); xpol(iipol(2)), ypol(iipol(2))];
celld.cellMeshMedialAxis = skel;
celld.cellMesh = cellmesh;
celld.cellMeshLength = sum(sqrt(sum(diff(skel).^2,2)));
celld.cellMeshRibsLength = sqrt((cellmesh(:,3)-cellmesh(:,1)).^2 + (cellmesh(:,4)-cellmesh(:,2)).^2);
celld.cellMeshRibsLengthMax = max(celld.cellMeshRibsLength);
celld.cellMeshRibsLengthMean = nanmean(celld.cellMeshRibsLength);
celld.cellMeshRibsLengthStd = nanstd(celld.cellMeshRibsLength);
celld.cellMeshSinuosity = celld.cellMeshLength / sqrt(sum(diff([xpol, ypol]).^2,2));

va = diff(skel);
vb = va(2:end,:);
va = va(1:end-1,:);
ang = sign(atan2d(va(:,1).*vb(:,2)-va(:,2).*vb(:,1),va(:,1).*vb(:,1)+va(:,2).*vb(:,2)));

skelRadius = zeros(size(skel,1)-2,1);
for k = 1:numel(skelRadius)
    [~,~,skelRadius(k),~] = circfit(skel(k:k+2,1),skel(k:k+2,2));
end
celld.cellMeshCurvature = [NaN; ang .* 1./skelRadius; NaN];
celld.cellMeshCurvatureMax = max(abs(celld.cellMeshCurvature));
celld.cellMeshCurvatureMean = nanmean(abs(celld.cellMeshCurvature));
celld.cellMeshCurvatureInflection = sum(ang(2:end) ~= ang(1:end-1));

celld.cellMeshVolume = celld.cellMeshLength * mean(celld.cellMeshRibsLength/2)^2 * pi;

imshow(celld.phase,[]);
hold on;
plot(imc(1,:),imc(2,:),'w-');
plot(skel(:,1),skel(:,2),'w-');
plot(xpol,ypol,'wo');
for k = 1:size(cellmesh,1)
    plot(cellmesh(k,[1 3]),cellmesh(k,[2 4]),'w:');
end
c = [NaN; 1./skelRadius / (max(1./skelRadius)); NaN];
scatter(skel(ang==1,1),skel(ang==1,2),50*c(ang==1)+1,'g');
scatter(skel(ang==-1,1),skel(ang==-1,2),50*c(ang==-1)+1,'r');
axis image;
hold off;
drawnow;
pause