function bw = SEMflagellaMask(im,cellmask,flagellawidth,FDR, minsize)


% flagellawidth = 0.75;
% minsize = 25;

F = imgaussfilt(im,flagellawidth);

[Dx, Dy] = gradient(F);
[Hxx, ~] = gradient(Dx);
[Hxy, Hyy] = gradient(Dy);

img = -0.5*(Hxx + Hyy - sqrt((Hxx - Hyy).^2 + 4*Hxy.^2));

%%
ev = fitgmdist(img(1:100:end)',2);
[cp,idx]=sort(ev.ComponentProportion);
f = @(x) FDR-1./(1+(cp(1)*(1-cdf('Normal',x,ev.mu(idx(1)),sqrt(ev.Sigma(:,:,idx(1))))))./(cp(2)*(1-cdf('Normal',x,ev.mu(idx(2)),sqrt(ev.Sigma(:,:,idx(2)))))));
t = fzero(f,0);
%     histogram(img(1:10:end),'Normalization','pdf'); hold on;
%     plot(-20:0.1:10,ev.pdf([-20:0.1:10]'));
%     figure, imshow(img>=t,[]);

%%
% close all
bw = img>t;
bw = bwmorph(bw,'clean');
bw = bwmorph(bw,'bridge');
bw = bwmorph(bw,'fill',Inf);
bw = bwmorph(bw,'spur');
bw = bwmorph(bw,'hbreak');
bw = imdilate(bw,[0 1 0;1 1 1; 0 1 0]);
% bw = imclose(bw,[1 1 1;1 1 1; 1 1 1]);
bw = bwmorph(bw,'open');
% bw = bwmorph(bw,'thin',1);
bw = bwareaopen(bw,minsize);
bw(cellmask>0) = 0;
bw = bwareaopen(bw,round(minsize/10));

% imshow(labeloverlay(uint8(im),bw+2*bw1),[])
%%

bw = bwlabel(bw);