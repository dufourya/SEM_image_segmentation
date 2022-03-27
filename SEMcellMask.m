function bw = SEMcellMask(im,cellradius,FDR,cellarea)

imk = imgaussfilt(im,1);
% imk2 = imgaussfilt(im,3);

% imk = imgradient(imk);
% imk = imk-min(imk(:));
% imk = imk/max(imk(:));
% bw = imbinarize(imk,'adaptive','Sensitivity',0);

ev = fitgmdist(imk(1:100:end)',2);
[cp,idx]=sort(ev.ComponentProportion);
f = @(x) FDR-1./(1+(cp(1)*(1-cdf('Normal',x,ev.mu(idx(1)),sqrt(ev.Sigma(:,:,idx(1))))))./(cp(2)*(1-cdf('Normal',x,ev.mu(idx(2)),sqrt(ev.Sigma(:,:,idx(2)))))));
t = fzero(f,0);
%     histogram(img(1:10:end),'Normalization','pdf'); hold on;
%     plot(-20:0.1:10,ev.pdf([-20:0.1:10]'));
%     figure, imshow(img>=t,[]);

bw = imk>t;

bw = activecontour(imk,bw);
bw = bwmorph(bw,'thin',1);

bw = imopen(bw,strel('disk',cellradius,8));

bw = bwmorph(bw,'hbreak');
bw = bwmorph(bw,'spur');
% bw = ~imclose(~bw,strel('disk',10,8));

% bw = ~bwmorph(~bw,'diag');
%
bw = ~bwmorph(~bw,'bridge');
bw = bwareaopen(bw,cellarea);

bw = bwmorph(bw,'thicken',2);
bw = bwlabel(bw);

% imshow(labeloverlay(imk,bwl),[])

end