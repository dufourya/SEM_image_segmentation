function fI = filterFlagella(I, sigma)

F = imgaussfilt(I,sigma);

[Dx, Dy] = gradient(F);
[Hxx, ~] = gradient(Dx);
[Hxy, Hyy] = gradient(Dy);

tmp = sqrt((Hxx - Hyy).^2 + 4*Hxy.^2);

% Compute the eigenvalues
% mu1 = 0.5*(Hxx + Hyy + tmp);
mu2 = 0.5*(Hxx + Hyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
% check=abs(mu1)>abs(mu2);

% Lambda2=mu2; Lambda2(check)=mu1(check);
% fI = -Lambda2;
% fI(fI<0)=0;

fI=-mu2;