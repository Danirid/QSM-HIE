function x = average(chi4D, mag, TEs)
R = 1/(TEs(1)-TEs(2))*log(mag(:,:,:,2)./mag(:,:,:,1));
R(isinf(R)|isnan(R)) = 0;
siz = size(chi4D);
n = zeros(siz(1:3)); % numerateur
d = zeros(siz(1:3)); %denominateur
for i = 1:length(TEs)
    w_i = TEs(i)*exp(-TEs(i)*R);
    n = n + w_i.^2.*chi4D(:,:,:,i);
    d = d + w_i.^2;
end
x = n./d;
end