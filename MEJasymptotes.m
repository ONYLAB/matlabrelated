function m = MEJasymptotes(t,D,k,s)

T = D.*t./s./s;
g = -psi(1);

dn = (log(4*T)-2*g+4*pi*D/k);
dn2 = dn.*dn;
dn3 = dn2.*dn;
dn4 = dn2.*dn2;

m = 4*pi*D.*(1./dn-g./dn2-1.311./dn3+0.25./dn4);

ind1 = intersect(find(t>0.001*s^2/D),find(t<10*s^2/D));
m(ind1) = nan;

ind2 = find(t<=0.001*s^2/D);

ts = t(ind2);

T = D.*ts./s./s;

if isinf(k)
    m(ind2) = 2*pi*D*(1./(sqrt(pi*T)) + 1/2 - 1/4*sqrt(T/pi) + 1/8*T - 100/384*T.*sqrt(T/pi));
else
    m(ind2) = 2*pi*D*(1/2 - k/(2*pi*D)*sqrt(T/pi) - k/(4*pi*D)*T.*sqrt(T/pi) + k/(2*pi*D)*erfcx(k/(2*pi*D)*sqrt(T)).*(1 - pi*D/k + k/(2*pi*D)*T + 1/2*(k/(2*pi*D)*T).^2)  );
end
%     m(ind2) = 2*pi*D*(1/2 - k/(2*pi*D)*sqrt(T/pi) - k/(4*pi*D)*T.*sqrt(T/pi) + k/(2*pi*D)*exp(T*(2*pi*D)^2).*erfc(k/(2*pi*D)*sqrt(T)).*(1 - pi*D/k + k/(2*pi*D)*T + 1/2*(k/(2*pi*D)*T).^2)  );
