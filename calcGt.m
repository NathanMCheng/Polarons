function [x] = calcGt(m,n,w)
%calcGt calculates hopping Gnm(w)
%   given w, calculates hoppin Gnm(w)

t = 1.0;
eta = 1e-3; %width of delta peaks


% A = 0.5*(-w/t+(sqrt(w/t-2).*sqrt(w/t+2)));

G0 = calcG0(w);

w = w+1i*eta;
G1 = (1-w.*G0)./(2*t);
if (n==m)
    Gmn = G0;
elseif abs(n-m)==1
    Gmn = G1;
else
    Gmn = G1;
    Gm2 = G0;
    Gm1 = G1;
    for i = 2:abs(n-m)
        Gt = Gmn;
        Gmn = -w.*Gm1./t-Gm2;
        Gm1 = Gmn;
        Gm2 = Gt;
    end
end

x = Gmn;

end

