function [G,Gj0] = calcGSSH(j,k,w,epsimp,lambda)
% Constants
Omega = 1.0;
t = 1.0;
g = sqrt(t*Omega*lambda/2); %positive only! (for -imag checking)



for N = 100
    
    n = N;
    
    wn = w-n*Omega;
    G0 = (calcGt(0,0,wn));
    
    an = g*G0;
    
    wp1 = w-(n+1)*Omega;
    wm1 = w-(n-1)*Omega;
    
    G0p1 = calcGt(0,0,wp1);
    G0m1 = calcGt(0,0,wm1);
    
    G2p1 = calcGt(0,2,wp1);
    G2m1 = calcGt(0,2,wm1);
    
    bp1 = 2*g*(G0p1-G2p1);
    bm1 = 2*g*(G0m1-G2m1);

% bp1 = g*(-1./(1-calcGt(0,0,wp1,0)+(wp1+1i*eta)./t^2))./(1-epsimp*calcGt(0,0,wp1,0));
% bm1 = g*(-1./(1-calcGt(0,0,wm1,0)+(wm1+1i*eta)./t^2))./(1-epsimp*calcGt(0,0,wm1,0));

    An = n*(n-1)*an.*bm1 ...
        ./(1-G0.*(g*(bm1*n)+epsimp));
    
    for n = n-2:-2:2
        wn = w-n*Omega;
        G0 = (calcGt(0,0,wn));
        
        an = g*G0;
        
        wp1 = w-(n+1)*Omega;
        wm1 = w-(n-1)*Omega;
        
        G0p1 = calcGt(0,0,wp1);
        G0m1 = calcGt(0,0,wm1);
        
        G2p1 = calcGt(0,2,wp1);
        G2m1 = calcGt(0,2,wm1);
        
        bp1 = 2*g*(G0p1-G2p1);
        bm1 = 2*g*(G0m1-G2m1);
%         bp1 = g*(-1./(1-calcGt(0,0,wp1,0)+(wp1+1i*eta)./t^2))./(1-epsimp*calcGt(0,0,wp1,0));
%         bm1 = g*(-1./(1-calcGt(0,0,wm1,0)+(wm1+1i*eta)./t^2))./(1-epsimp*calcGt(0,0,wm1,0));


        An = n*(n-1)*an.*bm1 ...
            ./(1-G0.*(g*((An+n+1).*bp1+n*bm1)+epsimp));
    end
    
    %         A(i,N) = An;
    
end
n = 0;
wn = w-n*Omega;
G0 = (calcGt(0,0,wn));

wp1 = w-(n+1)*Omega;

G0p1 = calcGt(0,0,wp1);

G2p1 = calcGt(0,2,wp1);

bp1 = 2*g*(G0p1-G2p1);
% bp1 = g*(-1./(1-calcGt(0,0,wp1,0)+(wp1+1i*eta)./t^2))./(1-epsimp*calcGt(0,0,wp1,0));
        
Gj0 = calcGt(j,0,wn)./(1-G0.*(g*(bp1.*(1+An))+epsimp));

if(k ==0)
    G = Gj0;
else
    G0k = calcGt(0,k,wn);
    G = calcGt(j,k,wn)+G0k.*Gj0.*(g*(bp1.*(1+An))+epsimp);
end

G = G/pi;

end