function [G] = calcGSSH2D(w,epsimpIn,lambdaIn)
% Constants
Omega = 1.0;
t = 1.0;

G = zeros(length(w),length(epsimpIn),length(lambdaIn));

disc = w(2)-w(1);
wlength = length(w);
wNlength = floor(Omega/disc);

N = 100;

wMin = w(1)-N*Omega;
wAll = wMin:disc:w(wlength);

[G0All,A1,A2] = calcGt2D(wAll);
G1 = A1.*G0All;
G2All(:,1) = transpose(A2(2,:)).*G1;

epsInd = 0;
for epsimp = epsimpIn
    epsInd = epsInd+1;
    
    lambInd = 0;
    for lambda = lambdaIn
        lambInd = lambInd+1;
        g = sqrt(t*Omega*lambda/4);
        
        n = N;
        
        diffN = N-n;
        diffNm1 = N-n+1;
        
        G0 = G0All((diffN*wNlength+1):(diffN*wNlength+wlength));
        
        an = g*G0;
        
        G0m1 = G0All((diffNm1*wNlength+1):(diffNm1*wNlength+wlength));
        
        G2m1 = G2All((diffNm1*wNlength+1):(diffNm1*wNlength+wlength));
        
        bm1 = 2*g*(G0m1-G2m1);
        
        % bp1 = g*(-1./(1-calcGt2D(0,0,wp1,0)+(wp1+1i*eta)./t^2))./(1-epsimp*calcGt2D(0,0,wp1,0));
        % bm1 = g*(-1./(1-calcGt2D(0,0,wm1,0)+(wm1+1i*eta)./t^2))./(1-epsimp*calcGt2D(0,0,wm1,0));
        
        An = n*(n-1)*an.*bm1 ...
            ./(1-G0.*(g*(bm1*n)+epsimp));
        
        for n = n-2:-2:2
            diffN = N-n;
            diffNp1 = N-n-1;
            diffNm1 = N-n+1;
            
            G0 = G0All((diffN*wNlength+1):(diffN*wNlength+wlength));
            
            an = g*G0;
            
            G0p1 = G0All((diffNp1*wNlength+1):(diffNp1*wNlength+wlength));
            G0m1 = G0All((diffNm1*wNlength+1):(diffNm1*wNlength+wlength));
            
            G2p1 = G2All((diffNp1*wNlength+1):(diffNp1*wNlength+wlength));
            G2m1 = G2All((diffNm1*wNlength+1):(diffNm1*wNlength+wlength));
            
            bp1 = 2*g*(G0p1-G2p1);
            bm1 = 2*g*(G0m1-G2m1);
            
            
            An = n*(n-1)*an.*bm1 ...
                ./(1-G0.*(g*((An+n+1).*bp1+n*bm1)+epsimp));
        end
        
        %         A(i,N) = An;
        
        
        n = 0;
        diffN = N-n;
        diffNp1 = N-n-1;
        G0 = G0All((diffN*wNlength+1):(diffN*wNlength+wlength));
        
        G0p1 = G0All((diffNp1*wNlength+1):(diffNp1*wNlength+wlength));
        
        G2p1 = G2All((diffNp1*wNlength+1):(diffNp1*wNlength+wlength));
        
        bp1 = 2*g*(G0p1-G2p1);
        
        G(:,epsInd,lambInd) = G0./(1-G0.*(2*g*(bp1.*(1+An))+epsimp));
        
        
    end
    
end
end

