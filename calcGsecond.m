function [G,Gj0,flag] = calcGsecond(j,k,w,epsimp,lambda)
% Constants
Omega = 1.0;
t = 1.0;
g = sqrt(t*Omega*lambda/2); %positive only! (for -imag checking)

G = zeros(length(lambda),length(w));

disc = w(2)-w(1);
wlength = length(w);
wNlength = floor(Omega/disc);

N = 50;

wMin = w(1)-N*Omega;
wAll = wMin:disc:w(wlength);

G0All = calcG0(wAll);
G1All = calcGt(1,0,wAll);
G2All = calcGt(2,0,wAll);

Gv = zeros(2,length(w));
Gj0 = zeros(2,length(w));
flag = 0;

alpha = eye(2);

for i = 1:length(w)
    for N = 50
        
        n = N;
        
        diffN = N-n;
        G0 = G0All(diffN*wNlength+i);
        G1 = G1All(diffN*wNlength+i);
        G2 = G2All(diffN*wNlength+i);
        
        beta = g*[2*G1, 2*G0; ...
            (G0+G2),2*G1];
        
        impM = epsimp*[G0,0; ...
            G1, 0];
        
        An = (alpha-impM)\(n*beta);
        
        for n = n-1:-1:1
            
            diffN = N-n;
            G0 = G0All(diffN*wNlength+i);
            G1 = G1All(diffN*wNlength+i);
            G2 = G2All(diffN*wNlength+i);
            
            impM = epsimp*[G0,0; ...
                G1, 0];
            
            beta = g*[2*G1, 2*G0; ...
                (G0+G2),2*G1];
            
            An = (alpha-impM-beta*An)\(n*beta);
            
        end
        
    end
    n = 0;
    wn = w(i);
    diffN = N-n;
    G0 = G0All(diffN*wNlength+i);
    G1 = G1All(diffN*wNlength+i);
    G2 = G2All(diffN*wNlength+i);
        
    impM = epsimp*[G0,0; ...
        G1, 0];
    
    beta = g*[2*G1, 2*G0; ...
        (G0+G2),2*G1];
    
    gamma = [calcGt(j,0,wn);calcGt(j,1,wn)];
    
    Gj0(:,i) = (alpha-impM-beta*An)\gamma;
    
    if(k ==0)
        Gv(:,i) = Gj0(:,i);
    else
        G0k = calcGt(0,k,wn);
        G1k = calcGt(1,k,wn);
        Gm1k = calcGt(-1,k,wn);
        
        alpha(1,1) = epsimp*G0;
        alpha(2,2) = epsimp*G1;
        gamma = [G0k;G1k+Gm1k];
        
        Gv(:,i) = gamma+(alpha-beta*An)*Gj0;
    end
    
    G(i) = Gv(1,i)/pi;
end
    
end