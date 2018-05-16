function [G,Gj0] = calcGHolstein(j,k,w,epsimp,lambda)
% Constants
Omega = 1.0;
t = 1.0;
g = sqrt(t*Omega*lambda/2); %positive only! (for -imag checking)



for N = 100
    
    n = N;
    
    wn = w-n*Omega;
    G0 = (calcGt(0,0,wn));
    
    an = g*G0;
    
    An = n*an ...
        ./(1-G0.*(epsimp));
    
    for n = n-1:-1:1
        wn = w-n*Omega;
        G0 = (calcGt(0,0,wn));
        
        an = g*G0;
        
        An = n*an ...
            ./(1-G0.*(g*An+epsimp));
    end
    
    %         A(i,N) = An;
    
end
n = 0;
wn = w-n*Omega;
G0 = (calcGt(0,0,wn));
      
Gj0 = calcGt(j,0,wn)./(1-G0.*(g*An+epsimp));

if(j ==0)
    G = Gj0;
else
    G0k = calcGt(0,k,wn);
    G = calcGt(j,k,wn)+G0k.*Gj0.*(g*An+epsimp);
end

G = G/pi;

end