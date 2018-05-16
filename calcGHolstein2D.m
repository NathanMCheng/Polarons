function [G00] = calcGHolstein2D(w,epsimpIn,lambdaIn)
% Constants
Omega = 1.0;
t = 1.0;

G00 = zeros(length(w),length(epsimpIn),length(lambdaIn));

disc = w(2)-w(1);
wlength = length(w);
wNlength = floor(Omega/disc);

Nmax = 50;

wMin = w(1)-Nmax*Omega;
wAll = wMin:disc:w(wlength);

[G0Allinit,~,~] = calcGt2D(wAll);
epsInd = 0;
for epsimp = epsimpIn
    epsInd = epsInd+1;
    G0All = G0Allinit./(1-epsimp*G0Allinit);
    
    lambInd = 0;
    for lambda = lambdaIn
        lambInd = lambInd+1;
        g = sqrt(t*Omega*lambda/4); %positive only! (for -imag checking)
        
        for N = Nmax:-1:1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make the (n-1) map
            Nalpha = N-1;
            i = 1;
            %p = 0
            m = floor(Nalpha/2);
            while m>0
                n = Nalpha-m;
                alphaKeySet{i} = int2str([m n]);
                alphaValueSet{i} = i;
                i = i+1;
                
                m = m-1;
            end
            alphaKeySet{i} = int2str([m Nalpha]);
            alphaValueSet{i} = i;
            
            alphaMap = containers.Map(alphaKeySet,alphaValueSet);
            clear alphaKeySet alphaValueSet
            
            if (N==Nmax)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Make the (n+1) map
                Nbeta = N+1;
                i = 1;
                %p = 0
                m = floor(Nbeta/2);
                while m>0
                    n = Nbeta-m;
                    betaKeySet{i} = int2str([m n]);
                    betaValueSet{i} = i;
                    i = i+1;
                    
                    m = m-1;
                end
                betaKeySet{i} = int2str([m Nbeta]);
                betaValueSet{i} = i;
                
                betaMap = containers.Map(betaKeySet,betaValueSet);
                clear betaKeySet betaValueSet
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make the alpha and beta matrices
            
            i = 1;
            jBeta = 1;
            jAlpha = 1;
            %p = 0
            m = floor(N/2);
            while m>0
                n = N-m;
                
                % Beta Matrix indexes
                
                if n>m
                    betaIndex1(jBeta) = i;
                    betaIndex2(jBeta) = betaMap(int2str([m+1 n]));
                    betaIndex3(jBeta) = 1;
                    jBeta = jBeta+1;
                end
                
                betaIndex1(jBeta) = i;
                betaIndex2(jBeta) = betaMap(int2str([m n+1]));
                if (n==m)
                    betaIndex3(jBeta) = 2;
                else
                    betaIndex3(jBeta) = 1;
                end
                jBeta = jBeta+1;
                
                % Alpha Matrix
                if n>m
                    alphaIndex1(jAlpha) = i;
                    alphaIndex2(jAlpha) = alphaMap(int2str([m n-1]));
                    alphaIndex3(jAlpha) = 1;
                    jAlpha = jAlpha+1;
                end
                
                alphaIndex1(jAlpha) = i;
                alphaIndex2(jAlpha) = alphaMap(int2str([m-1 n]));
                if (m==n)
                    alphaIndex3(jAlpha) = 2;
                else
                    alphaIndex3(jAlpha) = 1;
                end
                jAlpha = jAlpha+1;
                
                
                KeySet{i} = int2str([m n]);
                ValueSet{i} = i;
                i = i+1;
                
                m = m-1;
            end
            
            
            % p = 0, m = 0;
            n = N;
            
            % Beta Matrix indexes
            
            betaIndex1(jBeta) = i;
            betaIndex2(jBeta) = betaMap(int2str([m+1 n]));
            betaIndex3(jBeta) = 2;
            jBeta = jBeta+1;
            
            betaIndex1(jBeta) = i;
            betaIndex2(jBeta) = betaMap(int2str([m n+1]));
            betaIndex3(jBeta) = 1;
            jBeta = jBeta+1;
            
            % Alpha Matrix
            alphaIndex1(jAlpha) = i;
            alphaIndex2(jAlpha) = alphaMap(int2str([m n-1]));
            alphaIndex3(jAlpha) = 1;
            jAlpha = jAlpha+1;
            
            KeySet{i} = int2str([m n]);
            ValueSet{i} = i;
            
            Map = containers.Map(KeySet,ValueSet);
            clear KeySet ValueSet
            
            alpha = sparse(alphaIndex1,alphaIndex2,alphaIndex3);
            beta = sparse(betaIndex1, betaIndex2, betaIndex3);
            gamma = eye(i);
            
            clear alphaIndex1 alphaIndex2 alphaIndex3
            clear betaIndex1 betaIndex2 betaIndex3
            
            for j = 1:length(w)
                value = g*G0All((Nmax-N)*wNlength+1);
                
                if (N==Nmax)
                    temp = gamma\(value*alpha);
                    [sm,sn] = size(temp);
                    An(1:sm,1:sn,j) = temp;
                    %             [si,sj,sk] = find(temp);
                    %             l(j) = length(si);
                    %             An(1:l(j),1,j) = si;
                    %             An(1:l(j),2,j) = sj;
                    %             An(1:l(j),3,j) = sk;
                else
                    %             A = sparse(An(1:l(j),1,j),An(1:l(j),2,j),An(1:l(j),3,j),sm,sn);
                    temp = (gamma-value*beta*An(:,:,j))\(value*alpha);
                    [sm,sn] = size(temp);
                    An(1:sm,1:sn,j) = temp;
                    %             [si,sj,sk] = find(temp);
                    %             l(j) = length(si);
                    %             An(1:l(j),1,j) = si;
                    %             An(1:l(j),2,j) = sj;
                    %             An(1:l(j),3,j) = sk;
                end
                
            end
            [l,w] = size(temp);
            An = squeeze(An(1:l,1:w,:));
            
            betaMap = Map;
            Map = alphaMap;
            clear alphaMap
            
        end
        G00(:,epsInd,lambInd) = G0All((Nmax*wNlength+1):end)./(1-g*G0All((Nmax*wNlength+1):end).*An);
    end
end
end
