function [G,A] = calcGsecond2D(w,epsimpIn,lambdaIn)

Omega = 1.0;
t = 1.0;

G = zeros(length(w),length(epsimpIn),length(lambdaIn));

Gall = zeros(2,length(w));

Nmax = 10;

disc = w(2)-w(1);
wlength = length(w);
wNlength = floor(Omega/disc);

wMin = w(1)-Nmax*Omega;
wAll = wMin:disc:w(wlength);

[G0All,A1,A2] = calcGt2D(wAll);
G1All = G0All.*A1;
for i = 1:2
    G2All(:,i) = transpose(A2(i,:)).*G1All;
end
G2lateralAll = G2All(:,2);
G2diagonalAll = G2All(:,1);

epsInd = 0;
for epsimp = epsimpIn;
    epsInd = epsInd+1;
    
    lambdaInd = 0;
    for lambda = lambdaIn
        lambdaInd = lambdaInd+1;
        g = sqrt(t*Omega*lambda/4);
        
        for z = 1:length(w)
            
            % Set A_N = 0;
            An = 0;
            for N = Nmax:-1:1;
                Ndiff = Nmax-N;
                % Calculate Matrix Values
                G0 = G0All(Ndiff*wNlength+z);
                G1 = G1All(Ndiff*wNlength+z);
                G2lat = G2lateralAll(Ndiff*wNlength+z);
                G2diag = G2diagonalAll(Ndiff*wNlength+z);
                
                G0000 = G0;
                G1010 = G0;
                
                G0010 = G1;
                G1000 = G1;
                
                Gm1010 = G2lat;
                G0110 = G2diag;
                
                impValue = epsimp*[G0; G1];
                
                nValue = g*[2*G1000, 2*G0000; ...
                    (G1010+Gm1010), 2*G0010];
                
                mValue = 2*g*[G1000, G0000; ...
                    G0110, G0010];
                
                equalValue = g*[4*G1000, 4*G0000; ...
                    (G1010+Gm1010)+2*G0110, 4*G0010];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Make the (n-1) map
                Nalpha = N-1;
                i = 1;
                for site = 1:2
                    m = floor(Nalpha/2);
                    while m>0
                        n = Nalpha-m;
                        alphaKeySet{i} = int2str([m n site]);
                        alphaValueSet{i} = i;
                        i = i+1;
                        m = m-1;
                    end
                    %m = 0
                    n = Nalpha;
                    alphaKeySet{i} = int2str([m n site]);
                    alphaValueSet{i} = i;
                    i = i+1;
                end
                alphaMap = containers.Map(alphaKeySet,alphaValueSet);
                clear alphaKeySet alphaValueSet
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Make the (n+1) map
                if (N==Nmax)
                    Nbeta = N+1;
                    i = 1;
                    for site = 1:2
                        m = floor(Nbeta/2);
                        while m>0
                            n = Nbeta-m;
                            betaKeySet{i} = int2str([m n site]);
                            betaValueSet{i} = i;
                            i = i+1;
                            m = m-1;
                        end
                        %m = 0
                        n = Nbeta;
                        betaKeySet{i} = int2str([m n site]);
                        betaValueSet{i} = i;
                        i = i+1;
                    end
                    
                    betaMap = containers.Map(betaKeySet,betaValueSet);
                    clear betaKeySet betaValueSet
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Make the alpha and beta matrices
                
                i = 1;
                jBeta = 1;
                jAlpha = 1;
                
                for site = 1:2
                    m = floor(N/2);
                    while m>0
                        
                        n = N-m;
                        
                        % Beta Matrix indexes
                        if m<n
                            for site2 = 1:2
                                betaIndex1(jBeta) = i;
                                betaIndex2(jBeta) = betaMap(int2str([m+1 n site2]));
                                betaIndex3(jBeta) = mValue(site,site2);
                                jBeta = jBeta+1;
                            end
                        end
                        
                        for site2 = 1:2
                            betaIndex1(jBeta) = i;
                            betaIndex2(jBeta) = betaMap(int2str([m n+1 site2]));
                            if (m==n)
                                betaIndex3(jBeta) = equalValue(site,site2);
                            else
                                betaIndex3(jBeta) = nValue(site,site2);
                            end
                            jBeta = jBeta+1;
                        end
                        % Alpha Matrix indexes
                        if m<n
                            for site2 = 1:2
                                alphaIndex1(jAlpha) = i;
                                alphaIndex2(jAlpha) = alphaMap(int2str([m n-1 site2]));
                                alphaIndex3(jAlpha) = n*nValue(site,site2);
                                jAlpha = jAlpha+1;
                            end
                        end
                        
                        for site2 = 1:2
                            alphaIndex1(jAlpha) = i;
                            alphaIndex2(jAlpha) = alphaMap(int2str([m-1 n site2]));
                            if (m==n)
                                alphaIndex3(jAlpha) = m*equalValue(site,site2);
                            else
                                alphaIndex3(jAlpha) = m*mValue(site,site2);
                            end
                            jAlpha = jAlpha+1;
                        end
                        KeySet{i} = int2str([m n site]);
                        ValueSet{i} = i;
                        i = i+1;
                        
                        m = m-1;
                    end
                    
                    %%%%%%%%%%
                    
                    %m = 0;
                    n = N;
                    
                    % Beta Matrix indexes
                    if m<n
                        for site2 = 1:2
                            betaIndex1(jBeta) = i;
                            betaIndex2(jBeta) = betaMap(int2str([m+1 n site2]));
                            betaIndex3(jBeta) = mValue(site,site2);
                            jBeta = jBeta+1;
                        end
                    end
                    
                    for site2 = 1:2
                        betaIndex1(jBeta) = i;
                        betaIndex2(jBeta) = betaMap(int2str([m n+1 site2]));
                        if (m==n)
                            betaIndex3(jBeta) = equalValue(site,site2);
                        else
                            betaIndex3(jBeta) = nValue(site,site2);
                        end
                        jBeta = jBeta+1;
                    end
                    
                    % Alpha Matrix indexes
                    if m<n
                        for site2 = 1:2
                            alphaIndex1(jAlpha) = i;
                            alphaIndex2(jAlpha) = alphaMap(int2str([m n-1 site2]));
                            alphaIndex3(jAlpha) = n*nValue(site,site2);
                            jAlpha = jAlpha+1;
                        end
                    end
                    
                    
                    
                    KeySet{i} = int2str([m n site]);
                    ValueSet{i} = i;
                    i = i+1;
                    
                end
                i = i-1;
                
                
                gammaIndex1 = zeros(i,1);
                gammaIndex2 = zeros(i,1);
                gammaIndex3 = zeros(i,1);
                siteDivisor = i/2;
                for j = 1:siteDivisor
                    gammaIndex1(j) = j;
                    gammaIndex2(j) = j;
                    gammaIndex3(j) = impValue(1);
                    gammaIndex1(j+siteDivisor) = j+siteDivisor;
                    gammaIndex2(j+siteDivisor) = j;
                    gammaIndex3(j+siteDivisor) = impValue(2);
                end
                
                
                Map = containers.Map(KeySet,ValueSet);
                clear KeySet ValueSet
                
                alpha = sparse(alphaIndex1,alphaIndex2,alphaIndex3);
                beta = sparse(betaIndex1, betaIndex2, betaIndex3);
                gamma = speye(i)-sparse(gammaIndex1,gammaIndex2,gammaIndex3,i,i);
                
                clear alphaIndex1 alphaIndex2 alphaIndex3
                clear betaIndex1 betaIndex2 betaIndex3
                clear gammaIndex1 gammaIndex2 gammaIndex3
                
                betaMap = Map;
                Map = alphaMap;
                clear alphamap
                
                An = (gamma-beta*An)\alpha;
            end
            
            G0 = G0All(Nmax*wNlength+z);
            G1 = G1All(Nmax*wNlength+z);
            G2diag = G2diagonalAll(Nmax*wNlength+z);
            G2lat = G2lateralAll(Nmax*wNlength+z);
            
            alpha = [G0;G1];
            impM = epsimp*[G0,0; ...
                G1, 0];
            beta = g*[4*G1, 4*G0; ...
                (G0+G2lat+2*G2diag), 4*G1];
            
            Gall(:,z) = (eye(2)-impM-beta*An)\alpha;
            
        end
        G(:,epsInd,lambdaInd) = squeeze(Gall(1,:)./pi);
    end
end
end