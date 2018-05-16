function [G00,A1,A2] = calcGt3D(z)
%calcGt calculates hopping Gnm(w)
%   given w, calculates hoppin Gnm(w)

t = 1.0;
eta = 0.1; %width of delta peaks

z = z+1i*eta;

A1 = zeros(1,length(z));
A2 = zeros(2,1,length(z));
Nmax = 50;

for N = Nmax:-1:1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make the (n-1) map
    Nalpha = N-1;
    i = 1;
    p = floor(Nalpha/3);
    while p>0
        m = floor((Nalpha-p)/2);
        n = Nalpha-m-p;
        alphaKeySet{i} = int2str([p m n]);
        alphaValueSet{i} = i;
        i = i+1;
        while m>p
            m = m-1;
            n = n+1;
            alphaKeySet{i} = int2str([p m n]);
            alphaValueSet{i} = i;
            i = i+1;
        end
        p = p-1;
    end
    %p = 0
    m = floor(Nalpha/2);
    while m>0
        n = Nalpha-m;
        alphaKeySet{i} = int2str([p m n]);
        alphaValueSet{i} = i;
        i = i+1;
        
        m = m-1;
    end
    alphaKeySet{i} = int2str([p m Nalpha]);
    alphaValueSet{i} = i;
    
    alphaMap = containers.Map(alphaKeySet,alphaValueSet);
    clear alphaKeySet alphaValueSet
    
    if (N==Nmax)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the (n+1) map
        Nbeta = N+1;
        i = 1;
        p = floor(Nbeta/3);
        while p>0
            m = floor((Nbeta-p)/2);
            n = Nbeta-m-p;
            betaKeySet{i} = int2str([p m n]);
            betaValueSet{i} = i;
            i = i+1;
            while m>p
                m = m-1;
                n = n+1;
                betaKeySet{i} = int2str([p m n]);
                betaValueSet{i} = i;
                i = i+1;
            end
            p = p-1;
        end
        %p = 0
        m = floor(Nbeta/2);
        while m>0
            n = Nbeta-m;
            betaKeySet{i} = int2str([p m n]);
            betaValueSet{i} = i;
            i = i+1;
            
            m = m-1;
        end
        betaKeySet{i} = int2str([p m Nbeta]);
        betaValueSet{i} = i;
        
        betaMap = containers.Map(betaKeySet,betaValueSet);
        clear betaKeySet betaValueSet
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make the alpha and beta matrices
    
    i = 1;
    jBeta = 1;
    jAlpha = 1;
    p = floor(N/3);
    while p>0
        m = floor((N-p)/2);
        n = N-m-p;
        
        % Beta Matrix indexes
        if m>p
            betaIndex1(jBeta) = i;
            betaIndex2(jBeta) = betaMap(int2str([p+1 m n]));
            betaIndex3(jBeta) = 1;
            jBeta = jBeta+1;
        end
        if n>m
            betaIndex1(jBeta) = i;
            betaIndex2(jBeta) = betaMap(int2str([p m+1 n]));
            if(m==p)
                betaIndex3(jBeta) = 2;
            else
                betaIndex3(jBeta) = 1;
            end
            jBeta = jBeta+1;
        end
        
        betaIndex1(jBeta) = i;
        betaIndex2(jBeta) = betaMap(int2str([p m n+1]));
        if n==m
            if m == p
                betaIndex3(jBeta) = 3;
            else
                betaIndex3(jBeta) = 2;
            end
        else
            betaIndex3(jBeta) = 1;
        end
        jBeta = jBeta+1;
        
        % Alpha Matrix
        if n>m
            alphaIndex1(jAlpha) = i;
            alphaIndex2(jAlpha) = alphaMap(int2str([p m n-1]));
            alphaIndex3(jAlpha) = 1;
            jAlpha = jAlpha+1;
        end
        
        if m>p
            alphaIndex1(jAlpha) = i;
            alphaIndex2(jAlpha) = alphaMap(int2str([p m-1 n]));
            if (m==n)
                alphaIndex3(jAlpha) = 2;
            else
                alphaIndex3(jAlpha) = 1;
            end
            jAlpha = jAlpha+1;
        end
        
        alphaIndex1(jAlpha) = i;
        alphaIndex2(jAlpha) = alphaMap(int2str([p-1 m n]));
        if p == m
            if m == n
                alphaIndex3(jAlpha) = 3;
            else
                alphaIndex3(jAlpha) = 2;
            end
        else
            alphaIndex3(jAlpha) = 1;
        end
        jAlpha = jAlpha+1;
        
        KeySet{i} = int2str([p m n]);
        ValueSet{i} = i;
        i = i+1;
        while m>p
            m = m-1;
            n = n+1;
            
            % Beta Matrix indexes
            if m>p
                betaIndex1(jBeta) = i;
                betaIndex2(jBeta) = betaMap(int2str([p+1 m n]));
                betaIndex3(jBeta) = 1;
                jBeta = jBeta+1;
            end
            if n>m
                betaIndex1(jBeta) = i;
                betaIndex2(jBeta) = betaMap(int2str([p m+1 n]));
                if (m==p)
                    betaIndex3(jBeta) = 2;
                else
                    betaIndex3(jBeta) = 1;
                end
                jBeta = jBeta+1;
            end
            
            betaIndex1(jBeta) = i;
            betaIndex2(jBeta) = betaMap(int2str([p m n+1]));
            if n==m
                if m == p
                    betaIndex3(jBeta) = 3;
                else
                    betaIndex3(jBeta) = 2;
                end
            else
                betaIndex3(jBeta) = 1;
            end
            jBeta = jBeta+1;
            
            % Alpha Matrix
            if n>m
                alphaIndex1(jAlpha) = i;
                alphaIndex2(jAlpha) = alphaMap(int2str([p m n-1]));
                alphaIndex3(jAlpha) = 1;
                jAlpha = jAlpha+1;
            end
            if m>p
                alphaIndex1(jAlpha) = i;
                alphaIndex2(jAlpha) = alphaMap(int2str([p m-1 n]));
                if (n==m)
                    alphaIndex3(jAlpha) = 2;
                else
                    alphaIndex3(jAlpha) = 1;
                end
                jAlpha = jAlpha+1;
            end
            
            alphaIndex1(jAlpha) = i;
            alphaIndex2(jAlpha) = alphaMap(int2str([p-1 m n]));
            if p == m
                if m == n
                    alphaIndex3(jAlpha) = 3;
                else
                    alphaIndex3(jAlpha) = 2;
                end
            else
                alphaIndex3(jAlpha) = 1;
            end
            jAlpha = jAlpha+1;
            
            
            KeySet{i} = int2str([p m n]);
            ValueSet{i} = i;
            i = i+1;
        end
        p = p-1;
    end
    %p = 0
    m = floor(N/2);
    while m>0
        n = N-m;
        
        % Beta Matrix indexes
        
        betaIndex1(jBeta) = i;
        betaIndex2(jBeta) = betaMap(int2str([p+1 m n]));
        betaIndex3(jBeta) = 2;
        jBeta = jBeta+1;
        
        if n>m
            betaIndex1(jBeta) = i;
            betaIndex2(jBeta) = betaMap(int2str([p m+1 n]));
            betaIndex3(jBeta) = 1;
            jBeta = jBeta+1;
        end
        
        betaIndex1(jBeta) = i;
        betaIndex2(jBeta) = betaMap(int2str([p m n+1]));
        if (n==m)
            betaIndex3(jBeta) = 2;
        else
            betaIndex3(jBeta) = 1;
        end
        jBeta = jBeta+1;
        
        % Alpha Matrix
        if n>m
            alphaIndex1(jAlpha) = i;
            alphaIndex2(jAlpha) = alphaMap(int2str([p m n-1]));
            alphaIndex3(jAlpha) = 1;
            jAlpha = jAlpha+1;
        end
        
        alphaIndex1(jAlpha) = i;
        alphaIndex2(jAlpha) = alphaMap(int2str([p m-1 n]));
        if (m==n)
            alphaIndex3(jAlpha) = 2;
        else
            alphaIndex3(jAlpha) = 1;
        end
        jAlpha = jAlpha+1;
        
        
        KeySet{i} = int2str([p m n]);
        ValueSet{i} = i;
        i = i+1;
        
        m = m-1;
    end
    
    
    % p = 0, m = 0;
    n = N;
    
    % Beta Matrix indexes
    
    betaIndex1(jBeta) = i;
    betaIndex2(jBeta) = betaMap(int2str([p m+1 n]));
    betaIndex3(jBeta) = 4;
    jBeta = jBeta+1;
    
    betaIndex1(jBeta) = i;
    betaIndex2(jBeta) = betaMap(int2str([p m n+1]));
    betaIndex3(jBeta) = 1;
    jBeta = jBeta+1;
    
    % Alpha Matrix
    alphaIndex1(jAlpha) = i;
    alphaIndex2(jAlpha) = alphaMap(int2str([p m n-1]));
    alphaIndex3(jAlpha) = 1;
    jAlpha = jAlpha+1;
    
    KeySet{i} = int2str([p m n]);
    ValueSet{i} = i;
    
    Map = containers.Map(KeySet,ValueSet);
    clear KeySet ValueSet
    
    alpha = sparse(alphaIndex1,alphaIndex2,alphaIndex3);
    beta = sparse(betaIndex1, betaIndex2, betaIndex3);
    gamma = speye(i);
        
    clear alphaIndex1 alphaIndex2 alphaIndex3
    clear betaIndex1 betaIndex2 betaIndex3
    
    for j = 1:length(z)
        value = -t/z(j);
        
        if (N==Nmax)
            temp = gamma\(value*alpha);
            [sm,sn] = size(temp);
            [si,sj,sk] = find(temp);
            l(j) = length(si);
            An(1:l(j),1,j) = si;
            An(1:l(j),2,j) = sj;
            An(1:l(j),3,j) = sk;
        else
            A = sparse(An(1:l(j),1,j),An(1:l(j),2,j),An(1:l(j),3,j),sm,sn);
            temp = (gamma-value*beta*A)\(value*alpha);
            [sm,sn] = size(temp);
            [si,sj,sk] = find(temp);
            l(j) = length(si);
            An(1:l(j),1,j) = si;
            An(1:l(j),2,j) = sj;
            An(1:l(j),3,j) = sk;
        end

    end
%     [l,w,d] = size(temp);
%     An = An(1:l,1:w,1:d,:);
    if(N==1)
        A1 = squeeze(An(:,:,:));
    elseif(N==2)
        A2 = squeeze(An(:,:,:));
    end
    betaMap = Map;
    Map = alphaMap;
    clear alphaMap

end
G00 = 1./(z+6*t*A1);