function [G00,A1,A2] = calcGt2D(z)
%calcGt calculates hopping Gnm(w)
%   given w, calculates hoppin Gnm(w)

t = 1.0;
eta = 0.1; %width of delta peaks

%z = z+1i*eta;

A1 = zeros(1,length(z));
A2 = zeros(2,1,length(z));
Nmax = 100;

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
    
    for j = 1:length(z)
                
        if (abs(z(j)<4))
            NmaxZ = Nmax;
            value = -t/(z(j)+1i*eta);
        else
            NmaxZ = 30;
            value = -t/(z(j)+1i*eta*0.1);
        end
        
        if (N==NmaxZ)
            temp = gamma\(value*alpha);
            [sm,sn] = size(temp);
            An(1:sm,1:sn,j) = temp;
%             [si,sj,sk] = find(temp);
%             l(j) = length(si);
%             An(1:l(j),1,j) = si;
%             An(1:l(j),2,j) = sj;
%             An(1:l(j),3,j) = sk;
        elseif(N<NmaxZ)
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
    An = An(1:l,1:w,:);
    if(N==1)
        A1 = squeeze(An(:,:,:));
    elseif(N==2)
        A2 = squeeze(An(:,:,:));
    end
    betaMap = Map;
    Map = alphaMap;
    clear alphaMap

end
G00 = 1./(transpose(z)+4*t*A1);


% 
% t = 1.0;
% eta = 1e-1; %width of delta peaks
% 
% z = z+1i*eta;
% 
% 
% % A = 0.5*(-w/t+(sqrt(w/t-2).*sqrt(w/t+2)));
% A1 = zeros(1,length(z));
% A2 = zeros(2,1,length(z));
% for j=1:length(z)
%     
%     if(abs(z(j))>4)
%         N = 10;%50;
%     else
%         N = 100;%500;
%     end
%     An = sparse((N+2)/2,N/2);
%     for n = N:-1:2
%         
%         if mod(n,2) == 0
%             
%             % beta calculation
%             ind1 = zeros(n,1);
%             ind2 = zeros(n,1);
%             ind3 = -t/z(j)*ones(n,1);
%             
% %             ind3(1) = -2.0*t/z(j);
%             
%             
%             ind1(1) = 1;
%             ind2(1) = 1;
%             
%             ind1(2) = 1;
%             ind2(2) = 2;
%             
%             if (n>2)
%                 for i = 2:n/2
%                     
%                     ind1(2*i-1) = i;
%                     ind2(2*i-1) = i;
%                     
%                     ind1(2*i) = i;
%                     ind2(2*i) = i+1;
%                     
%                 end
%             end
%             ind3(n-1) = -2*t/z(j);
%             
%             beta = sparse(ind1,ind2,ind3);
%             
%             % alpha calculation
%             ind1 = zeros(n-1,1);
%             ind2 = zeros(n-1,1);
%             ind3 = -t/z(j)*ones(n-1,1);
%             ind1(1) = 1;
%             ind2(1) = 1;
% %             ind3(1) = -2*t/z(j);
% %             if(n==2)
% %                 ind3(1) = -4.0*t/z(j);
% %             end
%             
%             if (n>2)
%                 ind1(2) = 1;
%                 ind2(2) = 2;
% %                 ind3(n-2) = -2*t/z(j);
%                 
%                 for i = 2:(n/2)
%                     ind1(2*i-1) = i;
%                     ind2(2*i-1) = i;
%                     ind1(2*i) = i;
%                     ind2(2*i) = i+1;
%                 end
%                 ind1 = ind1(1:(n-1));
%                 ind2 = ind2(1:(n-1));
%             end
% 
%             
%             alpha = sparse(ind1,ind2,ind3);
%             
%             One = eye(n/2);
%             
%             
%         else
%             
%             % beta calculation
%             ind1 = zeros(n,1);
%             ind2 = zeros(n,1);
%             ind3 = -t/z(j)*ones(n,1);
%             
%             ind3(1) = -2*t/z(j);
%             ind3(n-1) = -2*t/z(j);
%                        
%             if n>2
%                 for i = 1:(n-1)/2
%                     
%                     ind1(2*i-1) = i;
%                     ind2(2*i-1) = i;
%                     
%                     ind1(2*i) = i+1;
%                     ind2(2*i) = i;
%                     
%                 end
%                 
%             end
%             ind1(n) = (n+1)/2;
%             ind2(n) = (n+1)/2;
%                        
%             beta = sparse(ind1,ind2,ind3);
%             
%             if n>1
%                 % alpha calculation
%                 ind1 = zeros(n-1,1);
%                 ind2 = zeros(n-1,1);
%                 ind3 = -t/z(j)*ones(n-1,1);
%                                
%                 ind1(1) = 1;
%                 ind2(1) = 1;
%                 
%                 ind1(2) = 2;
%                 ind2(2) = 1;
%                 ind3(1) = -2*t/z(j);
% %                 ind3(n-2) = -2*t/z(j);
%                 
%                 for i = 2:((n-1)/2)
%                     
%                     ind1(2*i-1) = i;
%                     ind2(2*i-1) = i;
%                     
%                     ind1(2*i) = i+1;
%                     ind2(2*i) = i;
%                     
%                 end
%                 alpha = sparse(ind1,ind2,ind3);
%             end
%             
%             One = eye((n+1)/2);
%             
%         end
%         
%         An = (One-beta*An)\alpha;
% 
%         if(n==3)
%             A2(:,:,j) = An;
%         end
%         
%                 
%         
%         
%         
%         
%         
%         
%         
%     end
%     A1(j) = full(An);
%     
%     
%     
% end
% 
% G00 = 1./(z+4*t*A1);
% end