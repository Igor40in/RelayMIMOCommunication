clc
clear 
% close all
Ntx=5;
Nrx=5;
Kinit = min(Ntx,Nrx);
nRelayStation = 0;
vecSNIR =0:1:12;
nSNIR = length(vecSNIR);
s02=1;
nRealiz = 30000;
mode = 'WF Separate Optimization';%'WF Independent Optimization';%'WF Independent NormSym Optimization';%'uniform';%  'aligning';%
nSym = 4;
% typeReciver = 1;%0-ZF; 1-MSME
rotate = 1;
BitErr = zeros(nSNIR,1);
symStdVec = zeros(nSNIR,1);
Pvec = zeros(nSNIR,nRelayStation+1,nRealiz);
nChannel = zeros(nRealiz,nSNIR);
for iSNIR = 1:nSNIR
    tic
    iSNIR
    P0 = 10^(vecSNIR(iSNIR)/10);
    symStd = zeros(nRealiz*Kinit,1);
    cntSym = 0;
    cntErr = 0;
    cntBit = 0;
    for iRealiz = 1:nRealiz
        [F, G, P, Q, K, H, U, lamda, V] = distrPower_v2(nRelayStation,Ntx,Nrx,vecSNIR(iSNIR),mode);%
        nChannel(iRealiz,iSNIR)=K;
        
%         symTx = [2 2 1 2 2 0 3 1 0]';
%         symTx = [2 2 1 2 2 0 3]';
%         z = [-0.578842622902202 + 0.726858584797751i;0.445391314292817 - 0.527746512540682i;0.823327699613980 - 0.00640439739626086i;1.08033888019668 - 0.684508703934331i;-0.0589423821151435 + 0.348209815641474i;-1.15559810465270 + 0.786754217735540i;0.386147838697824 + 1.13051475696330i;0.0847255638774373 - 1.52772429201155i;-0.127012169400862 - 0.663466739468551i;0.137007132622476 - 0.141770099207851i];
%         P = P*P0;
%         F = sqrt(P0)*F;
%         iRelay = 1;
%         G = inv( diag(lamda(1:K,iRelay))*diag(P(:,iRelay))+inv(Q(:,:,iRelay)) )*...
%                         sqrt(diag(P(:,iRelay)))*diag(sqrt(lamda(1:K,iRelay)) )*U(:,1:K,iRelay)';
        
        
        
        if rotate==1
            Urot = 1/sqrt(K)*exp( j*2*pi/K*repmat(0:K-1,K,1).*repmat( (0:K-1)',1,K) );
        else
            Urot = eye(K);
        end
        symTx = randi(nSym,K,1)-1;
        symMod = pskmod(symTx,nSym,pi/4,'gray');
        S = zeros(K,nRelayStation+2);
        S(:,1) = Urot*symMod;
        for iRelay = 1:nRelayStation+1
%             if  strfind(mode,'Independent') & strfind(mode,'NormSym')
%                 S(:,iRelay) = S(:,iRelay)./abs(S(:,iRelay));
%                 indNan = find(isnan(S(:,iRelay)));
%                 S(indNan,iRelay) = 0;
%                 if sum(isnan(S(:,iRelay)))
%                     1;
%                 end
%                 1;
%             elseif  strfind(mode,'Independent')
%                 Pcoeff = sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 );
%                 P(:,iRelay) = P(:,iRelay)./Pcoeff*P0;
%                 F(:,:,iRelay) = F(:,:,iRelay)*sqrt(P0/Pcoeff);
%             end
            Pvec(iSNIR,iRelay,iRealiz) = sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 );
            z = sqrt(s02/2)*randn(Nrx,1)+j*sqrt(s02/2)*randn(Nrx,1);
            S(:,iRelay+1) = G(:,:,iRelay)*H(:,:,iRelay)*F(:,:,iRelay)*S(:,iRelay)+G(:,:,iRelay)*z;
            
            if sum(isnan(S(:,iRelay)))
                1;
            end
            %             S(:,iRelay+1) = U(:,1:K,iRelay)'*H(:,:,iRelay)*V(:,1:K,iRelay)*sqrt(p)*S(:,iRelay)+U(:,1:K,iRelay)'*z;
            %             ind = 3;w = sqrt(Q(ind,ind,iRelay)*P(ind,iRelay)*lamda(ind,iRelay))*S(ind,iRelay)+U(:,ind,iRelay)'*z;
            %             w = w./(sqrt( diag(P(ind,iRelay)).*lamda(1:K,iRelay) )+typeReciver );
            
            
            %             if iRelay==1
            %                 q(:,:,iRealiz) = S(:,iRelay+1)*S(:,iRelay+1)';
            %                  sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 )
            %                  1;
            %             end
            1;
        end
        
        symDeMod = Urot'*S(:,end);
        symStd(cntSym+1:cntSym+K) = (S(:,end)-S(:,1));
        mean(abs(S(:,end)-S(:,1)).^2);
        
%         sum(1./(P.*lamda(1:K)+1));
        symStd(cntSym+1:cntSym+K)'*symStd(cntSym+1:cntSym+K);
        cntSym = cntSym+K;
        %% demodulator
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRx,symTx);
        cntErr = cntErr+nErr;
        cntBit = cntBit+2*K;
    end
    mean(abs(symStd(1:cntSym)).^2)
    BitErr(iSNIR) = cntErr./cntBit;
    symStdVec(iSNIR) = mean( abs( symStd(1:cntSym) ).^2 );
    toc
end
% figure
subplot(2,4,[1 2])
hold on
semilogy(vecSNIR,BitErr);grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
subplot(2,4,[3 4]);hold on;plot(vecSNIR,symStdVec);grid on;xlabel('SNIR, dB');ylabel('MCKO');xlim([vecSNIR(1) vecSNIR(end)])
nCh = mean(nChannel,1);
subplot(2,4,[6 7]);hold on;plot(vecSNIR,nCh);grid on;xlabel('SNIR, dB');ylabel('Average SubChannels');xlim([vecSNIR(1) vecSNIR(end)])

h = figure('Name','Statistcal Power Mean')
PavVec = mean(Pvec,3);
P0 = 10.^(vecSNIR'/10);
yFig = 2;
xFig = round((nRelayStation+1)/yFig);
for iRelay = 1:nRelayStation+1
    subplot(yFig,xFig,iRelay);
    plot(vecSNIR,PavVec(:,iRelay)./P0);grid on
    grid on;xlabel('SNIR, dB');ylabel('average Power, P/P_0');title(strcat('Power P_',num2str(iRelay-1)) );xlim([vecSNIR(1) vecSNIR(end)])
end

h = figure('Name','Statistcal Power Std')
PstdVec = std(Pvec,0,3);
for iRelay = 1:nRelayStation+1
    subplot(yFig,xFig,iRelay);
    plot(vecSNIR,PstdVec(:,iRelay)./P0);grid on
    grid on;xlabel('SNIR, dB');ylabel('Std Power, P/P_0');title(strcat('Power P_',num2str(iRelay-1)) );xlim([vecSNIR(1) vecSNIR(end)])
end
return
for iRelay=1:nRelayStation+1
    Pav(:,iRelay) = mean( Pvec(:,iRelay,:),3 )./10.^(vecSNIR'/10);
end


