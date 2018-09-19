clc
clear 
Ntx=4;
Nrx=4;
Kinit = min(Ntx,Nrx);
nRelayStation = 0;
vecSNIR = 0:2:35;
nSNIR = length(vecSNIR);
s02=0;
nRealiz = 10000;
mode ='waterFilling';%'uniform';%  'aligning';%
nSym = 4;
typeReciver = 1;%0-ZF; 1-MSME
rotate = 0;
BitErr = zeros(nSNIR,1);
symStdVec = zeros(nSNIR,1);
Pvec = zeros(nSNIR,nRelayStation+1,nRealiz);
for iSNIR = 1:nSNIR
    iSNIR
    P0 = 10^(vecSNIR(iSNIR)/10);
    symStd = zeros(nRealiz*Kinit,1);
    cntSym = 0;
    cntErr = 0;
    cntBit = 0;
    for iRealiz = 1:nRealiz
        [F,G,P, Q,K,H,U,lamda,V] = distrPower(nRelayStation,Ntx,Nrx,vecSNIR(iSNIR),typeReciver,mode);
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
            z = sqrt(s02/2)*randn(Nrx,1)+j*sqrt(s02/2)*randn(Nrx,1);
            p = diag(P(:,iRelay));
            
            coeff = sum( abs(sqrt(P(:,iRelay)).*S(:,iRelay)).^2 );
            p=p./coeff*P0;
            %sum(diag(p))
%             sum(abs(sqrt(diag(p)).*S(:,iRelay)).^2)


%             sum( abs(sqrt(P(:,iRelay)).*S(:,iRelay)).^2 )
%             P2 = P(:,iRelay);
%             P2 = P2./abs(S(:,iRelay)).^2;
%             sum( abs(sqrt(P2).*S(:,iRelay)).^2 )
%             p = diag(P2);
%             sum( abs(sqrt(p)*S(:,iRelay)).^2 )
            Pvec(iSNIR,iRelay,iRealiz) = sum( abs( sqrt(diag(p)).*S(:,iRelay) ).^2 );
            S(:,iRelay+1) = U(:,1:K,iRelay)'*H(:,:,iRelay)*V(:,1:K,iRelay)*sqrt(p)*S(:,iRelay)+U(:,1:K,iRelay)'*z;
            ind = 2;w = sqrt(p(ind,ind)*lamda(ind,iRelay))*S(ind,iRelay)+U(:,ind,iRelay)'*z;
            w = w./(sqrt( p(ind,ind)*lamda(ind,iRelay) )+typeReciver );
            S(:,iRelay+1) = S(:,iRelay+1)./(sqrt( diag(p).*lamda(1:K,iRelay) )+typeReciver );
            1;
        end
        symDeMod = Urot'*S(:,end);
        symStd(cntSym+1:cntSym+K) = (S(:,end)-S(:,1));
        cntSym = cntSym+K;
        %% demodulator
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRx,symTx);
        cntErr = cntErr+nErr;
        cntBit = cntBit+2*K;
    end
    BitErr(iSNIR) = cntErr./cntBit;
    symStdVec(iSNIR) = mean( abs( symStd(1:cntSym) ).^2 );
end
% figure
subplot(1,2,1)
hold on
semilogy(vecSNIR,BitErr);grid on;xlabel('dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
subplot(1,2,2);hold on;plot(vecSNIR,symStdVec);grid on;xlabel('dB');ylabel('MCKO');xlim([vecSNIR(1) vecSNIR(end)])
h = figure('Name','Statistcal Power')
PavVec = mean(Pvec,3);
P0 = 10.^(vecSNIR'/10);
xFig = 5;
yFig = 2;
for iRelay = 1:nRelayStation+1
    subplot(yFig,xFig,iRelay);
    plot(vecSNIR,PavVec(:,iRelay)./P0);grid on
    grid on;xlabel('dB');ylabel('average Power, P/P_0');title(strcat('Power P_',num2str(iRelay-1)) );xlim([vecSNIR(1) vecSNIR(end)])
end