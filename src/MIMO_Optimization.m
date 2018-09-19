clc
clear
Ntx=15;
Nrx=15;
Kinit = min(Ntx,Nrx);
nRelayStation = 0;
vecSNIR = 0:2:20;
nSNIR = length(vecSNIR);
s02=1;
nRealiz = 200;
mode = 'WF Separate Optimization';%'WF Independent Optimization';%'WF Independent NormSym Optimization';%'uniform';%  'aligning';%
nSym = 4;
rotate = 1;
BitErr = zeros(nSNIR,1);
BitErr_Opt = zeros(nSNIR,1);
symStdVec = zeros(nSNIR,1);
symStdVec_Opt = zeros(nSNIR,1);
Pvec = zeros(nSNIR,nRelayStation+1,nRealiz);
Pvec_Opt = zeros(nSNIR,nRealiz);

nChannel = zeros(nSNIR,nRealiz);
K_opt = zeros(nSNIR,nRealiz);
MMSE = zeros(nSNIR,nRealiz);
MMSE_Opt = zeros(nSNIR,nRealiz);
optBudget = 300;
for iSNIR = 1:nSNIR
    tic
    iSNIR
    P0 = 10^(vecSNIR(iSNIR)/10);
    symStd = zeros(nRealiz*Kinit,1);
    symStd_Opt = zeros(nRealiz*Kinit,1);
    cntSym = 0;
    cntSym_Opt = 0;
    cntErr = 0;
    cntErr_Opt = 0;
    cntBit = 0;
    cntBit_Opt = 0;
    Pmin = 0*P0*ones(Kinit,1);
    Pmax = P0*ones(Kinit,1);
    Pdelta = P0/10*ones(Kinit,1);
    for iRealiz = 1:nRealiz
        [F, G, P, Q, K, H, U, lamda, V] = distrPower_v2(nRelayStation,Ntx,Nrx,vecSNIR(iSNIR),mode);%
        nChannel(iSNIR,iRealiz)=K;
        %         MMSE(iSNIR,iRealiz)  = sum(1./(P(1:K).*lamda(1:K)+1));
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
            Pvec(iSNIR,iRelay,iRealiz) = sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 );
            z = sqrt(s02/2)*randn(Nrx,1)+j*sqrt(s02/2)*randn(Nrx,1);
            S(:,iRelay+1) = G(:,:,iRelay)*H(:,:,iRelay)*F(:,:,iRelay)*S(:,iRelay)+G(:,:,iRelay)*z;
            1;
        end
        symDeMod = Urot'*S(:,end);
        symStd(cntSym+1:cntSym+K) = (S(:,end)-S(:,1));
%         sum(1./(P.*lamda(1:K)+1))
%         symStd(cntSym+1:cntSym+K)'*symStd(cntSym+1:cntSym+K)
        cntSym = cntSym+K;
        %% demodulator
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRx,symTx);
        cntErr = cntErr+nErr;
        cntBit = cntBit+2*K;
        
        % NMDS Optimization
        K_opt(iSNIR,iRealiz) = Kinit;
        while 1
            k = K_opt(iSNIR,iRealiz);
            Pinit = P0/(Kinit*1.5)*ones(k,1);
            [metric,Parameters] = NMDS_Optimization(Pinit,Pmin(1:k ),Pmax(1:k),Pdelta(1:k),lamda(1:k),P0,optBudget);
%             return
            P_opt = Parameters(:,end);
            sum(1./(P_opt.*lamda(1:k) +1));
            if length( find(P_opt<0) )%P_opt(end)<1e-2%
                K_opt(iSNIR,iRealiz) = K_opt(iSNIR,iRealiz)-1;
                continue
            else
                if P_opt(end,1)<0.001*P0
                    K_opt(iSNIR,iRealiz) = K_opt(iSNIR,iRealiz)-1;
                    continue
                end
            end
            %             MMSE_Opt(iSNIR,iRealiz) = sum( 1./(P_opt(1:k).*lamda(1:k)+1) );
            break
        end
        
        if rotate==1
            Urot = 1/sqrt(k)*exp( j*2*pi/K*repmat(0:k-1,k,1).*repmat( (0:k-1)',1,k) );
        else
            Urot = eye(k);
        end
        symTx = randi(nSym,k,1)-1;
        symMod = pskmod(symTx,nSym,pi/4,'gray');
        S_Opt = zeros(k,nRelayStation+2);
        S_Opt(:,1) = Urot*symMod;
        Pvec_Opt(iSNIR,iRealiz) = sum( abs( sqrt(P_opt).*S_Opt(:,1) ).^2 );
%         z = sqrt(s02/2)*randn(Nrx,1)+j*sqrt(s02/2)*randn(Nrx,1);
        G_Opt = inv( diag(lamda(1:k,1))*diag(P_opt)+eye(k) )*...
            sqrt(diag(P_opt))*diag(sqrt(lamda(1:k)))*U(:,1:k)';
        F_Opt = V(:,1:k,1)*sqrt( diag(P_opt) );
        S_Opt(:,iRelay+1) = G_Opt(:,:,1)*H(:,:,1)*F_Opt(:,:,1)*S_Opt(:,1)+G_Opt(:,:,1)*z;

        symDeMod = Urot'*S_Opt(:,end);
        symStd_Opt(cntSym_Opt+1:cntSym_Opt+k) = (S_Opt(:,end)-S_Opt(:,1));
        cntSym_Opt = cntSym_Opt+k;
        %% demodulator
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr_Opt,q] = biterr(symRx,symTx);
        cntErr_Opt = cntErr_Opt+nErr_Opt;
        cntBit_Opt = cntBit_Opt+2*k;
        if nErr~=nErr_Opt
            1;
        end
        1;
    end
    BitErr(iSNIR) = cntErr./cntBit;
    symStdVec(iSNIR) = mean( abs( symStd(1:cntSym) ).^2 );
    BitErr_Opt(iSNIR) = cntErr_Opt./cntBit_Opt;
    symStdVec_Opt(iSNIR) = mean( abs( symStd(1:cntSym_Opt) ).^2 );
    toc
end
% MMSE = mean(MMSE,2);
% MMSE_Opt = mean(MMSE_Opt,2);
% MMSE_delta = abs(MMSE-MMSE_Opt);
% MMSE_delta = mean(MMSE_delta,2);
h = subplot(2,4,[1 2]);
semilogy(vecSNIR,BitErr,'LineWidth',2);
hold on
semilogy(vecSNIR,BitErr_Opt,'--','LineWidth',2);grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
title(h,strcat('MIMO=',num2str(Ntx),'x',num2str(Nrx)))
set(h,'FontSize',12)

h = subplot(2,4,[3 4]);
plot(vecSNIR,symStdVec,'LineWidth',2)
hold on
plot(vecSNIR,symStdVec_Opt,'--','LineWidth',2)
xlabel('SNIR, dB');ylabel('CKO');grid on
legend('WF','Optimization')
set(h,'FontSize',12)

h = subplot(2,4,[6 7]);
nCh = mean(nChannel,2);
K_opt = mean(K_opt,2);
plot(vecSNIR,nCh,'LineWidth',2)
hold on
plot(vecSNIR,K_opt,'--','LineWidth',2)
xlabel('SNIR, dB');ylabel('<Number Channels>');grid on
set(h,'FontSize',12)
return


nChannel = mean(nChannel,2);
K_opt = mean(K_opt,2);

set(h,'FontSize',12)
h = subplot(2,4,[3 4]);
plot(vecSNIR,MMSE_delta,'LineWidth',2)
xlabel('dB');ylabel('\Delta CKO');grid on
set(h,'FontSize',12)
h = subplot(2,4,[6 7]);
plot(vecSNIR,nChannel,'LineWidth',2)
hold on
plot(vecSNIR,K_opt,'--','LineWidth',2)
xlabel('dB');ylabel('<Number Channel>');grid on
set(h,'FontSize',12)

