clc
clear
close all
Ntx=10;
Nrx=20;
Kinit = min(Ntx,Nrx);
nRelayStation = 2;
vecSNIR = 0:2:16;
nSNIR = length(vecSNIR);
s02=1;
nRealiz = 3000;
nSym = 4;
mode1 = 'WF Separate Optimization';%'WF Independent Optimization';%'WF Independent NormSym Optimization';%'uniform';%  'aligning';%
mode2 = 'WF Independent Optimization';%'WF Separate Optimization';%
% optBudget = 2*1600;
rotate = 1;
% P_relay = zeros(nSNIR,nRelayStation+1,nRealiz);
nChannelSep = zeros(nSNIR,nRealiz);
nChannelIndep = zeros(nSNIR,nRealiz);
PvecSep = zeros(nSNIR,nRelayStation+1,nRealiz);
PvecIndep = zeros(nSNIR,nRelayStation+1,nRealiz);
K_opt = zeros(nSNIR,nRealiz);
% nRepeatOptim = 3;
norm = 1;
for iSNIR = 1:nSNIR
    tic
    iSNIR
    P0 = 10^(vecSNIR(iSNIR)/10);
    %     Pmin = 0*P0*ones(Kinit,1);
    %     Pmax = P0*ones(Kinit,1);
    %     Pdelta = P0/(10*Kinit)*ones(Kinit,1);%Pdelta = P0/(20*Kinit)*ones(Kinit,1);
    cntErrSep = 0;
    cntBitSep = 0;
    cntSymSep = 0;
    cntErrInd = 0;
    cntBitInd = 0;
    cntSymInd = 0;
    
    cntErrIndNorm = 0;
    cntBitIndNorm = 0;
    cntSymIndNorm = 0;
    for iRealiz = 1:nRealiz
        symTx = randi(nSym,Kinit,1)-1;
        symMod = pskmod(symTx,nSym,pi/4,'gray');
        
        H = sqrt(0.5/norm)*randn(Nrx,Ntx,nRelayStation+1)+j*sqrt(0.5/norm)*randn(Nrx,Ntx,nRelayStation+1);
        Z = sqrt(s02/2)*randn(Nrx,nRelayStation+1)+j*sqrt(s02/2)*randn(Nrx,nRelayStation+1);
        
%         [F, G, P, Q, K, H, U, lamda, V] = distrPower_v3(H,vecSNIR(iSNIR),mode1);%
%         q1 = F;
%         nChannelSep(iSNIR,iRealiz)=K;
%         if rotate==1
%             Urot = 1/sqrt(K)*exp( j*2*pi/K*repmat(0:K-1,K,1).*repmat( (0:K-1)',1,K) );
%         else
%             Urot = eye(K);
%         end
%         S = zeros(K,nRelayStation+2);
%         S(:,1) = Urot*symMod(1:K);
%         for iRelay = 1:nRelayStation+1
%             PvecSep(iSNIR,iRelay,iRealiz) = sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 );
%             S(:,iRelay+1) = G(:,:,iRelay)*H(:,:,iRelay)*F(:,:,iRelay)*S(:,iRelay)+G(:,:,iRelay)*Z(:,iRelay);
%             1;
%         end
%         symDeMod = Urot'*S(:,end);
%         symStdSep(cntSymSep+1:cntSymSep+K) = (S(:,end)-S(:,1));
%         cntSymSep = cntSymSep+K;
%         symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
%         [nErr,q] = biterr(symRx,symTx(1:K));
%         cntErrSep = cntErrSep+nErr;
%         cntBitSep = cntBitSep+2*K;
        
        [F, G, P, Q, K, H, U, lamda, V] = distrPower_v3(H,vecSNIR(iSNIR),mode2);%
        q2 = F;
        nChannelIndep(iSNIR,iRealiz)=K;
        if rotate==1
            Urot = 1/sqrt(K)*exp( j*2*pi/K*repmat(0:K-1,K,1).*repmat( (0:K-1)',1,K) );
        else
            Urot = eye(K);
        end
        S = zeros(K,nRelayStation+2);
        S(:,1) = Urot*symMod(1:K);
        Snorm = zeros(K,nRelayStation+2);
        Snorm(:,1) = Urot*symMod(1:K);
        for iRelay = 1:nRelayStation+1
            PvecIndep(iSNIR,iRelay,iRealiz) = sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 );
            coeff = P0/PvecIndep(iSNIR,iRelay,iRealiz);
            PvecIndepNorm(iSNIR,iRelay,iRealiz) = PvecIndep(iSNIR,iRelay,iRealiz)*coeff;

            S(:,iRelay+1) = G(:,:,iRelay)*H(:,:,iRelay)*F(:,:,iRelay)*S(:,iRelay)+G(:,:,iRelay)*Z(:,iRelay);
            Snorm(:,iRelay+1) = G(:,:,iRelay)*H(:,:,iRelay)*F(:,:,iRelay)*Snorm(:,iRelay)*coeff+G(:,:,iRelay)*Z(:,iRelay);
            1;
        end
        symDeMod = Urot'*S(:,end);
        symStdIndep(cntSymInd+1:cntSymInd+K) = (S(:,end)-S(:,1));
        cntSymInd= cntSymInd+K;
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRx,symTx(1:K));
        cntErrInd = cntErrInd+nErr;
        cntBitInd = cntBitInd+2*K;
        
        symDeModnorm = Urot'*Snorm(:,end);
        symStdIndepNorm(cntSymIndNorm+1:cntSymIndNorm+K) = (Snorm(:,end)-Snorm(:,1));
        cntSymIndNorm= cntSymIndNorm+K;
        symRxNorm = pskdemod(symDeModnorm,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRxNorm,symTx(1:K));
        cntErrIndNorm = cntErrIndNorm+nErr;
        cntBitIndNorm = cntBitIndNorm+2*K;

    end
%     BitErrSep(iSNIR) = cntErrSep./cntBitSep;
%     symStdVecSep(iSNIR) = mean( abs( symStdSep(1:cntSymSep) ).^2 );
    BitErrIndep(iSNIR) = cntErrInd./cntBitInd;
    BitErrIndepNorm(iSNIR) = cntErrIndNorm./cntBitIndNorm;

    symStdVecInd(iSNIR) = mean( abs( symStdIndep(1:cntSymInd) ).^2 );
    symStdVecIndNorm(iSNIR) = mean( abs( symStdIndepNorm(1:cntSymIndNorm) ).^2 );
    toc
end
h = subplot(2,4,[1 2]);
% semilogy(vecSNIR,BitErrSep,'-','LineWidth',1);grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
hold on
y1 = plot(vecSNIR,BitErrIndep,'--','LineWidth',1)%;grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
y2 = plot(vecSNIR,BitErrIndepNorm,'--','LineWidth',1)%;grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
title(h,strcat('nRelayStation=',num2str(nRelayStation),'[',num2str(Ntx),'x',num2str(Nrx),']'))
set(h,'FontSize',12)
% legend(h,{'независимая оптимизация','независимая оптимизация\n(нормализация)'})
set(y1,'DisplayName','независимая оптимизация');
set(y2,'DisplayName',['независимая оптимизация',sprintf('\n'),'(нормализация)']);
legend(h,'show');

h = subplot(2,4,[3 4]);
% plot(vecSNIR,symStdVecSep,'-','LineWidth',1)
hold on
plot(vecSNIR,symStdVecInd,'--','LineWidth',1)
plot(vecSNIR,symStdVecIndNorm,'--','LineWidth',1)
xlabel('SNIR, dB');ylabel('CKO');grid on
set(h,'FontSize',12)

h = subplot(2,4,[6 7]);
% Ksep = mean(nChannelSep,2);
Kind = mean(nChannelIndep,2);
% plot(vecSNIR,Ksep,'-','LineWidth',1)
hold on
plot(vecSNIR,Kind,'--','LineWidth',1)
xlabel('SNIR, dB');ylabel('<Number Channels>');grid on
set(h,'FontSize',12)

h = figure('Name','Statistcal Power Mean')
% PavVecSep = mean(PvecSep,3);
PavVecIndep = mean(PvecIndep,3);
PavVecIndepNorm = mean(PvecIndepNorm,3);
P0 = 10.^(vecSNIR'/10);
yFig = 2;
xFig = round((nRelayStation+1)/yFig);
for iRelay = 1:nRelayStation+1
    hh = subplot(yFig,xFig,iRelay);
%     plot(vecSNIR,PavVecSep(:,iRelay)./P0);grid on;hold on
    hold on;
    y1 = plot(vecSNIR,PavVecIndep(:,iRelay)./P0)
    y2 = plot(vecSNIR,PavVecIndepNorm(:,iRelay)./P0)
    xlabel('SNIR, dB');ylabel('average Power, P/P_0');title(strcat('Power P_',num2str(iRelay-1)) );xlim([vecSNIR(1) vecSNIR(end)])
    if iRelay==1
%         legend({'независимая оптимизация',strcat('независимая оптимизация',sprintf('\n'),'(нормализация)')})
    set(y1,'DisplayName','независимая оптимизация');
    set(y2,'DisplayName',['независимая оптимизация',sprintf('\n'),'(нормализация)']);
    legend(hh,'show');
    end
end

h = figure('Name','Statistcal Power Std')
% PstdVecSep = std(PvecSep,0,3);
PstdVecIndep = std(PvecIndep,0,3);
PstdVecIndepNorm = std(PvecIndepNorm,0,3);
for iRelay = 1:nRelayStation+1
    hh = subplot(yFig,xFig,iRelay);
%     plot(vecSNIR,PstdVecSep(:,iRelay)./P0);grid on;
    hold on
    y1 = plot(vecSNIR,PstdVecIndep(:,iRelay)./P0)
    y2 = plot(vecSNIR,PstdVecIndepNorm(:,iRelay)./P0)
    xlabel('SNIR, dB');ylabel('Std Power, P/P_0');title(strcat('Power P_',num2str(iRelay-1)) );xlim([vecSNIR(1) vecSNIR(end)])
    if iRelay==1
%         legend({'независимая оптимизация','независимая оптимизация\n(нормализация)'})
    set(y1,'DisplayName','независимая оптимизация');
    set(y2,'DisplayName',['независимая оптимизация',sprintf('\n'),'(нормализация)']);
    legend(hh,'show');
    end
end