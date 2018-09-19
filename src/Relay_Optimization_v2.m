clc
clear
Ntx=10;
Nrx=20;
Kinit = min(Ntx,Nrx);
nRelayStation = 4;
vecSNIR = 0:3:24;
nSNIR = length(vecSNIR);
s02=1;
nRealiz = 100;
nSym = 4;
mode = 'WF Separate Optimization';%'WF Independent Optimization';%'WF Independent NormSym Optimization';%'uniform';%  'aligning';%
optBudget = 2*1600;
rotate = 1;
P_relay = zeros(nSNIR,nRelayStation+1,nRealiz);
nChannel = zeros(nSNIR,nRealiz);
Pvec = zeros(nSNIR,nRelayStation+1,nRealiz);
K_opt = zeros(nSNIR,nRealiz);
nRepeatOptim = 2;
for iSNIR = 1:nSNIR
    tic
    iSNIR
    P0 = 10^(vecSNIR(iSNIR)/10);
    Pmin = 0*P0*ones(Kinit,1);
    Pmax = P0*ones(Kinit,1);
    Pdelta = P0/(10*Kinit)*ones(Kinit,1);%Pdelta = P0/(20*Kinit)*ones(Kinit,1);
    cntErr = 0;
    cntBit = 0;
    cntSym = 0;
    cntErr_opt = 0;
    cntBit_opt = 0;
    cntSym_opt = 0;
    for iRealiz = 1:nRealiz
        symTx = randi(nSym,Kinit,1)-1;
        symMod = pskmod(symTx,nSym,pi/4,'gray');
        Z = sqrt(s02/2)*randn(Nrx,nRelayStation+1)+j*sqrt(s02/2)*randn(Nrx,nRelayStation+1);
        
        [F, G, P, Q, K, H, U, lamda, V] = distrPower_v2(nRelayStation,Ntx,Nrx,vecSNIR(iSNIR),mode);%
        nChannel(iSNIR,iRealiz)=K;
        K1 = K;
        %         sum(1./(P.*lamda(1:K1)+1))
        if rotate==1
            Urot = 1/sqrt(K)*exp( j*2*pi/K*repmat(0:K-1,K,1).*repmat( (0:K-1)',1,K) );
        else
            Urot = eye(K);
        end
        S = zeros(K,nRelayStation+2);
        S(:,1) = Urot*symMod(1:K);
        for iRelay = 1:nRelayStation+1
            Pvec(iSNIR,iRelay,iRealiz) = sum( abs( sqrt(P(:,iRelay)).*S(:,iRelay) ).^2 );
            S(:,iRelay+1) = G(:,:,iRelay)*H(:,:,iRelay)*F(:,:,iRelay)*S(:,iRelay)+G(:,:,iRelay)*Z(:,iRelay);
            1;
        end
        %         x = zeros(Ntx,nRelayStation+1);
        %         for iRelay = 1:nRelayStation+1
        %             x(:,iRelay) = F(:,:,iRelay)*S(:,iRelay);
        %             if iRelay==1
        %                 y{iRelay} = S(:,1);
        %             else
        %                 y{iRelay} = H(:,:,iRelay-1)*x(:,iRelay-1)+Z(:,iRelay);
        %             end
        %             LAMDA(:,:,iRelay) = V(:,1:K)'*x(:,iRelay)*y{iRelay}'*inv(y{iRelay}*y{iRelay}');
        %             q1 =V(:,1:K)*LAMDA(:,:,iRelay)*y{iRelay}-x(:,iRelay);
        %             1;
        %         end
        
        
        
        
        symDeMod = Urot'*S(:,end);
        symStd(cntSym+1:cntSym+K) = (S(:,end)-S(:,1));
        cntSym = cntSym+K;
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRx,symTx(1:K));
        cntErr = cntErr+nErr;
        cntBit = cntBit+2*K;
        
%         K_opt(iSNIR,iRealiz) = Kinit;
%         K2 = Kinit;
%                 K_opt(iSNIR,iRealiz) = K;
%                 K2 = K;
%         X = repmat( (P0/K2-Pdelta(1,1)*1)*ones(K2,1),1,nRelayStation+1);
        X = P(1:K_opt(iSNIR,iRealiz),:);

%         X(1:K2,:) =  rand(K2,nRelayStation+1);
%         X = X./repmat(sum(X),K2,1)*P0;
        for iRep = 1:nRepeatOptim
            iRelay = 1;
            while iRelay <=nRelayStation+1
                while 1
                    if (iRelay == 1)&(K2~=K_opt(iSNIR,iRealiz))
                        %                     X = repmat(P0/(K_opt(iSNIR,iRealiz)*1.5)*ones(K_opt(iSNIR,iRealiz),1),1,nRelayStation+1);
                        %                     X = repmat( (P0/K_opt(iSNIR,iRealiz)-Pdelta(1,1)*0)*ones(K_opt(iSNIR,iRealiz),1),1,nRelayStation+1);
                        X = P(1:K_opt(iSNIR,iRealiz),:);
                        
%                         X(1:K2,:) =  repmat(P0/K2-Pdelta(1,1)*1,K2,nRelayStation+1);
                        
%                         X(1:K2,:) =  rand(K2,nRelayStation+1);
%                         X = X./repmat(sum(X),K2,1)*P0;
                        1;
                    end
                    %                     K = K_opt(iSNIR,iRealiz);
                    Beta = zeros(K2,1);
                    for k=1:K2
                        Beta(k,1) = 1;
                        for i=1:nRelayStation+1
                            if i==iRelay
                                continue
                            end
                            Beta(k,1) = Beta(k,1)*(lamda(k,i)*X(k,i)./(lamda(k,i)*X(k,i)+1));
                        end
                    end
                    %                     X(1:K2, iRelay ) = P0/(K2*1.0)*ones(K,1);
                    %                 X(1:K, iRelay ) = P(1:K,iRelay)*0.8;
                    %%only mimo
                    %                 X(1:K1, iRelay ) = P;
                    [metric,Parameters] = NMDS_RelayOptimization(X,iRelay,Beta,Pmin(1:K2 ),Pmax(1:K2),Pdelta(1:K2),lamda(1:K2,iRelay),P0,optBudget);
                    x_opt = Parameters(:,end);
                    
                    if length( find(x_opt<0) )%P_opt(end)<1e-2%
%                         K_opt(iSNIR,iRealiz) = K_opt(iSNIR,iRealiz)-1;
                        K2 = K2-1;
                        X(end,:) = [];
                        iRelay = 1;
                        continue
                    else
                        if x_opt(end,1)<0.001*P0
%                             K_opt(iSNIR,iRealiz) = K_opt(iSNIR,iRealiz)-1;
                            K2 = K2-1;
                            X(end,:) = [];
                            iRelay = 1;
                            continue
                        end
                    end
                    K_opt(iSNIR,iRealiz) = K2;
                    %                 sum(x_opt)
                    X(1:K2,iRelay) = x_opt;
                    break;
                end
                iRelay = iRelay+1;
                1;
            end
        end
        P_opt = zeros(K2,nRelayStation+1);
        for iRelay = 1:nRelayStation+1
            if iRelay==1
                P_opt(:,iRelay) = X(:,iRelay);
            else
                P_opt(:,iRelay) = X(:,iRelay)./(lamda(1:K2,iRelay-1).*X(:,iRelay-1)+1);
            end
        end
%         K2 = K;
        %         sum(1./(P(:,1).*lamda(1:K1)'+1))
        %         sum(1./(X(:,1).*lamda(1:K2)'+1))
        %         if ( sum(1./(P_opt.*lamda(1:K2)+1))-sum(1./(P.*lamda(1:K1)+1)) )<-0.01
        %             (sum(1./(P_opt.*lamda(1:K2)+1))-sum(1./(P.*lamda(1:K1)+1)) )/sum(1./(P_opt.*lamda(1:K2)+1));
        %             [sum(P) sum(P_opt)];
        %             1;
        %         end
        %         sum(1./(P_opt.*lamda(1:K2)+1))
        %         symTx = randi(nSym,K,1)-1;
        %         symMod = pskmod(symTx,nSym,pi/4,'gray');
        if rotate==1
            Urot = 1/sqrt(K2)*exp( j*2*pi/K*repmat(0:K2-1,K2,1).*repmat( (0:K2-1)',1,K2) );
        else
            Urot = eye(K2);
        end
        S = zeros(K2,2);
        S(:,1) = Urot*symMod(1:K2);
        
        %         Z = sqrt(s02/2)*randn(Nrx,nRelayStation+1)+j*sqrt(s02/2)*randn(Nrx,nRelayStation+1);
        for iRelay = 1:nRelayStation+1
            if iRelay==1
                FF{1} = V(:,1:K2,iRelay)*sqrt(diag(P_opt(1:K2,iRelay)) );
                x(:,1) =FF{1}*S(:,1) ;
            else
                FF{iRelay} = V(:,1:K2,iRelay)*sqrt(diag(P_opt(1:K2,iRelay)))*U(:,1:K2,iRelay-1)';
                x(:,iRelay) = FF{iRelay}*y(:,iRelay-1);
            end
            y(:,iRelay) = H(:,:,iRelay)*x(:,iRelay)+Z(:,iRelay);
        end
        
        %%Hcommon
        %         Hcommon = 1;
        %         for iRelay = nRelayStation+1:-1:1
        %             Hcommon = Hcommon*H(:,:,iRelay)*FF{iRelay};
        %         end
        %         Hcommon2 = U(:,1:K,end)*diag(prod(sqrt(P_opt.*lamda(1:K,:))'));
        Hcommon2 = U(:,1:K2,end)*diag(prod(sqrt(P_opt.*lamda(1:K2,:)),2));
        %         Hcommon2 = Hcommon;
        %         Hcommon2 = U(:,1:K,end)*diag(sqrt(P_opt.*lamda(1:K,:)));
        
        noiseCommmon = 0;
        for l=2:nRelayStation+1
            mult = 1;
            for i=nRelayStation+1:-1:l
                mult = mult*H(:,:,i)*FF{i};
            end
            noiseCommmon = noiseCommmon+mult*Z(:,l-1);
        end
        noiseCommmon = noiseCommmon+Z(:,end);
        1;
        yest = Hcommon2*symMod(1:K2)+noiseCommmon;
        
        Cv = 0;
        for l=2:nRelayStation+1
            mult1 = 1;
            for i=nRelayStation+1:-1:l
                mult1 = mult1*H(:,:,i)*FF{i};
            end
            
            mult2 = 1;
            for i=l:nRelayStation+1
                mult2 = mult2*FF{i}'*H(:,:,i)';
            end
            Cv = Cv+mult1*mult2;
        end
        Cv = Cv+eye(Nrx);
        1;
        W = (Hcommon2*Hcommon2'+Cv)\Hcommon2;
        S(:,2) = (W'*y(:,end));
        symDeMod = Urot'*S(:,2);
        symRx = pskdemod(symDeMod,nSym,pi/4,'gray');
        [nErr,q] = biterr(symRx,symTx(1:K2));
        cntErr_opt = cntErr_opt+nErr;
        cntBit_opt = cntBit_opt+2*K2;
        
        symStd_opt(cntSym_opt+1:cntSym_opt+K2) = (S(:,end)-S(:,1));
        cntSym_opt = cntSym_opt+K2;
        1;
        %         P_relay(iSNIR,:,iRealiz) = diag(x'*x);
        %         x(:,2)
        %         x2 = V(:,1:K,2)*diag(sqrt(P_opt(:,2)))*diag(sqrt(lamda(1:K,1)))*diag(sqrt(P_opt(:,1)))*symMod
        %         x22 = V(:,1:K,2)*diag(sqrt(P_opt(:,2)))*U(:,1:K,1)'*U(:,1:K,1)*diag(sqrt(lamda(1:K,1)))*V(:,1:K,1)'*V(:,1:K,1)*diag(sqrt(P_opt(:,1)))*symMod
        %           x222 = diag(P_opt(:,2))*diag(P_opt(:,1))*diag(lamda(1:K,1))
        
        %%check p_i^2
        %         for iRelay=1:nRelayStation+1
        %             if iRelay==1
        %                 xx(:,:,iRelay) = diag(P_opt(:,iRelay));
        %             else
        %                 xx(:,:,iRelay) = diag(P_opt(:,iRelay)).*(xx(:,:,iRelay-1).*diag(lamda(1:K,iRelay-1))+eye(K)) ;
        %             end
        %                 sum(diag(xx(:,:,iRelay)))
        %         end
        %         sum(xx)
        1;
        
        
        
        
        
        continue;
        if  K_opt(iSNIR,iRealiz) == nChannel(iSNIR,iRealiz)
            for k=1:K_opt(iSNIR,iRealiz)
                value1 = 1;
                for iRelay = 1:nRelayStation+1
                    value1 = value1*P(k,iRelay)*lamda(k,iRelay);
                end
                
                value2 = 0;
                for iRelay = 2:nRelayStation+1
                    value3 = 1;
                    for i=iRelay:nRelayStation+1
                        value3 = value3*P(k,i)*lamda(k,i)
                    end
                    value2 = value2 + value3;
                end
                E(k) = (1+value1/(1+value2) )^(-1);
            end
        end
        
        1;
    end
    BitErr(iSNIR) = cntErr./cntBit;
    BitErr_opt(iSNIR) = cntErr_opt./cntBit_opt;
    symStdVec(iSNIR) = mean( abs( symStd(1:cntSym) ).^2 );
    symStdVec_opt(iSNIR) = mean( abs( symStd_opt(1:cntSym_opt) ).^2 );
    toc
end
h = subplot(2,4,[1 2]);
semilogy(vecSNIR,BitErr,'-','LineWidth',1);grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
hold on
plot(vecSNIR,BitErr_opt,'--','LineWidth',1)%;grid on;xlabel('SNIR, dB');ylabel('BER');xlim([vecSNIR(1) vecSNIR(end)])
title(h,strcat('nRelayStation=',num2str(nRelayStation),'[',num2str(Ntx),'x',num2str(Nrx),']'))
set(h,'FontSize',12)
legend(h,{'раздельная оптимизация','совместная оптимизация'})

h = subplot(2,4,[3 4]);
plot(vecSNIR,symStdVec,'-','LineWidth',1)
hold on
plot(vecSNIR,symStdVec_opt,'--','LineWidth',1)
xlabel('SNIR, dB');ylabel('CKO');grid on
% legend('WF','Optimization')
set(h,'FontSize',12)

h = subplot(2,4,[6 7]);
Kav = mean(nChannel,2);
K_optAv = mean(K_opt,2);
plot(vecSNIR,Kav,'-','LineWidth',1)
hold on
plot(vecSNIR,K_optAv,'--','LineWidth',1)
xlabel('SNIR, dB');ylabel('<Number Channels>');grid on
set(h,'FontSize',12)