clc
clear all
Ntx=4;
Nrx=8;
Kinit = min(Ntx,Nrx);
nRealiz = 50000;
vecSNIR = 0:2:40;
nSNIR = length(vecSNIR);
s02=0;
mode ='waterFilling';%  'aligning';%'uniform';%



BitErr = zeros(nSNIR,1);
symTest = [0:3];
symModTest = exp(j*(pi/2*symTest+pi/4) );
symTable = repmat(symModTest,min(Ntx,Nrx),1);
indGray=[0 1 3 2]';
nChannelWaterFilling = zeros(nRealiz,nSNIR);
typeReciver = 1;%0-ZF; 1-MSME
Q = zeros(Kinit);
rotate = 0;

symStdVec = zeros(nSNIR,1);

for iSNIR=1:nSNIR
    iSNIR
    P0 = 10^(vecSNIR(iSNIR)/10);
    cntErr = 0;
    cntBit = 0;
    symStd = zeros(nRealiz*Kinit,1);
    cntSym = 0;
    for iRealiz = 1:nRealiz
        z = sqrt(s02/2)*randn(Nrx,1)+j*sqrt(s02/2)*randn(Nrx,1);
        H = sqrt(0.5)*randn(Nrx,Ntx,1)+j*sqrt(0.5)*randn(Nrx,Ntx,1);
        [u,singValue,v]=svd(H);
        singValue = diag(singValue);
        eigValue = singValue.^2;
        eigValue = [14.5216 4.5444 0.8235 0.2344];
        singValue = sqrt(eigValue);
        %% distibution Power
        %         lamda = sqrt([9.6 4.2 1.3 0.1])';
        switch mode
            case 'uniform'
                K = Kinit;
%                 K = 2;
                p = P0/K*ones(K,1);
            case 'waterFilling'
                K = Kinit;
                while 1
                    metka=1;
                    p = zeros(K,1);
                    for iChannel=1:K
                        p(iChannel) = ( P0+typeReciver*sum(1./eigValue(1:K)) )/( singValue(iChannel)*sum(1./singValue(1:K)) )-typeReciver/eigValue(iChannel);
                        p(iChannel) = max( [p(iChannel) 0] );
                        if ~p(iChannel)
                            K = K-1;
                            metka = 0;
                            break
                        end
                    end
                    p;
                    sum(p);
                    if metka
                        break
                    end
                end
                
                
                %                 while 1
                %                     p = zeros(K,1);
                %                     C0=( P0+sum(s02./lamda(1:K).^2) )/K;
                %                     metka = 1;
                %                     for iChannel=1:K
                %                         p(iChannel) = C0-s02./lamda(iChannel).^2;
                %                         if ( p(iChannel)<0 )
                %                             K=K-1;
                %                             metka = 0;
                %                             break;
                %                         end
                %                     end
                %                     p
                % %                     C = p+1./lamda.^2
                %                     sum(p)
                %                     1;
                %                     if metka ==1
                %                         break
                %                     end
                %                     1;
                %                 end
            case 'aligning'
                K=Kinit;
                p = zeros(K,1);
                coeff = 1/sum(1./eigValue);
                for iChannel=1:K
                    p(iChannel) = P0*coeff/eigValue(iChannel);
                end
                p;
                sum(p);
                1;
        end
        sum(p);
        P = diag(p);
        
        %% QPSK
        symTx = randi(4,K,1)-1;
        symMod = exp(j*(pi/2*symTx+pi/4) );
        symTxGray = indGray(symTx+1);
        
        %         q1 = [lamda(1:K).*sqrt(p) lamda(1:K)*sqrt(P0/K)]
        %         P0*( 1/sum(1./lamda(1:K).^2)./lamda(1:K).^2 )
%         L = 4;
        if rotate==1
            Urot = 1/sqrt(K)*exp( j*2*pi/K*repmat(0:K-1,K,1).*repmat( (0:K-1)',1,K) );
        else
            Urot = eye(K);
        end
%         Urot'*Urot
        y = u(:,1:K)'*H*v(:,1:K)*sqrt(P)*Urot*symMod+u(:,1:K)'*z;
%         y = Urot'*(u(:,1:K)'*H*v(:,1:K)*sqrt(P)*Urot*symMod);
        %         u(:,1:K)'*H*v(:,1:K)*sqrt(P)*symMod+u(:,1:K)'*z;
        %         sqrt(P(1,1))*singValue(1)*symMod(1)+u(:,1)'*z;
        %         sqrt(P(2,2))*singValue(2)*symMod(2)+u(:,2)'*z;
        %         sqrt(P(3,3))*singValue(3)*symMod(3)+u(:,3)'*z
        %         sqrt(P(4,4))*singValue(4)*symMod(4)+u(:,4)'*z
        
        symEst = y./(sqrt(p).*singValue(1:K)+typeReciver );
        symStd(cntSym+1:cntSym+K) = (symEst-Urot*symMod);
        cntSym = cntSym+K;
        symEst = Urot'*symEst;
        %% demodulator
        [~,symDemod] = min(abs(repmat(symEst,1,4)-symTable(1:K,:))');
        symRxGray = indGray(symDemod);
        [nErr,q] = biterr(symRxGray,symTxGray);
        cntErr = cntErr+nErr;
        cntBit = cntBit+2*K;
        nChannelWaterFilling(iRealiz,iSNIR)=K;
        %% 
%         S = Urot*symMod;
%         Q = Q+S*S';
        
        
    end
    BitErr(iSNIR) = cntErr./cntBit;
    symStdVec(iSNIR) = mean( abs( symStd(1:cntSym) ).^2 );
end

subplot(2,4,1:2)
hold on
semilogy(vecSNIR,BitErr);
grid on
subplot(2,4,3:4)
hold on
semilogy(vecSNIR,symStdVec);
grid on
subplot(2,4,6:7)
hold on
plot(vecSNIR,symStdVec);
grid on
% y = (u'*H*v*sym)./diag(la)
% abs(v*sym)
return
figure
plot(vecSNIR,mean(nChannelWaterFilling,1))
grid on
