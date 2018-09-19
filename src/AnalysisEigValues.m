clc
clear all
nRealiz = 300000;
Ntx = 30;
Nrx = 30;
Kinit = min(Ntx,Nrx);
eigValues = zeros(Kinit*nRealiz,1);
norm = Kinit;
for iRealiz=1:nRealiz
    H = sqrt(0.5/norm)*randn(Nrx,Ntx)+j*sqrt(0.5/norm)*randn(Nrx,Ntx);
    [~,value,~]=svd(H);
    eigValues((iRealiz-1)*Kinit+1:iRealiz*Kinit) = (diag(value)).^2;
end

nSmpl = length(eigValues);
valMin = min(eigValues);
valMax = max(eigValues);
nIterv = 50;
d = (valMax-valMin)/nIterv;
vectorVal = valMin+d:d:valMax;

Prob=zeros(nIterv,1);
subplot(1,2,1)
for i = 1:nIterv
    thLeft = valMin+(i-1)*d;
    thRight = valMin+i*d;
    n = length(find( ( eigValues>thLeft).*(eigValues<=thRight)>0 ));
    Prob(i) = n;
end
normDistr = sum(Prob)*d;
Prob = Prob./normDistr;
% Prob = Prob/nSmpl;
plot(vectorVal,Prob)

da = 0.0001;
a = 0.00001:0.0001:4;
% a = vectorVal;
Ya = real(1/(2*pi)*sqrt((4-a)./a));
hold on
plot(a,Ya)

subplot(1,2,2)
for i = 1:nIterv
    cdfProb(i) = sum(Prob(1:i))*d;
end
for i = 1:length(a)
    cdfYa(i) = sum(Ya(1:i))*da;
end
plot(vectorVal,cdfProb)
hold on
plot(a,cdfYa)