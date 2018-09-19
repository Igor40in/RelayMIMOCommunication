function [F,G,P, Q,K,H,U,eigValue,V]=distrPower(nRelay,Ntx,Nrx,SNIR,typeReciver,mode)
Kinit = min(Ntx,Nrx);
P0 = 10^(SNIR/10);
H = sqrt(0.5)*randn(Nrx,Ntx,nRelay+1)+j*sqrt(0.5)*randn(Nrx,Ntx,nRelay+1);
% H = [0.372300254368754 - 0.00675457414572699i -0.0367628650280381 + 0.892245418745627i -0.760876958863598 - 1.10726349704768i -1.72388208668673 + 0.706074776888852i;-0.986952719841068 + 0.334433704450663i 0.936619451754528 + 1.08016086111540i 0.491492193445486 - 0.775955169208214i 0.923755293462155 + 0.477690264955616i;-0.160688556853085 + 0.400735711575352i 0.882514874400525 - 0.234131559555941i -0.772839422224093 - 0.0288852440663898i -0.179028213721407 + 0.0184294193752975i;0.272926047160276 + 1.60962553611882i 1.06555596956601 + 0.640158647576080i -0.183052936742923 + 0.750211819023124i 0.0463113922983262 - 0.609176575902673i]
% H = [0.347065911749090 + 0.448680138520898i -0.0923433294174835 - 0.844227320738246i 1.92256138294173 - 0.894893435781267i -1.43551371311988 - 0.944779630607836i 0.649429197355367 + 0.677401121007361i 0.101698878614314 + 0.747140022404322i -0.0493776348095990 + 0.0978266745662167i 0.550792223699840 + 0.580560405734954i -0.573092981468529 + 0.902153492324553i -0.431568250984491 + 0.141286964208784i;-0.495257923257252 + 0.815513879029497i -0.528475145670871 - 1.13488078085302i 0.00673857699002712 + 0.425710372668068i 0.0278459534934303 - 0.366279402374708i 0.388369336796214 + 0.266570024727337i 0.144750502134817 + 0.200941050892603i 0.670263481830292 - 0.110945503263869i 1.48279637885568 + 0.211282464737902i 0.167261709390836 + 1.30175755251161i 0.581555118342853 - 0.0998641812473246i;1.18313791622376 + 0.837625709183349i 0.576971282964319 + 0.598279040835383i -0.287757151125083 + 1.13449710416679i 0.458366464458715 - 0.585967022621084i 0.167891482772170 + 0.491885010582092i 1.06622652422905 - 1.19353865138968i 0.208571988101931 + 1.19422099431444i 0.509829598534678 - 0.364548700510952i -0.385071316375026 + 0.251995468159808i 1.10649372667141 + 0.113404038485887i;-0.460947640059068 + 0.261646367543188i 1.13664499050747 + 0.456405579528999i -0.196282325516195 + 0.443154843079336i 0.341519417866849 - 0.335114969920296i 0.854022178458493 - 0.381294339173429i 0.860336261140855 - 0.268135898630528i -1.67561969052154 + 0.229967877699710i -0.637317935878853 + 0.700634678038804i 0.209343241452974 - 0.0722198914199157i 0.650350031825359 + 0.140721973829524i;1.03791572908060 - 0.692393647982220i -1.64989190824666 + 0.192851643969447i -0.955671978911341 - 0.576663980155574i -1.54191320761167 - 0.144886810525326i 0.219349249891475 - 0.958078261403366i 0.736235255664792 - 0.328209719901179i 0.0695022819502740 - 0.792855836387654i 0.789607311286492 + 1.19554847953884i -0.0529107440849235 - 0.209602202195590i -1.35526992144490 + 0.612344105039080i;-1.68024660934791 - 0.478222791339813i 0.164191036340043 + 0.568825341005564i 0.601779264125291 + 0.207614280980248i -0.661620789629777 + 1.16330497969824i 0.704802431974308 + 0.655266367031554i 0.374308393048227 - 0.180414885061811i 0.746077975727588 + 0.0826177976372925i 0.0412972731913097 - 0.497602770951078i 2.02924995997643 - 0.940257076736530i 0.514841314998012 + 0.532581155258042i;-0.211424953342217 + 0.185962114010602i 0.651839374225012 - 1.29972255710074i -0.141439011010883 - 0.826271289687212i 1.78802892621144 + 0.567347039854983i -0.305873070806384 + 0.194486853706926i -0.299885834242861 + 0.277915004003000i 1.63019379598899 + 0.0391154425749193i -0.197233685168507 + 0.177778224200280i 1.89341528498122 + 0.925664575995708i -0.197137547711058 - 1.27016053449527i;0.366879066554910 - 2.50094085211860i 0.696462178017664 + 0.130173608063395i -0.189434560414857 + 0.199190913956420i 1.80145294053095 + 0.0777333099277575i -0.625629950451812 - 0.647765130898620i 1.20872973401845 - 0.650781152776072i -0.514032396912790 - 1.18857418500625i 0.910092955689122 + 0.485066719426446i 0.333274507140209 + 0.179931580698528i -0.389585439142311 - 0.254033839999090i;-0.148120072206945 - 0.0268997372057953i -0.296875063739360 - 1.38789705957577i -1.52644531814743 - 0.0506913111230935i 0.269413352174774 - 0.212182493452369i 0.309511517930694 + 0.00476865139037810i 0.926860485108819 + 0.180628045642815i -0.276962261344992 + 0.110247234055784i 0.179848617710254 + 0.840460924818070i 0.585704791683789 + 0.785858594571912i 0.179308054905278 - 1.80558264268594i;0.953521145878068 + 0.161113165307411i -0.857162450463999 + 0.575738698525164i -0.664327352780252 + 0.0980024021617333i -0.186115083518577 - 0.657323425786544i -1.27457745613167 + 0.511193187114342i -0.0924152709558573 - 0.235425685808785i 0.366788370091362 - 0.690358386453712i 0.780457281799533 + 0.258478513523699i -0.747780139624982 + 0.383752590766452i 1.18056443162047 - 0.0137711370248039i]
% H = repmat(H,1,1,3);
singValue = zeros(Kinit,nRelay+1);
eigValue = zeros(Kinit,nRelay+1);
U = zeros(Nrx,Nrx,nRelay+1);
V = zeros(Ntx,Ntx,nRelay+1);
for iRelay = 1:nRelay+1
    [u,value,v]=svd(H(:,:,iRelay));
    U(:,:,iRelay) = u;
    V(:,:,iRelay) = v;
    singValue(:,iRelay) = diag(value);
    eigValue(:,iRelay) = singValue(:,iRelay).^2;
end
% eigValue(:,1) = [9.216509064996606;6.238802890683078;2.457554530763014;0.184353637022665];
% singValue(:,1) = sqrt(eigValue(:,1));

% eigValue = [14.5216 4.5444 0.8235 0.2344]';
% singValue = sqrt(eigValue);
K = Kinit;
isReadyChannels = 0;
isBadChannel = 0;
while ~isReadyChannels
    P = zeros(K,nRelay+1);
    iRelay = 1;
    while iRelay<=nRelay+1
        
        
        switch mode
            case 'uniform'
                P(:,iRelay) = P0/K*ones(K,1);
                iRelay = iRelay+1;
            case 'waterFilling'
                while 1
                    if iRelay == 1
                        Q(:,:,iRelay) = eye(K);
                    else
%                         F = V(:,1:K,iRelay-1)*sqrt( diag(P(1:K,iRelay-1)) );
                        Q(:,:,iRelay) = G(:,:,iRelay-1)*H(:,:,iRelay-1)*F(:,:,iRelay-1)*Q(:,:,iRelay-1)*F(:,:,iRelay-1)'*H(:,:,iRelay-1)'*G(:,:,iRelay-1)'+G(:,:,iRelay-1)*G(:,:,iRelay-1)';
                    end
                    P(1:K,iRelay) = zeros(K,1);
                    isBadChannel = 0;
                    for iChannel=1:K
                        P(iChannel,iRelay) = ( P0+sum(1./eigValue(1:K,iRelay)) )/( sqrt(Q(iChannel,iChannel,iRelay))*singValue(iChannel,iRelay)*sum( sqrt(diag(Q(:,:,iRelay)))./singValue(1:K,iRelay)) )-...
                            1/( Q(iChannel,iChannel,iRelay)*eigValue(iChannel,iRelay) );
                        P(iChannel,iRelay) = real(P(iChannel,iRelay));
                        P(iChannel,iRelay) = max( [P(iChannel,iRelay) 0] );
                        
                        if ~P(iChannel,iRelay)
                            K = K-1;
                            isBadChannel = 1;
                            break
                        end
                    end
                    if ~isBadChannel %iChannel == K
                        F(:,:,iRelay) = V(:,1:K,iRelay)*sqrt( diag(P(1:K,iRelay)) );
                        G(:,:,iRelay) = inv( diag(eigValue(1:K,iRelay))*diag(P(:,iRelay))+inv(Q(:,:,iRelay)) )*...
                                        sqrt(diag(P(:,iRelay)))*diag(singValue(1:K,iRelay))*U(:,1:K,iRelay)';
%                         sum(diag(F(:,:,iRelay)*Q(:,:,iRelay)*F(:,:,iRelay)'))
                        iRelay = iRelay+1;
                        break;
                    else
                        iRelay = 1;
                        Q = [];
                        G = [];
                        F = [];
                        P = zeros(K,nRelay+1);
                    end
                end
                %                 if ~isBadChannel
                %                     break;
                %                 end
            case 'aligning' %% not chacking
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
        if iRelay==nRelay+2
            isReadyChannels = 1;
            break
        end
    end
end
P(K+1:end,:) = [];
return

for iRelay = 1:nRelay
    
end
