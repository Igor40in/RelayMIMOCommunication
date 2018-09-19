function [F,G,P, Q,K,H,U,eigValue,V]=distrPower(nRelay,Ntx,Nrx,SNIR,mode)
F=[];G=[];P=[];Q=[];K=[];H=[];U=[];eigValue=[];V=[];

Kinit = min(Ntx,Nrx);
P0 = 10^(SNIR/10);
% P0 = 5;
norm = 1;%Kinit;%
H = sqrt(0.5/norm)*randn(Nrx,Ntx,nRelay+1)+j*sqrt(0.5/norm)*randn(Nrx,Ntx,nRelay+1);
% H = [-0.610710445329405 + 0.613567962680578i 0.853773080642614 + 0.401821553612705i -0.479681376045679 + 1.00531898850104i 1.14310287532886 - 0.117885524292688i;0.178246645851937 - 1.72208790285052i 0.305956954838553 - 0.312359382294011i -0.871974205950409 - 0.705563669353526i -0.998407932329956 - 1.08381078397046i;0.492684720867831 + 0.732042232546813i -0.851334162541907 + 0.110288365304141i 0.633536086860089 - 0.196817668059379i 1.16880843416264 + 0.680162319651726i;0.294620330350429 - 0.571855872996243i 0.930679850111546 + 0.701270409058833i -0.112364863358604 + 0.317726719211475i 0.592422553619285 + 0.372076545727887i];
% H = [0.936630064473899 - 0.912072269577213i 0.178557175611568 - 0.550044537958612i 0.579318237145130 - 0.354160686046319i 0.452235079738647 - 0.0531907362860184i;0.00720734199208105 - 0.642714631387622i -0.969596383080590 + 0.406597866363869i 0.245912824849761 + 0.150507859414785i 0.594257373589646 + 0.554690736467412i;1.06158501391750 + 0.550813813374430i 0.971424269091315 + 0.566198651710452i -1.26250845424579 + 0.703573560156512i -0.305906129183943 - 0.0172433790643602i;-0.629823633962719 - 0.756186183884771i 0.490304433274532 + 0.820205780753551i 0.449865609498062 - 0.614778562371777i -0.414106579513048 + 0.0337627349747282i];
% H = [-0.0239914655081684 + 0.194649530252264i -1.43215804518476 + 0.0450745302172133i -0.841779247976635 + 0.658401091971142i -0.0512188427914512 - 0.651949954227535i -0.524152031002036 + 1.08916097497483i -0.814553318813407 - 0.736110682954024i 0.385891183884599 - 0.278256516602899i -0.470821553371007 + 0.923559033269876i 0.815000490157656 - 0.400223727294702i 1.11322267418524 - 0.416223394195067i;0.431448292512599 - 0.0332839921141366i 0.628419136969722 + 1.63869431123794i -0.0408685019070174 + 1.25531429398936i 1.32520302595386 - 0.566109554180205i 0.329772232052704 - 0.762130190259398i -0.266952873276953 + 0.247682836555073i -0.290386510657718 + 1.26447260350458i -0.233566261805794 - 0.108199644505171i 0.810598460606792 - 0.911569679151346i -0.970593161169570 - 0.453850125771472i;-0.399933045003963 - 1.34707924449659i -0.249513908480544 + 0.663034440889015i -1.34218578838604 + 0.612167541987309i -0.132988424020837 - 0.00388902888764533i 0.823594764847632 - 0.0205097150546492i -0.0674552072572293 - 0.00843379546823615i -0.0118596929851988 + 0.0461742859942588i 1.05393565941984 + 0.108798293070356i -0.0959336977750042 - 0.00572643977861973i -1.33301777390874 + 0.938679593540409i;-0.388486163731857 - 0.127452063254815i 1.05555980585154 + 0.703262108424320i -0.297566966890022 + 0.375354612737671i 0.266137774133298 + 0.533931295827120i -0.661949622055118 - 0.568212311069179i -1.23353416243525 + 0.249822094335605i 0.0556623525079984 - 0.599512555073500i -1.19053648951469 - 1.00742765247239i -0.814698576456900 + 0.276186664830798i 0.201793266891921 + 1.03013098927613i;0.541236061689640 + 0.141610736510078i -0.368691917654018 - 0.313195797869534i -0.0202637981543274 - 0.857796438809843i 0.110160263468872 + 0.799550948110619i 0.184473087070863 + 0.294637593461297i 0.335708909620362 - 0.241074151880072i 0.772755799750676 - 0.205171550033189i 0.879929455137810 - 0.939372150401960i 0.516471419802982 + 0.738550157150376i -0.627232100268169 - 0.614795986406037i;0.00366565586244431 - 0.288336387097609i 0.708100832395293 + 0.288582819972268i 0.146643571442741 - 0.692633844888844i -1.01164930594226 + 0.0702313780176091i 1.06692388669230 - 0.352098686590135i 0.186086787577025 + 0.364667872908212i 0.757196969973072 + 0.945297255202440i -0.588078534830282 - 0.0562316112377496i -0.432338413945014 + 0.456221455115524i -0.587959154666459 + 0.631312630793623i;-0.536459504104915 - 0.923943737067690i -1.00953167426777 - 0.605603932088439i 0.0833618286421073 + 0.607733982414243i -0.445087831530378 + 0.300478691529943i 0.652759623025972 - 0.274208967478296i 0.661040484520370 + 0.283084987021165i 0.486548762655323 + 0.640381518459735i -0.252603033429638 + 0.830235253593659i 0.856930275636910 - 1.17292859408674i -0.186633380019648 + 0.822808147366479i;-0.403696711364464 + 1.03809987344836i -0.0830251300939353 - 0.144767335822325i -0.0216541369948180 + 0.124435521574517i 0.808746242810276 + 1.12549343700381i 0.466595045182514 + 0.333151739275127i -0.833716501339321 - 0.683531692303479i -0.714868192194049 - 0.0941424108389944i -0.0149727938195659 + 0.628047138606443i -0.995164874918881 + 0.251518193126372i 0.544963054696021 - 0.994819777084080i;-0.256218688649324 - 0.0912815682805835i 1.05183915342211 - 0.284594288658755i 0.0322741440191761 - 0.00736942090904399i 0.660457877641500 - 0.122211412985660i -0.792158312332550 + 0.221403913780932i -0.804984232132473 + 0.269787074528931i 0.325331433864884 - 1.02616154967149i 0.572531562491188 - 1.51887455054477i -0.0456256228699999 - 0.239273762164524i 1.37444182629057 - 0.708388272645869i;-0.348914014411850 + 0.216646679389860i -1.61049074425415 - 0.587046463359394i -0.814096759088603 - 0.193248501973912i -0.301302731464639 - 0.305815272223345i 1.83754805255708 - 0.320089730367297i -1.00342820628160 + 0.980259777535161i -0.430171749566233 - 0.487736130661123i -0.391454364811382 + 0.721645856521439i -1.46604856542110 - 0.421668841751681i 0.806764666640648 - 1.02141450793055i];
% H = [0.979836965402952 - 0.288473477604005i 0.299955305809267 - 1.12302971003183i 0.0777241214195054 - 0.622915673294521i 0.00848047892110321 - 0.532793238363008i;0.161274883964821 + 0.850097605374137i 0.507624224520579 + 0.0904340490278794i 0.762622732529626 + 0.636102162515930i 0.627070978815529 - 0.811139159330881i;-1.14898287263612 + 0.239039741634630i 0.0586675058679010 - 0.356613731661867i 0.444157599650455 - 0.666471382382458i 0.439853378882175 - 0.698485212863829i;0.221279201503905 + 0.767600933143939i 1.00413849732678 + 0.975942601516493i 0.0443986693492905 - 1.98448068398915i 0.174048853649218 + 0.0112119189769208i];
% H = [-0.485709649428944 - 0.443035449123101i -0.683877690479277 + 0.0625847243937286i 1.10497528104600 + 0.232615776714010i 0.548863136329041 - 0.214567328807759i;0.698599476873045 + 0.215158029994060i 0.889432317061025 + 1.17789953177362i -0.778191695472233 + 0.531085594495243i 0.547083342824163 - 0.323070548773543i;-0.0353267881771935 + 0.0136458588359820i 0.608364876377415 + 0.873133059209399i 0.627255032094848 - 0.272981793644085i -0.304160136306646 + 0.907819102383105i;1.46855905346172 - 0.433273646197841i -0.170005266999019 + 0.964420659175331i -0.856717617713001 + 0.446638209795896i -0.263524882990195 + 0.410472387366076i];
% H=[[-0.523178297241429 - 0.960626063813382i -1.26819896740364 - 0.133496370721499i -0.604958342028620 - 0.258561461857661i 0.482495513776905 - 0.356008959962932i 0.622724841307276 - 0.444781698583983i 0.302433545752687 - 1.17594469240269i 0.366447191024225 + 0.497066896657707i -0.213578037616910 - 0.539572793104541i -0.663153253933177 - 0.483047826475256i -0.423235216610257 - 0.164464469142183i;-0.529499356017282 - 0.414860895080231i 1.41210506818524 - 1.45603446748661i 1.77090167518602 - 0.977646500748111i 0.912330299215993 + 0.633267803485333i 0.198971171778270 - 0.401322690356326i -0.574423781535317 + 0.263373335527370i 0.814546875809817 - 0.417895986482605i 0.295037496877922 - 1.64964683751242i 0.799291617506264 - 0.738658129007354i -0.500114450399464 - 0.00310800112087012i;-0.952922775393404 - 0.168819729463332i 0.0238152422646001 + 0.592869095810843i -0.670038614464918 - 1.22745287890084i -0.476449796262484 + 0.798929140170980i -0.478750635593925 + 0.446440427003831i -0.225806619646297 + 0.00783542895030736i 0.267144239484859 + 1.19753575925681i 0.389427529552528 + 0.00340655136937681i 0.603471073403849 - 0.530511402329159i -0.241231941157954 - 0.487650892873229i;-0.104332715406209 + 0.309424156220941i 0.475443334697664 - 1.11412948808272i 1.17085281902098 - 0.762137795406681i -0.107034194356776 - 1.16964048366533i 0.172279467992998 - 0.284337300585645i -0.141385168801651 + 0.690607656355700i -0.369008322240325 - 0.516131999041804i 0.967407338422865 + 0.0542297007634695i -0.739974096030759 - 0.517803101056452i 0.377731229694409 + 0.827117082254196i;0.445879316390665 + 1.17494972437187i -0.584785191870063 - 0.724894819602419i -1.57753740199674 - 0.610629703594081i -0.0324703769754375 - 0.667433119704392i -0.0971748993090017 + 0.322322019702816i 0.542857003888012 + 0.00685421915218053i -0.600575167330994 - 1.03224376418064i -0.249725036097672 - 0.472917407827417i -0.807797503412611 - 0.965088385307907i -0.140766995950466 + 0.283996218354065i;0.0176616122858007 + 1.09492246264706i -0.0921958479692223 + 0.959252819002027i 1.25130633874525 - 0.164176825658639i 1.14430042989449 + 0.309564190259775i 0.719880601097874 - 1.81113873766141i 0.591215875442923 - 0.461058235294519i 1.13500792471961 - 0.407416539550690i -0.0378429544473133 - 1.36000578993466i -0.361613580043176 - 0.467935095172244i 0.0814706666561990 - 1.50943662004110i;-1.05285120462059 - 0.320942707838713i 0.143143097713817 + 0.928309927632521i 0.413929083076663 - 0.972436454173446i 0.0790557891325688 - 0.527080810913327i 0.117050595633534 - 0.761684406088064i -0.575264097668221 - 0.132615907264356i -0.0515084301222391 + 0.179667570804390i -1.12867053533527 - 1.27535761847903i 1.57716906283405 - 0.393724470059327i -0.193324286292313 + 0.578346242640949i;0.869716083764414 + 0.0267333985108625i -0.415134342008354 + 0.0473602543655001i 0.629325247714338 - 1.13653653129346i -0.144157092224316 + 0.313795406221448i 0.130482409675458 + 0.200206634603401i -0.0555226395124413 - 1.24874607248829i 0.440696028765043 - 0.628923331249478i 0.387654412345996 + 0.794050903841801i -1.16174312135247 - 0.957250317098100i 0.557744856135060 - 1.03031194217093i;-0.105218544404042 - 0.528360140772702i -1.05936937992179 + 1.05966626768449i -0.0840811164294009 + 0.418468972220684i -0.222870766418406 - 0.483426463264424i 0.335726864377729 + 0.529988638607328i -0.672205937274449 - 0.116899442855666i 0.0279051161285221 - 0.0505804727228947i -0.559679647271218 + 0.121303530949032i 0.191350605906208 - 0.205089689505836i -0.507139374650879 + 0.235609109918127i;-0.630718999903788 - 0.377745304222722i -0.305055468540781 + 0.297807869430721i 0.447405992949762 - 0.602640756135306i -0.567808567088105 - 0.852878516649231i 0.568004114471823 + 1.01638130478738i -0.172812701668263 - 0.916816290398520i -0.822782965261340 + 0.0147694051144036i -0.221776997239493 + 1.07373285444925i -0.602660848775378 - 0.541700017831029i -0.176672522525631 - 0.811804899080195i]];
% H = [[-0.384992822782392 - 1.57888977615430i 0.0568850266658465 - 1.83749534341371i 0.731793127779554 + 1.20828073021787i -0.258838611223545 - 0.628273415512718i;-0.170035611727804 + 1.36171659500041i -0.534830808273811 - 0.964053061495748i -0.398198624829146 - 0.261021122353217i -0.129048794082311 + 0.0369134313914847i;-0.296618106307074 - 0.0766088553573626i 1.73253969953794 - 0.688116405661226i -0.0104520253350722 + 1.45248433475641i -0.531395679783236 + 0.364910132239963i;-0.997707406320021 - 0.755736654301369i 0.515656756233553 + 0.959570134256285i -0.526172144480472 + 0.281256751902015i 0.531018591169660 - 0.0859518796577853i]];
% H = [-0.542882751792271 - 0.304809049165547i 0.154246747708710 + 0.827119100964630i 1.05571623623117 + 0.832595381389597i -0.536408262069608 - 0.613965714946982i 0.472296725884552 - 0.402718824281848i;0.509605376268358 + 0.132229247651890i 1.98639275696402 - 0.238476586014611i 0.00624254386138586 - 0.689557680850944i 1.70363801619731 - 0.889812155543451i 0.607702137687825 - 0.978051942146919i;0.138576143627612 + 0.0739989779481696i -1.38197744138537 + 0.217572588061087i -0.400056223719364 - 0.407547147718321i 1.18303673507476 - 0.706510366073804i -1.02279541860476 - 0.170075591953002i;0.831213129570521 + 1.05325761892170i -0.605367070683721 + 1.63522999781038i -0.208527608307871 - 0.164665615788134i -0.754215633778035 - 0.360299966185866i -1.78589081082087 - 1.31277703600896i;-0.980612278226565 + 0.345058341829933i -0.114297512018638 - 0.228946249660042i 0.184344504522392 - 0.255942245310648i -0.729590097016352 + 0.659454736206250i -0.283689703854611 - 0.507409884937575i];
% H = [0.412246458790627 + 0.353406262591349i -0.286724963452034 - 0.157984753755752i -0.292012036774348 - 0.187617052190734i -0.559724125991516 - 0.213995204565634i -0.147739782742775 - 0.206832282406856i -0.451568769435229 - 0.240942414304644i -0.258072094226792 - 0.0114111603100531i 0.326872394536588 + 1.85259088938607i 0.941207114914868 - 1.14102529191406i 0.362832159161409 + 1.62674388566850i;0.269479012769564 + 0.808517194392612i 0.362590304576200 + 0.992538059000057i 0.494682628187656 - 0.172176938487343i -0.175107600860572 - 0.133272684643150i -0.845058637989381 + 0.123245955534721i -0.811462324884678 + 0.822912890962759i -0.771104841405928 - 0.425731607399610i -1.12877296054910 + 0.958793913816756i 0.115548412834946 - 0.205857124256646i -1.37652913127261 + 0.319088324793716i;0.342586698144773 + 0.0738278890126121i -0.833747930588939 + 0.993904088228871i 0.773323487482860 - 0.0687535683434592i -1.00932842927428 + 1.27296773010432i 0.112142840120187 + 0.0525101220827832i -0.457700700771783 - 0.744954061366906i -1.16946550056771 - 0.932250027541390i -0.719719424645080 + 0.643652974050636i -0.115929859151005 + 0.289806337984472i -0.218573889069872 - 1.04970673546006i;-0.00335112183440466 - 0.516393715962562i 0.327135910470642 - 0.386071991904642i -0.409103886355127 + 0.260989321536051i 0.0207102855897127 + 0.258053330501774i 0.593402651650213 + 0.247603943618532i -0.314882377405081 + 1.42330457767519i 0.113416383911172 + 0.887655259498316i -0.402248763600146 - 0.339214191627547i -0.663960873556412 + 0.418032186277265i 1.32474728266637 - 0.0986663240651964i;-1.15067085339206 - 0.0741920191754055i 0.339902398522088 + 0.434041662060957i -0.574315478390875 + 1.17622085008619i 0.532049075927418 - 0.609561322999834i -0.138001888221301 + 0.550777494223502i 0.167592015283854 - 0.653725818098379i 0.923882340010243 + 0.268893330780900i 0.0946356242447649 + 0.391087730392493i 1.26283549057319 + 0.640358211655134i 0.874959797651889 + 0.839917496390940i;-0.0490412170639579 - 0.778967204137249i 0.335263758272431 - 0.198389768979120i -0.179124114892756 + 0.737538545696026i 0.0527327977628948 - 0.102693241444766i 0.535486285636594 - 0.170954588201088i -0.276567707681411 + 1.11906678661444i -0.322205158038962 + 0.109591001600964i -0.566218375610338 + 0.0398465691821864i 0.489646022269407 + 0.106612869342907i 0.104313419941040 - 0.394793619711118i;0.418874291439889 - 0.515043728068719i 0.790583122401083 - 1.14020480219202i 0.375930911699241 + 0.105117336088667i 1.22629975970944 + 0.498281494430148i 0.568218344078594 - 0.236992341547193i 0.363350528991052 - 1.25642988245845i 0.489985315750614 - 0.419643927077922i -0.487509585809037 + 0.250369344784144i 0.224065963496744 - 1.46405068666783i -0.822073389829111 + 0.795166768736863i;0.210615773380314 + 0.403330963151998i -0.610075687670280 + 0.000530532634592698i -0.184628598192023 + 0.0786007709098886i 0.340390821853655 + 1.63732584344081i 0.446189016792203 - 0.843491749962551i 0.104220093552750 + 0.605847448342868i 1.30147563819753 + 1.11354961565162i -0.384791160240710 + 0.123110286254831i 0.531273373185248 + 0.329513390007900i -1.02734996818236 - 0.132876159279368i;0.646455825861981 - 1.83004226049428i 0.245414883756264 - 0.657934044566066i 0.455501426938201 + 0.134084483220587i 1.01387887642235 - 0.240419219948592i -0.239741072148086 + 0.161831870505540i -1.25775912413851 - 1.10494607997661i 0.724115395627145 + 0.162845856051067i 0.316519241773284 - 0.302296728518881i -0.536668798915112 - 0.0607102103936126i 0.633432278078489 - 0.720776382534473i;-0.316428848898081 - 0.146405764175596i -0.343828939611140 - 0.285524433288874i 0.447122571525987 + 0.0274471475297182i 0.558009285970865 - 0.0152604507753512i -0.178347212264207 + 0.495413872920281i -0.685789978141891 + 0.168818688721364i -0.165768164400746 + 0.206719940807860i -0.101442784023576 - 0.0982019554219953i 0.960375985648100 + 0.398157291281277i -0.333982729824137 + 0.419479572600230i];
% H = [0.412246458790627 + 0.353406262591349i -0.286724963452034 - 0.157984753755752i -0.292012036774348 - 0.187617052190734i -0.559724125991516 - 0.213995204565634i -0.147739782742775 - 0.206832282406856i -0.451568769435229 - 0.240942414304644i -0.258072094226792 - 0.0114111603100531i 0.326872394536588 + 1.85259088938607i 0.941207114914868 - 1.14102529191406i 0.362832159161409 + 1.62674388566850i;0.269479012769564 + 0.808517194392612i 0.362590304576200 + 0.992538059000057i 0.494682628187656 - 0.172176938487343i -0.175107600860572 - 0.133272684643150i -0.845058637989381 + 0.123245955534721i -0.811462324884678 + 0.822912890962759i -0.771104841405928 - 0.425731607399610i -1.12877296054910 + 0.958793913816756i 0.115548412834946 - 0.205857124256646i -1.37652913127261 + 0.319088324793716i;0.342586698144773 + 0.0738278890126121i -0.833747930588939 + 0.993904088228871i 0.773323487482860 - 0.0687535683434592i -1.00932842927428 + 1.27296773010432i 0.112142840120187 + 0.0525101220827832i -0.457700700771783 - 0.744954061366906i -1.16946550056771 - 0.932250027541390i -0.719719424645080 + 0.643652974050636i -0.115929859151005 + 0.289806337984472i -0.218573889069872 - 1.04970673546006i;-0.00335112183440466 - 0.516393715962562i 0.327135910470642 - 0.386071991904642i -0.409103886355127 + 0.260989321536051i 0.0207102855897127 + 0.258053330501774i 0.593402651650213 + 0.247603943618532i -0.314882377405081 + 1.42330457767519i 0.113416383911172 + 0.887655259498316i -0.402248763600146 - 0.339214191627547i -0.663960873556412 + 0.418032186277265i 1.32474728266637 - 0.0986663240651964i;-1.15067085339206 - 0.0741920191754055i 0.339902398522088 + 0.434041662060957i -0.574315478390875 + 1.17622085008619i 0.532049075927418 - 0.609561322999834i -0.138001888221301 + 0.550777494223502i 0.167592015283854 - 0.653725818098379i 0.923882340010243 + 0.268893330780900i 0.0946356242447649 + 0.391087730392493i 1.26283549057319 + 0.640358211655134i 0.874959797651889 + 0.839917496390940i;-0.0490412170639579 - 0.778967204137249i 0.335263758272431 - 0.198389768979120i -0.179124114892756 + 0.737538545696026i 0.0527327977628948 - 0.102693241444766i 0.535486285636594 - 0.170954588201088i -0.276567707681411 + 1.11906678661444i -0.322205158038962 + 0.109591001600964i -0.566218375610338 + 0.0398465691821864i 0.489646022269407 + 0.106612869342907i 0.104313419941040 - 0.394793619711118i;0.418874291439889 - 0.515043728068719i 0.790583122401083 - 1.14020480219202i 0.375930911699241 + 0.105117336088667i 1.22629975970944 + 0.498281494430148i 0.568218344078594 - 0.236992341547193i 0.363350528991052 - 1.25642988245845i 0.489985315750614 - 0.419643927077922i -0.487509585809037 + 0.250369344784144i 0.224065963496744 - 1.46405068666783i -0.822073389829111 + 0.795166768736863i;0.210615773380314 + 0.403330963151998i -0.610075687670280 + 0.000530532634592698i -0.184628598192023 + 0.0786007709098886i 0.340390821853655 + 1.63732584344081i 0.446189016792203 - 0.843491749962551i 0.104220093552750 + 0.605847448342868i 1.30147563819753 + 1.11354961565162i -0.384791160240710 + 0.123110286254831i 0.531273373185248 + 0.329513390007900i -1.02734996818236 - 0.132876159279368i;0.646455825861981 - 1.83004226049428i 0.245414883756264 - 0.657934044566066i 0.455501426938201 + 0.134084483220587i 1.01387887642235 - 0.240419219948592i -0.239741072148086 + 0.161831870505540i -1.25775912413851 - 1.10494607997661i 0.724115395627145 + 0.162845856051067i 0.316519241773284 - 0.302296728518881i -0.536668798915112 - 0.0607102103936126i 0.633432278078489 - 0.720776382534473i;-0.316428848898081 - 0.146405764175596i -0.343828939611140 - 0.285524433288874i 0.447122571525987 + 0.0274471475297182i 0.558009285970865 - 0.0152604507753512i -0.178347212264207 + 0.495413872920281i -0.685789978141891 + 0.168818688721364i -0.165768164400746 + 0.206719940807860i -0.101442784023576 - 0.0982019554219953i 0.960375985648100 + 0.398157291281277i -0.333982729824137 + 0.419479572600230i];
% H = [0.422908024187963 - 1.09670721038489i 1.36852474771961 + 0.299257967508388i 1.20779955067076 + 0.598916182787320i -0.205515684504568 - 1.22031966193830i -0.416905376221829 + 0.799693967125371i;0.497248379187511 - 1.06600915296125i 0.0327137893222458 - 0.746668655614259i 0.682654346999734 - 0.101451306490223i -0.569649301050501 + 0.851325624855875i -0.850752958101867 + 0.865419080333930i;0.262144618517377 + 0.318059851692015i -0.248576304331525 - 0.668927549074541i 0.394989427122515 - 0.602125051887423i 1.22627156922223 - 0.728167373821871i -0.637648773857439 - 0.277990667677902i;-0.0560791073734730 - 1.44248755699904i 0.984177489533024 - 0.409334219126562i 0.983012602205002 + 0.356609383012099i 1.66065745273484 - 1.37520439362181i 0.0405186686695813 + 0.487644665125851i;0.0872489445817198 + 0.347843092400993i 0.818261134526316 - 0.450901880644853i -0.117235194579301 + 0.693043144287049i 0.719304536374384 + 0.0529269275335984i -0.909558484346915 - 0.0877003047814881i];
% H = [1.33269263810989 + 0.408555345215781i 0.269433272773026 + 0.0569282151662241i 0.653812991303650 + 0.177125556417211i 0.266779110729647 - 0.244063842439694i -0.410396799503188 - 0.892270689948754i 0.424772007796446 + 0.0869737994199469i 0.940688028881620 - 0.463505828304910i -0.130235835863900 + 0.666454276581022i -1.00410375984609 + 0.249773690607900i 0.315828215786470 - 0.527907336115938i;-0.690467824787123 - 0.424833009664706i 0.337453206476451 - 0.103742198130051i 1.59189567892011 + 0.949291025376143i 0.272324262574478 + 0.644675889900994i -0.646631554534903 + 0.0979381654493152i -0.982409660296009 + 0.570730460025372i -0.427323801464746 - 0.117311970247474i 0.0395218213057966 + 0.691501532426420i 0.538851749002549 + 0.180694202385139i -0.793872937112927 - 0.180022247736220i;0.113544417496756 + 1.59052071193629i 0.777644262192577 + 0.893513410568357i -1.25375574331602 - 0.447987318738497i -0.247953564370057 + 1.06401694410135i -0.103200583779596 - 0.471141329044659i -0.373689577878358 + 0.0862313775142984i -0.960136183165212 - 1.09727933085160i 0.0168302194418440 - 0.350753620541534i -1.56360418190802 + 0.572839042538491i 0.368074889814394 + 0.183969330420280i;0.906962360369035 - 0.817369487882631i 0.785766969543866 + 0.745983373693680i -1.30245462423229 - 0.237579277013020i 0.0473430552531045 - 0.355725713195805i 0.345910323523835 + 0.471204704033887i -0.158331344150551 - 0.384399232937529i 1.12373593039877 + 0.745903210640407i 0.243478057497025 - 0.484283008997666i 0.634207898644753 - 0.460541908077982i 0.125928621262555 - 0.0924098379504425i;0.631274584485841 - 0.763394476606790i -1.14547972125294 - 1.19868646648457i -0.0226574196612785 - 0.461340558908347i 0.0801055710815921 + 0.137816812108863i 0.110289870097596 + 0.138517056138363i -0.135518627502674 + 0.954725286962484i 0.387189049863655 - 0.0113123402268469i 0.955080738112588 + 0.804826804717663i -0.644409335172225 - 0.753340016936977i -0.0562709932266179 + 0.449772817847265i;-0.437747288806351 - 0.433104182718948i 0.219788112446059 + 0.567975994500800i 0.275109847365351 - 0.743763798294337i -0.473950601160555 + 0.232610953270901i -0.683628093213313 + 0.968319148839914i -0.0343374634508152 - 0.201137722123267i 0.570843253454773 + 0.395973586881842i -0.332508394827185 - 0.420491321183684i -0.472621876822735 + 0.475591972650846i -0.661727659725195 + 0.569178162954534i;-0.455153650024527 - 1.04084842891983i -0.0713803719408492 + 0.281072166949436i -0.191860935968712 - 0.265983052420725i -0.748256918754431 - 0.775421871291425i 0.619851069212114 + 0.00958930434559191i -0.0446420558758532 + 0.991874943015270i 0.157933025277541 + 0.199780111539474i -1.12668990935256 - 0.680069371802344i 0.321578766210009 + 0.203328485929610i 0.501838758427329 + 0.539874435543446i;-0.125752460587519 - 0.0653143340311745i 0.0955390902079189 + 0.978577679971439i 0.758923913297168 - 0.487456338587646i 0.259291488589421 + 0.332135080316020i -0.650518564581885 + 0.332312604249619i -0.260093145813208 + 0.376987463777276i 0.826948517323638 - 0.0239097485794937i 0.660321736069386 - 0.217164491056055i -0.612482344616219 - 0.525346240502409i -1.29771979275782 + 1.02228943461204i;-0.100759816561332 + 0.488818505915324i 0.896336801741806 + 0.272038449017886i 0.698393494267256 + 0.897114795366478i 0.637084269781616 - 0.304870648358032i 1.44250556432177 + 1.11019865597052i -0.561089313202071 + 0.372799241981488i -0.774272010354240 + 0.0783936097343152i -0.245707434311301 + 0.665585507771284i -0.174030535402054 - 0.0549338366098522i 0.606359134650009 - 0.428501702086666i;0.601128553900310 + 0.0992321030415227i -0.693114748510136 + 0.230923551438579i -1.32021388771105 - 1.28216821626371i 0.187915719156324 - 0.323472079985223i -0.472601166628645 - 1.09774560673423i 0.565710248423319 - 0.190999211503153i 0.204194218017173 + 0.569469088878071i 0.655868667509088 - 0.0325621609686216i 0.688523326377659 + 0.887244002728576i 0.770500176063290 - 0.245443863656992i];
% H = [0.356030251232995 + 0.176366131797958i -0.553549499225168 - 0.148502413357783i -0.839623160283537 + 0.171871829436033i -0.560213258441509 - 1.85264773076677i -0.530288557100176 - 0.604488771760839i;-0.304041505134952 - 0.343228102047720i -0.000216532334654556 + 0.388002742535461i -0.984565925088192 - 0.125998726256707i 0.116344534800558 + 0.433767771120836i -0.0416536159225009 - 0.210809740504384i;-0.472917824798041 - 0.173193804640230i -0.135303152805121 + 0.389506137629916i 0.327970103937801 + 0.502198357735333i -1.04869954182585 - 0.647165502184525i 0.794282158981899 + 1.02437349650657i;-0.306241751345337 + 0.418997258555094i 0.223257516771372 + 0.374451752159181i 0.141729800956831 + 0.151602885324946i 0.662185281616636 - 1.42587322780420i 0.528954955230633 + 0.739464441472211i;0.297047058072866 - 1.01692168854246i 0.735563739775095 - 1.34648392265127i 0.541284988646917 + 1.07997969638610i 0.786652564215026 - 1.33883820148693i 1.28268459920033 - 0.366375571087866i];
% H = [-0.918108099175772 + 0.258956811676495i 0.561239759310838 - 0.659499219841714i -0.440754596815474 - 0.400867492587136i -0.592130391128383 - 0.755059898679473i -0.619899930658449 + 0.592720692461061i 0.219712790709775 + 0.865950539546335i 0.0405401550110543 - 0.444076030660745i 0.788803603095980 - 0.0514681858160361i 0.106645886476731 - 1.43904806028899i -0.0875302061800762 + 0.510656951950410i;1.42587240569684 + 0.583432329598388i -0.339421631359256 - 0.167582319779919i -0.351188135952944 + 0.331716041247333i -0.699329781933597 + 0.0846516784829853i 0.495809105739074 + 0.696038292566358i -0.666821535013628 + 0.236506784425628i -0.124813877961377 - 0.956123278882529i -0.685567592905650 + 0.284015164206481i 0.508502859118905 + 0.184243646471843i -0.0329245306325335 - 1.18624339953247i;0.322319538498723 + 0.894124792424777i -0.226030087855838 - 0.235599838079481i 0.732132598943456 + 0.472767877461633i 0.833508681910382 - 0.445916142876023i -0.434854830088828 + 0.397595753972340i 0.395702310035311 - 0.648044862268243i 0.179151364996762 - 0.405049817472229i 0.531789109811438 + 0.848504656650422i 0.350598607359591 - 0.733672720127338i -0.174130992199138 + 0.352461665339853i;0.683657540615195 - 1.52649590003463i 0.0450087189856489 - 0.770190835130551i -0.535688420884718 - 0.852490674175692i 1.05968801870817 - 0.802816383371322i -0.125236000671164 + 0.792412472415800i -0.505216917852659 - 0.148910877434014i 0.00165046203018977 - 0.134726904938056i -0.640505914017623 - 0.111912503946224i -1.70271712449650 - 0.998169429581627i 1.12282865504888 - 0.406715903246433i;-0.518613287945126 + 0.507170480849967i 0.802926232978638 + 0.799435200483512i 1.30979498833584 + 1.59673411241685i 0.457116297095648 - 0.169169157310416i 0.974388483543074 + 1.59687598969844i 0.375815635396421 + 0.478514589169743i 1.25818473389950 + 0.248633314349890i 0.486074044508707 + 0.533971841328369i 0.896094096665626 + 0.432085233863151i -0.223378479427238 - 1.73402593643675i;-0.146532712668339 + 1.59912306052106i 1.13000474472464 - 0.00935579343327753i -0.115043436377111 + 0.0249584292235459i -0.164939957504828 + 1.25869867813806i 0.515960783993968 + 0.534233498329647i -0.630133002122477 + 1.14911563816874i -0.288217012704462 + 0.0719291163317987i -0.257005896742044 - 0.344793499550368i -0.139392903669537 - 0.239207070870824i -0.385058931620986 + 0.336515765742659i;-0.601222285195663 + 0.432087387627704i 0.624216499812947 - 0.0124387535359494i 0.574795715772037 + 0.309141885979891i -1.05312164375059 + 0.810359419325752i -0.904615968864665 + 0.367295843712598i -0.832298509387351 - 0.252517780193735i -0.0790497266857378 + 1.05691903223872i 0.00785167732584626 + 0.208941489744235i -0.229982159753931 + 0.684904934930192i -2.03346058178249 + 0.942258934335588i;-0.370175239299181 - 0.146853988985475i 0.112486244717315 - 1.30085007901809i -0.454623813992251 - 0.302213701575057i -1.24599235161847 - 0.311325058854938i -0.249734858688756 + 0.0337575667269029i -1.11319995103832 + 0.610294346620625i 0.0331689441125845 + 0.0954804953162591i 1.13240299661194 + 0.259063649123722i -0.231707967795963 + 1.28539034611004i -1.14480227763020 - 0.521700108753708i;0.146880912383789 + 0.330231966917253i 0.440677826877988 + 0.405577337948797i -0.0822172142478361 + 0.144641160025988i 1.37757546639081 - 0.840667057540610i -0.711237358362912 + 0.614037124394649i 0.00466092250067576 - 0.608015349125285i -0.380462943601383 + 0.441795883510352i 0.853617456752362 - 0.0479623944230664i -0.537805245005439 - 0.637834610832023i 0.380394667749686 - 1.13099497503607i;0.289407771372229 + 1.27278090322634i 0.644816015816626 + 0.376246429918693i -0.175933269072172 + 0.456867641440277i -0.476379468999254 + 0.640348379277749i -0.0777972676728269 - 1.21047682386273i 0.546177172194797 - 0.142459462988153i 0.385853377748800 - 0.233480048526528i -1.38078406688673 + 0.443761688165279i 0.116752762973992 - 0.423776796894050i 0.155873005651947 + 0.451734600576380i];
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
% return
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
        if  strfind(mode(1:2),'WF')
            while 1
                if strfind(mode,'Separate');
                    if iRelay == 1
                        Q(:,:,iRelay) = eye(K);
                    else
                        %                         F = V(:,1:K,iRelay-1)*sqrt( diag(P(1:K,iRelay-1)) );
                        Q(:,:,iRelay) = G(:,:,iRelay-1)*H(:,:,iRelay-1)*F(:,:,iRelay-1)*Q(:,:,iRelay-1)*F(:,:,iRelay-1)'*H(:,:,iRelay-1)'*G(:,:,iRelay-1)'+G(:,:,iRelay-1)*G(:,:,iRelay-1)';
                    end
                elseif strfind(mode,'Independent');
                     Q(:,:,iRelay) = eye(K);
                end
                P(1:K,iRelay) = zeros(K,1);
                isBadChannel = 0;
                for iChannel=1:K
                    P(iChannel,iRelay) = ( P0+sum(1./eigValue(1:K,iRelay)) )/( sqrt(Q(iChannel,iChannel,iRelay))*singValue(iChannel,iRelay)*sum( sqrt(diag(Q(:,:,iRelay)))./singValue(1:K,iRelay)) )-...
                        1/( Q(iChannel,iChannel,iRelay)*eigValue(iChannel,iRelay) );
                    P(iChannel,iRelay) = real(P(iChannel,iRelay));
                    
                    if P(iChannel,iRelay)<0
                        P;
                        sum(P);
%                         sum(1./( P.*eigValue(1:size(P,1))+1) );
                        1;
                    end
                    
                    P(iChannel,iRelay) = max( [P(iChannel,iRelay) 0] );
                    
                    if ~P(iChannel,iRelay)
                        P;
                        K = K-1;
                        isBadChannel = 1;
                        break
                    end
                end
%                 sum(1./( P(1:K).*eigValue(1:K)+1) )
                if ~isBadChannel %iChannel == K
                    F(:,:,iRelay) = V(:,1:K,iRelay)*sqrt( diag(P(1:K,iRelay)) );
                    G(:,:,iRelay) = inv( diag(eigValue(1:K,iRelay))*diag(P(:,iRelay))+inv(Q(:,:,iRelay)) )*...
                        sqrt(diag(P(:,iRelay)))*diag(singValue(1:K,iRelay))*U(:,1:K,iRelay)';
%                     sum(diag(F(:,:,iRelay)*Q(:,:,iRelay)*F(:,:,iRelay)'))
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
        elseif strfind(mode,'uniform')
            P(:,iRelay) = P0/K*ones(K,1);
            iRelay = iRelay+1;
        elseif strfind(mode,'aligning')
            % no Testing
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
        % no Testing
        
    end
    
    
end
P(K+1:end,:) = [];
return

for iRelay = 1:nRelay
    
end

