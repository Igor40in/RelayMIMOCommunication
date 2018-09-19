function [metric, Parameters] = NMDS_Optimization(parIn,parMin,parMax,parDelta,lamda,P0,optBudget)
nParam = size(parIn,1);
alpha = 1;
beta = 0.5;
gamma = 2;
delta = 0.5;
Parameters(:,1) = LogTrsf(parIn,parMin,parMax);
Vertex(:,1) = Parameters(:,1);
for iParam = 1:nParam
    Vertex_new = parIn;
    Vertex_new(iParam) = Vertex_new(iParam)+parDelta(iParam);
    Vertex(:,iParam+1) = LogTrsf(Vertex_new,parMin,parMax);
end
Interation = 1;
Steps(Interation) = 1;
Param_curr = Vertex(:,1);
Param_real = InvLogTrsf(Param_curr,parMin,parMax);
metricStep = metricEvaluation(Param_real, lamda,P0);
metric(1) = metricStep;
ObFun(1) = metric;

Iteration = 2;
metric(Iteration,1) = metric(Iteration-1,1);
value = 0;
for iVertx = 1:nParam
    Param_curr = Vertex(:,iVertx+1);
    Param_real = InvLogTrsf(Param_curr,parMin, parMax);
    metricStep = metricEvaluation(Param_real, lamda,P0);
    ObFun(iVertx+1,1) = metricStep;
    if ObFun(iVertx+1,1)<metric(Iteration,1)
        value = 1;
        metric(Iteration,1) = metricStep;
        Parameters(:,Iteration) = Param_curr;
    end
end
if value==0
    Parameters(:,Iteration) = Parameters(:,Iteration-1);
end
Steps(Iteration) = nParam+1;
% 'Start'
while Steps(Iteration)<optBudget
    Iteration = Iteration+1;
    if Iteration==8
        1;
    end

    Steps(Iteration) = Steps(Iteration-1);
    metric(Iteration,1) = metric(Iteration-1,1);
    
    [ObFunSort, indSort] = sort(ObFun);
    bestObFun =ObFunSort(1);
    worstObFun =ObFunSort(end);
    secWorstObFun =ObFunSort(end-1);
    indWorst = indSort(end);
    paramMean =zeros(size(parIn));
    
    paramMean = ( sum(Vertex,2)-Vertex(:,indWorst) )./nParam;
    % Reflaction
%     'Reflection',Iteration
    paramRef = paramMean+alpha*( paramMean-Vertex(:,indWorst) );
    paramReal = InvLogTrsf(paramRef,parMin,parMax);
    Steps(Iteration) = Steps(Iteration)+1;
    metricRef = metricEvaluation(paramReal,lamda,P0);
%     metricRef = metricStep;
%     ObFunRef = ObFunSort(1)+0.001
    if metricRef>=bestObFun & metricRef<secWorstObFun
%         strcat('Iteration=',num2str(Iteration),'_Reflection')
        Vertex(:,indWorst) = paramRef;
        ObFun(indWorst) = metricRef;
        if Iteration==14
            1;
        end
        if Iteration-1>size(Parameters,2)  
            'Bad'
        end
        Parameters(:,Iteration) = Parameters(:,Iteration-1);
        metric(Iteration) = metric(Iteration-1);
        continue;
    end
    % Expansion
    if metricRef<=bestObFun
        
        metric(Iteration) = metricRef;
        Parameters(:,Iteration) = paramRef;
        paramExp = paramMean+gamma*(paramRef-paramMean);
        paramReal = InvLogTrsf(paramExp,parMin,parMax);
        Steps(Iteration) = Steps(Iteration)+1;
        metricExp = metricEvaluation(paramReal,lamda,P0);
        if metricExp<metricRef
%             strcat('Iteration=',num2str(Iteration),'_Expansion +')
            Vertex(:,indWorst) = paramExp;
            ObFun(indWorst,1) = metricExp;
            metric(Iteration) = metricExp;
            Parameters(:,Iteration) = paramExp;
        else
%             strcat('Iteration=',num2str(Iteration),'_Expansion -')
            Vertex(:,indWorst) = paramRef;
            metric(Iteration) = metricRef;
            ObFun(indWorst,1) = metricRef;
        end
        continue;
    end
    % Outside contraction
    if metricRef>=secWorstObFun & metricRef<worstObFun
%         strcat('Iteration=',num2str(Iteration),'Outside contraction')
        Parameters(:,Iteration) = Parameters(:,Iteration-1);
        paramConOut = paramMean+beta*(paramRef-paramMean);
        paramReal = InvLogTrsf(paramConOut,parMin,parMax);
        Steps(Iteration) = Steps(Iteration)+1;
        metricConOut =  metricEvaluation(paramReal,lamda,P0);
        if metricConOut<=metricRef
            Vertex(:,indWorst) = paramConOut;
            ObFun(indWorst,1) = metricConOut;
            if metricConOut<bestObFun
                metric(Iteration) = metricConOut;
                Parameters(:,Iteration) = paramConOut;
            end
            continue;
        end
    end %go to Shrink
    % Inside contraction
    if metricRef>=worstObFun
%     strcat('Iteration=',num2str(Iteration),'Inside contraction')
        Parameters(:,Iteration) = Parameters(:,Iteration-1);
        paramConIns = paramMean+beta*(Vertex(:,indWorst)-paramMean);
        paramReal = InvLogTrsf(paramConIns,parMin,parMax);
        Steps(Iteration) = Steps(Iteration)+1;
        metricConIns = metricEvaluation(paramReal,lamda,P0);
        if metricConIns<worstObFun
            Vertex(:,indWorst) = paramConIns;
            ObFun(indWorst,1) = metricConIns;
            if metricConIns<bestObFun
                metric(Iteration) = metricConIns;
                Parameters(:,Iteration) = paramConIns;
            end
            continue;
        end
    end%go to Shrink
    %Shrink
%     strcat('Iteration=',num2str(Iteration),'Shrink')
    indBest = find(ObFun==min(ObFun));
    indBest(2:end) = [];
    for iVertx = 1:nParam+1
        if iVertx==indBest
            continue
        end
        Vertex(:,iVertx) = Vertex(:,indBest)+delta*(Vertex(:,iVertx)-Vertex(:,indBest));
        paramReal = InvLogTrsf(Vertex(:,iVertx),parMin,parMax);
        Steps(Iteration)=Steps(Iteration)+1;
        metricShr = metricEvaluation(paramReal,lamda,P0);
        if metricShr<ObFun(indBest)
            metric(Iteration) = metricShr;
            Parameters(:,Iteration) = Vertex(:,iVertx);
        end
    end
end

for i=1:Iteration
    Parameters(:,i) = InvLogTrsf(Parameters(:,i),parMin,parMax);
end
1;

function parOut = LogTrsf(parIn,parMin,parMax)
% parOut = log( (parIn-parMin)./(parMax-parIn) );
parOut = parIn;

function parOut = InvLogTrsf(parIn,parMin,parMax)
% parOut = parMin+(parMax-parMin)./(ones(size(parIn))+exp(-parIn));
parOut = parIn;

function metric = metricEvaluation(parIn,lamda,P0)
p = sum(parIn);
% if p>P0
% %     'More P0'
% %     parIn = parIn./p*P0;
% %     penaltyAdd = 0;
%     penaltyAdd = 2*(p-P0);
% else
%     penaltyAdd = 0;
% end
% metric = sum( 1./(parIn.*lamda+1) )+penaltyAdd;
% return
%%--------------------------------
indEdgeMore = find( parIn>P0 );
indEdgeLow = find( parIn<0 );
l = length(indEdgeMore)+length(indEdgeLow);

penaltyAdd = 0;
if l>0
    penaltyAdd = penaltyAdd+20*P0/1*sum( abs( repmat(P0,length(indEdgeMore),1)-parIn(indEdgeMore) ) );
    penaltyAdd = penaltyAdd+20*P0/1*sum( abs(parIn(indEdgeLow)) );
    penaltyAdd = penaltyAdd+2*(p-P0);
elseif p>P0
    penaltyAdd = 20*(p-P0);
else
    penaltyAdd = 0;
end
metric = sum( 1./(parIn.*lamda+1) )+penaltyAdd;

if (metric-penaltyAdd)<0
    metric = Inf;
end
% metric = sum( 1./(parIn.*lamda+1) )+penaltyAdd;

