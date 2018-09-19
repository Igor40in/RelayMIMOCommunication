clc
clear all
close all
% h = openfig('tmp2.fig');
h = openfig('F:\MIMO_Communication\result\CollaborativeOptimization\for paper\compare_Relay5x10_UniformMorePointsWithoutMIMO.fig');
axes = get(h,'Children');

% iAxes = 1;
% ax = axes(iAxes)
% data = get( ax,'Children');
% Xdata = cell2mat(get(data,'XData') );
% Ydata = cell2mat(get(data,'YData') );
% Xdata(:,17:end)=[];
% Ydata(:,17:end)=[];


figure;

for iAxes=1:4
    ax = axes(iAxes);
    data = get( ax,'Children');
    Xdata = cell2mat(get(data,'XData') );
    Ydata = cell2mat(get(data,'YData') );
    Xdata(:,13:end)=[];
    Ydata(:,13:end)=[];
    nLine =size( Xdata,1);
    if iAxes==4
        h = subplot(2,4,[1 2]);
%         h = subplot(1,1,1);
    elseif(iAxes==2)
        h = subplot(2,4,[3 4]);
    elseif (iAxes==3)
        leg = ax;
    else
        h = subplot(2,4,[6 7]);
    end
    
    for iLine = 5:6%nLine
        hold on
        if iAxes==4
            Y(iLine,:) = log10(Ydata(iLine,:));
%             Y(iLine,:) = Ydata(iLine,:);
            p = polyfit(Xdata(iLine,1:1:end),Y(iLine,:),4);
            f = polyval(p,Xdata(iLine,1:1:end));
            YdataAppr(iLine,:) = 10.^(f);
%             YdataAppr(iLine,:) = f;
%             semilog
        elseif iAxes==3
            continue;
        else
            p = polyfit(Xdata(iLine,1:1:end),Ydata(iLine,1:1:end),4);
            YdataAppr(iLine,:) = polyval(p,Xdata(iLine,1:1:end));
        end
        plot(h,Xdata(iLine,1:end),Ydata(iLine,:),'LineWidth',2,'LineStyle','--')
        hold on
        plot(h,Xdata(iLine,1:end),YdataAppr(iLine,:),'LineWidth',2)
        grid on
        1;
    end
    
    ind = iLine/2;
    if iAxes==4
        title('Relay MIMO=5x5')
        legend(leg.String{end-(2*(ind-1)+0)},'',leg.String{end-(2*(ind-1)+1)},'')
        1;
    end
end