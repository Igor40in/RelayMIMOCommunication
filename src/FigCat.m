clc
clear
close all
h1 = openfig('F:\MIMO_Communication\result\CollaborativeOptimization\for paper\compare_nRelay=4_5x10_Uniform.fig');
axes{1} = get(h1,'Children');
h2 = openfig('F:\MIMO_Communication\result\CollaborativeOptimization\for paper\compare_nRelay=4_5x10_shft.fig');
axes{2} = get(h2,'Children');
figure;
for iAxes = 1:4
    if iAxes==4
        h = subplot(2,4,[1 2]);
%         h = subplot(1,1,1);

    elseif(iAxes==2)
        h = subplot(2,4,[3 4]);
    elseif (iAxes==3)
        continue
    else
        h = subplot(2,4,[6 7]);
    end
    hold on
    grid on
    Xdata = [];
    Ydata = [];
    for iFig=1:2
        ax = axes{iFig}(iAxes);
        data = get( ax,'Children');
        Xdata = [Xdata cell2mat(get(data,'XData') )];
        Ydata = [Ydata cell2mat(get(data,'YData') )];
    end
    nLine = size(Xdata,1);
    X = zeros(size(Xdata));
    Y = zeros(size(Xdata));
    for iLine=nLine:-1:1
        [a indSort] = sort(Xdata(iLine,:));
        X(iLine,:) = a;
        Y(iLine,:) = Ydata(iLine,indSort);
        plot(h,X(iLine,:),Y(iLine,:));
        textX = get( ax,'Xlabel');
        xlabel(h,textX.String);
        textY = get( ax,'Ylabel');
        ylabel(h,textY.String);
    end
    if iAxes ==4
        leg = axes{1}(3);
        legend(leg.String);
        tit = get( ax,'Title');
        title(h,tit.String);
    end
   

    
end
% close all
