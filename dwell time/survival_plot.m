% swi4 bounnd time compare 
% finalize the datasets 
close all;
plot_log_exp(boundtime_extr1,'r');
hold on;
plot_log_exp(boundtime_extr3,'g');
hold on;
plot_log_exp(boundtime_extr4,'b');
axis([2 15 -7 0])
legend('Mother cells','Daughter cells (10fL -- 35fL )','Daughter cells (35fL -- 60fL )','Location','southwest')
title('Whi5-Halo');

function adjust_fig
set(gca, ...
     'Box'         , 'off'     , ...
     'TickDir'     , 'out'     , ...
     'XMinorTick'  , 'off'      , ...
     'YMinorTick'  , 'off'      , ...
      'YGrid'       , 'off',...
      'XGrid'       , 'off');
set(gca,'LineWidth',3,'FontSize',16,'box','off');
end


function plot_log_exp(datasets,datacolor)

[f,x] = ecdf((datasets),'Function','survivor','Alpha',0.01,'Bounds','on');
f2 = log(f);
aaa = find(f2 == 0);
x(aaa) = [];
f2(aaa) = [];
plot(x,f2,'.-','Color',datacolor,'MarkerSize',15,'LineWidth',2);
% hold on 
% plot(x,log(f),'-','Color',datacolor,'LineWidth',2);

adjust_fig
xlabel('Track duration (sec)');
ylabel('Survival probability (Log)');

set(gcf,'position',[200 200 380 320])

end