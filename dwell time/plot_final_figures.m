% read results and final plot
swi4 = load('Swi4_bound_time_results.mat').results;
Whi5 = load('Whi5_bound_time_results.mat').results;
YHZ68 = load("YHZ68_bound_time_results.mat").results;
YHZ75 = load('YHZ75_bound_time_results.mat').results;

figure 
plotboundtime('Swi4-Halo',swi4{1,1},swi4{3,1},swi4{4,1});

figure
plotboundtime('Whi5-Halo',Whi5{1,1},Whi5{3,1},Whi5{4,1});

figure
plotboundtime('Swi4-Halo in Whi5 overexpression',YHZ68{1,1},YHZ68{3,1},YHZ68{4,1});

figure
plotboundtime('Swi6-Halo in Mbp1-deletion',YHZ75{1,1},YHZ75{3,1},YHZ75{4,1});

figure
plotboundtime2('Pre-Start daughter cells ',Whi5{2,1},swi4{2,1},YHZ68{2,1},YHZ75{2,1});

figure
plotboundtime2('Mother cells',Whi5{1,1},swi4{1,1},YHZ68{1,1},YHZ75{1,1});



















function plotboundtime(strainnname,boundtime_extr1,boundtime_extr3,boundtime_extr4)

plot_log_exp(boundtime_extr1,'r');
hold on;
plot_log_exp(boundtime_extr3,'g');
hold on;
plot_log_exp(boundtime_extr4,'b');
axis([2 15 -7 0])
legend('Mother cells','Small pre-Start daughter cells','Large pre-Start daughter cells','Location','southwest')
title(strainnname);

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

end

function plotboundtime2(strainnname,boundtime_extr1,boundtime_extr2,boundtime_extr3,boundtime_extr4)

plot_log_exp(boundtime_extr1,'r');
hold on;
plot_log_exp(boundtime_extr2,'k');
hold on;
plot_log_exp(boundtime_extr3,'g');
hold on;
plot_log_exp(boundtime_extr4,'b');
axis([2 15 -7 0])
legend('Whi5-Halo','Swi4-Halo','Swi4-Halo in Whi5-overexpression','Swi6-Halo in Mbp1-deletion','Location','southwest')
title(strainnname);

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

end