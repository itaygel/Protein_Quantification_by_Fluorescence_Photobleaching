
%מציג רעש יריה והתיקון שלו
% רעש רקע והתיקון שלו
%רעש יריה ורעש רקע והתיקון
xaxes=[100 1200];

close all
figure(11)
hold on

set(0,'DefaultFigureWindowStyle' , 'normal')
set(gcf,'PaperUnits','cent','PaperPosition',[0 0 23.5 11])

subplot(2,3,1)
load('ex2_no_fix_noise_bg_amp_0_noise_amp100.mat')

errorbar(v_m(2:2:length(v_m)),v_var_m(2:2:length(v_var_m)),sqrt(mean_var_m(2:2:length(mean_var_m))),'k','MarkerSize',4,'Marker','square','LineStyle','none')
%plot(v_m,v_var_m,'MarkerSize',4,'Marker','.')
hold on
errorbar(v_m(1:2:length(v_m)),v_var_m2(1:2:length(v_var_m2)),sqrt(mean_var_m2(1:2:length(mean_var_m2))),'r','MarkerSize',4,'Marker','diamond','LineStyle','none')

%plot(v_m,v_var_m2,'MarkerSize',4,'Marker','.')
int_start=damp;
int_end=damp^23; 




%int_limit=1/((int_start)^2/2-(int_start)^3/3-((int_end)^2/2-(int_end)^3/3));

%bg_error1=(1/(0.5-int_start/3))-(1/(0.5-int_end/3));
%bg_error1=int_limit*(int_start^2-int_end^2);

bg_error1=2*log((int_end-1)/(int_start-1));
plot(v_m,(v_m+bg_error1)./v_m,'r:','LineWidth',1)

%bg_error1=2*log((int_start)/(int_end))-2*log((1-int_start)/(1-int_end))+(int_start^2)/2-(int_end^2)/2;
% indepedend integral
% int_start=1;
% int_end=damp^23; 
% up=int_start^2/2-int_end^2/2;
% int_start=(damp+1)/2;
% int_end=damp^23; 
% down=log(int_start)-log(1-int_start)-(log(int_end)-log(1-int_end));
% bg_error1=2*up/down;


plot(v_m,(v_m+bg_error1)./v_m,'r:','LineWidth',1)


bg_error=2/(1-damp);
plot(v_m,(v_m+bg_error)./v_m,'k:','LineWidth',1)

xlim(xaxes);
ylim([0.88 1.25]);
set(gca, 'YGrid', 'on', 'XGrid', 'on')
ylabel({'\nu_{error}/\nu_{real}'},'fontsize',12);

subplot(2,3,4)



load('ex2_fix_noise_bg_amp_0_noise_amp100.mat')

errorbar(v_m(2:2:length(v_m)),v_var_m(2:2:length(v_var_m)),sqrt(mean_var_m(2:2:length(mean_var_m))),'k','MarkerSize',4,'Marker','square','LineStyle','none')
%plot(v_m,v_var_m,'MarkerSize',4,'Marker','.')
hold on
errorbar(v_m(1:2:length(v_m)),v_var_m2(1:2:length(v_var_m2)),sqrt(mean_var_m2(1:2:length(mean_var_m2))),'r','MarkerSize',4,'Marker','diamond','LineStyle','none')


%plot(v_m,v_var_m2,'MarkerSize',4,'Marker','.')

% Create xlabel
%xlabel({'\nu'},'fontsize',12);
ylim([0.88 1.25]);

% Create ylabel
ylabel({'\nu_{error}/\nu_{real}'},'fontsize',12);
set(gca, 'YGrid', 'on', 'XGrid', 'on')

xlabel({'\nu'},'fontsize',12);




xlim(xaxes);


subplot(2,3,2)



load('ex2_no_fix_noise_bg_amp_10_noise_amp0.mat')

errorbar(v_m(2:2:length(v_m)),v_var_m(2:2:length(v_var_m)),sqrt(mean_var_m(2:2:length(mean_var_m))),'k','MarkerSize',4,'Marker','square','LineStyle','none')
%plot(v_m,v_var_m,'MarkerSize',4,'Marker','.')
hold on
errorbar(v_m(1:2:length(v_m)),v_var_m2(1:2:length(v_var_m2)),sqrt(mean_var_m2(1:2:length(mean_var_m2))),'r','MarkerSize',4,'Marker','diamond','LineStyle','none')

%plot(v_m,v_var_m2,'MarkerSize',4,'Marker','.')




op_damp=1/damp;
bg_error=(noise_bg_amp*(1+damp)/(1-damp))*(op_damp^frame_num-1)/((op_damp-1)*frame_num);
% int_limit=1/((1-serviv)^2/2-(1-serviv)^3/3);
% bg_error1=((1-serviv)+(1-serviv)^2/2)*int_limit*noise_bg_amp;
%bg_error1=(1.5-(1-serviv)-(1-serviv)^2/2)*int_limit*noise_bg_amp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_start=damp;
int_end=serviv;

%int_limit=1/((int_start)^2/2-(int_start)^3/3-((int_end)^2/2-(int_end)^3/3));

%bg_error1=(int_start+int_start^2/2-(int_end+int_end^2/2))*noise_bg_amp*int_limit;
bg_error1=noise_bg_amp*log((int_start)/(int_end))-2*noise_bg_amp*log((1-int_start)/(1-int_end));
plot(v_m,(v_m+bg_error1)./v_m,'r:','LineWidth',1)
%   seperate integrals
%bg_error1=noise_bg_amp*log((int_start)/(int_end))-noise_bg_amp*log((1-int_start)/(1-int_end))+(int_start^2)/2-(int_end^2)/2+int_start-int_end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot(v_m,(v_m+bg_error)./v_m,'k:','LineWidth',1)
plot(v_m,(v_m+bg_error1)./v_m,'r:','LineWidth',1)


ylim([0.9 1.15]);
xlim(xaxes);

set(gca, 'YGrid', 'on', 'XGrid', 'on')

subplot(2,3,5)

load('ex2_fix_noise_bg_amp_10_noise_amp0.mat')

errorbar(v_m(2:2:length(v_m)),v_var_m(2:2:length(v_var_m)),sqrt(mean_var_m(2:2:length(mean_var_m))),'k','MarkerSize',4,'Marker','square','LineStyle','none')
%plot(v_m,v_var_m,'MarkerSize',4,'Marker','.')
hold on
errorbar(v_m(1:2:length(v_m)),v_var_m2(1:2:length(v_var_m2)),sqrt(mean_var_m2(1:2:length(mean_var_m2))),'r','MarkerSize',4,'Marker','diamond','LineStyle','none')

%plot(v_m,v_var_m2,'MarkerSize',4,'Marker','.')
xlim(xaxes);
ylim([0.9 1.15]);


% Create xlabel
%xlabel({'\nu'},'fontsize',12);

% Create ylabel
set(gca, 'YGrid', 'on', 'XGrid', 'on')
xlabel({'\nu'},'fontsize',12);






subplot(2,3,3)


set(gca, 'YGrid', 'on', 'XGrid', 'on')

load('ex2_no_fix_noise_bg_amp_10_noise_amp100.mat')

errorbar(v_m(2:2:length(v_m)),v_var_m(2:2:length(v_var_m)),sqrt(mean_var_m(2:2:length(mean_var_m))),'k','MarkerSize',4,'Marker','square','LineStyle','none')
%plot(v_m,v_var_m,'MarkerSize',4,'Marker','.')
hold on
errorbar(v_m(1:2:length(v_m)),v_var_m2(1:2:length(v_var_m2)),sqrt(mean_var_m2(1:2:length(mean_var_m2))),'r','MarkerSize',4,'Marker','diamond','LineStyle','none')

ylim([0.85 1.35]);

xlim(xaxes);


int_start=damp;
int_end=serviv;

bg_error1=-noise_bg_amp*log((int_end)/(int_start))+2*noise_bg_amp*log((1-int_end)/(1-int_start))+2*log((int_end-1)/(int_start-1));


plot(v_m,(v_m+bg_error1)./v_m,'r:','LineWidth',1)

op_damp=1/damp;
bg_error=2/(1-damp)+(noise_bg_amp*(1+damp)/(1-damp))*(op_damp^frame_num-1)/((op_damp-1)*frame_num);

plot(v_m,(v_m+bg_error)./v_m,'k:','LineWidth',1)

set(gca, 'YGrid', 'on', 'XGrid', 'on')

%plot(v_m,v_var_m2,'MarkerSize',4,'Marker','.')
subplot(2,3,6)
load('ex2_fix_noise_bg_amp_10_noise_amp100.mat')

errorbar(v_m(2:2:length(v_m)),v_var_m(2:2:length(v_var_m)),sqrt(mean_var_m(2:2:length(mean_var_m))),'k','MarkerSize',4,'Marker','square','LineStyle','none')
%plot(v_m,v_var_m,'MarkerSize',4,'Marker','.')
hold on
errorbar(v_m(1:2:length(v_m)),v_var_m2(1:2:length(v_var_m2)),sqrt(mean_var_m2(1:2:length(mean_var_m2))),'r','MarkerSize',4,'Marker','diamond','LineStyle','none')








%plot(v_m,v_var_m2,'MarkerSize',4,'Marker','.')
xlim(xaxes);
ylim([0.85 1.35]);

% Create xlabel
xlabel({'\nu'},'fontsize',12);

% Create ylabel
set(gca, 'YGrid', 'on', 'XGrid', 'on')










%set(gcf,'PaperUnits','cent','PaperPosition',[0 0 8.25 10])

set(gcf,'PaperUnits','cent','PaperPosition',[0 0 17.3 10])
%print(['my_figure1'], '-dtiff', '-r900')