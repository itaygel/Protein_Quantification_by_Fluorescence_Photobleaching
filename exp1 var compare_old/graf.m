%v_m=v_m/100;
close all
figure(1)
set(gcf,'PaperUnits','cent','PaperPosition',[0 0 8 6])

hold on
t=1:5:40;

tt=3:5:40;
errorbar(v_m(t),v_var_m(t),mean_var_m(t))
errorbar(v_m(tt),v_var_m2(t),mean_var_m2(t))
xlabel('Percentage of survivors [%]','fontsize',12)
ylabel('\nu^2_{est}/\nu^2','fontsize',12)
plot(v_m(t),v_var_m(t),'rx','MarkerSize',5)
plot(v_m(tt),v_var_m2(tt),'b+','MarkerSize',5)
plot(v_m,1+(2-6./1000+1./(1000*p*(1-p)))./(frame_numa_mean-1),':k')
plot(v_m,1-(2-6./1000+1./(1000*p*(1-p)))./(frame_numa_mean-1),':k')

plot([0 20],[1 1])


figure(2)
hold on
set(gcf,'PaperUnits','cent','PaperPosition',[0 0 8.5 6])


plot(v_m(t),1+mean_var_m(t),'.r','MarkerSize',8)
plot(v_m(tt),1+mean_var_m2(tt),'.b','MarkerSize',8)
%plot(v_m,1+mean_var_m2_artical_res_no_integral,'.k','MarkerSize',8)

plot(v_m(t),v_var_m(t),'rx','MarkerSize',5)


plot(v_m(t),1-mean_var_m(t),'.r','MarkerSize',8)


plot(v_m(tt),v_var_m2(tt),'b+','MarkerSize',5)
hold on

plot(v_m(tt),1-mean_var_m2(tt),'.b','MarkerSize',8)


%plot(v_m,v_var_m2_artical_res_no_integral,'g')
hold on

%plot(v_m,1-mean_var_m2_artical_res_no_integral,'.k','MarkerSize',8)

xlabel('Percentage of survivors [%]','fontsize',12)
ylabel('\nu^2_{est}/\nu^2','fontsize',12)
%legend('uncorelated','corelated','corelated no integral')


plot(v_m,1+(2-6./1000+1./(1000*p*(1-p)))./(frame_numa_mean-1),'k')
plot(v_m,1-(2-6./1000+1./(1000*p*(1-p)))./(frame_numa_mean-1),'k')

%set(gcf,'PaperUnits','cent','PaperPosition',[0 0 8.5 11.5])

% print(['my_figure1'], '-dtiff', '-r600')