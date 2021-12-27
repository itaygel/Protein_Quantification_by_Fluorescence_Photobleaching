close all
 tic
format % long
clearvars -except  yt  peter_Swain_data dd

%figure(10)
hold on

namber_of_experiments=19;

for j=1:namber_of_experiments %123
    %for j=8:12
j

f=1;
nf=1;




yt=peter_Swain_data(j).yt;
%yt=round(yt,2);
%yt_big=yt;
%yt = zip_yt( yt,5 );



f_name=peter_Swain_data(j).name;
% start_frame=6; %123
% end_frame=36;%40;

h=yt(:,3:5)./yt(:,2:4);
p_test(j)=mean(h(:));
p_var_test(j)=var(h(:));




if j<7
    start_frame=2;
end_frame=40;%40;
end

if j>=7 && j<13
start_frame=2;
end_frame=20;%40;
end

if j>12 && j<19
start_frame=2;
end_frame=10;%40;
end

if j==19
start_frame=2;
end_frame=40;%40;
end

frame_num=end_frame-start_frame+1;
[cell_num,~]=size(yt);
x=zeros(cell_num,frame_num);
xx=zeros(cell_num,frame_num);
xxx=zeros(cell_num,frame_num-1);
xxxx=zeros(cell_num,frame_num-1);
ni=zeros(cell_num,frame_num-1);



for i=1:cell_num
    
    
    
    y=yt(i,start_frame:end_frame)';
   %y=yt_big(i,:)';
    %[fitresult, gof] = two_exp_fit(y);
    [fitresult, gof] = two_exp_fit_t0_fix(y);
    %y=yt(i,start_frame:end_frame)';
    if gof.rsquare>0.95
   
        %ey=feval(fitresult, 1:frame_num);
        ey=feval(fitresult, 1:frame_num);
        
%         if compare_to_prevus==1
%             p=ey(2:end)./ey(1:end-1);
%             ey(2:end)=y(1:end-1).*p;
%             
%          end   
            
            p=ey(2:end)./ey(1:end-1);
            ey_ctp=ey;      
            ey_ctp(2:end)=y(1:end-1).*p; 
            
       
    
    
%     y_noise_fix=zeros(1,length(y));
%     
%     for k=1:length(y)
%         if (y(k)-ey(k))^2>ey(k)
%         y_noise_fix(k)=(y(k)-ey(k))^2-ey(k);
%          f=f+1;
%          else
%             y_noise_fix(k)=(y(k)-ey(k))^2;
%             nf=nf+1;
%          end
%         
%     end

y_noise_fix=(y-ey).^2-ey;
%y_noise_fix_ctp=(y-ey).^2-2*ey;
y_noise_fix_ctp=(y-ey_ctp).^2-2*ey;





   %y_noise_fix=y_noise_fix';
    y_noise=(y-ey).^2;
    y_noise_ctp=(y-ey_ctp).^2;

    
       
    % y_noise_fix=y_noise_fix';
    
    
       % ni(i,1:frame_num-1)=y(1)./(y_noise_ctp(2:end)./(y(1:end-1).*p.*(1-p)));

    ni(i,1:frame_num-1)=y(1)./(y_noise_fix_ctp(2:end)./(y(1:end-1).*p.*(1-p)));
    
    
    
   x(i,1:frame_num)=ey.*(y(1)-ey)./y_noise;
   xx(i,1:frame_num)=ey.*(y(1)-ey)./y_noise_fix;
   %xxx(i,1:frame_num)=ey_ctp.*(y(1)-ey_ctp)./y_noise_ctp;
   %xxxx(i,1:frame_num)=ey.*(y(1)-ey_ctp)./y_noise_fix_pfc;
     xxx(i,1:frame_num-1)=(y(1)./y(2:end)).*ey_ctp(2:end).*(y(1:end-1)-ey_ctp(2:end))./y_noise_ctp(2:end);
   xxxx(i,1:frame_num-1)=(y(1)./y(2:end)).*ey_ctp(2:end).*(y(1:end-1)-ey_ctp(2:end))./y_noise_fix_ctp(2:end); 
    end
end

% f
% nf

% x(x<=0)=[];
% xx(xx<=0)=[];
x=real(x);
xx=real(xx);
xxx=real(xxx);
xxxx=real(xxxx);

xni=real(ni);
vv1=log10(real(xni(:)));
vv1=real(vv1);
vv1(vv1<-10000000000000)=[];



v=log10(real(x(:)));
v=real(v);
v(v<-10000000000000)=[];

subplot(2,round(namber_of_experiments/2),j) %123
%subplot(1,namber_of_experiments,j-7)
hold on
%a=histogram(real(v),'DisplayStyle','stairs','LineWidth',2,'Normalization','probability');
a=histogram(real(v),'LineWidth',0.5,'Normalization','probability');
mean(xx(:))
vv=log10(real(xx(:)));
vv=real(vv);
vv(vv<-10000000000000)=[];

%histogram(real(vv))


vvv=log10(real(xxx(:)));
vvv=real(vvv);
vvv(vvv<-10000000000000)=[];

%histogram(real(vvv),'DisplayStyle','stairs','LineWidth',2,'Normalization','probability')
histogram(real(vvv),'LineWidth',0.5,'Normalization','probability')
%plot_19ex_data(peter_Swain_data(j).ex19_plot,a.Values)
vvvv=log10(real(xxxx(:)));
vvvv=real(vvvv);
vvvv(vvvv<-10000000000000)=[];


%legend('standart','fix')
%title([f_name ' ' num2str(median(v(:))) ' ' num2str(median(vv(:)))])
title([f_name])

xlim([1 8])

plotyy(nan,nan,peter_Swain_data(j).ex19_plot(:,1),peter_Swain_data(j).ex19_plot(:,2))
if j==1
    xlabel('Log_{10}(number of molecules)')
    ylabel('Frequency')
end

res(j,1)=median(v(:));
res(j,2)=median(vv(:));
res(j,3)=var(v(:));
res(j,4)=var(vv(:));
res(j,5)=median(vvv(:));
res(j,6)=var(vvv(:));
res(j,7)=median(vvvv(:));
res(j,8)=var(vvvv(:));


res(j,9)=median(vv1(:));
res(j,10)=var(vv1(:));

end


for d=1:namber_of_experiments
    %t(d)=peter_Swain_data(d).EX_19_res;
    t(d)=peter_Swain_data(d).true_res;
end
toc
res(:,1)./t'
res(:,2)./t'
res(:,5)./t'
 toc
 
 res(:,1)-t'
res(:,2)-t'
res(:,5)-t'

%T = array2table(res,'VariableNames',{'ffc','ffcfix','ffcvar','ffc_fix_var','ffc_fix','ffc_var','dfw','thwsr'},'RowNames',{'sdhg','sdgds','sdg','brtg','sgs','bdt'})


%set(gcf,'PaperUnits','cent','PaperPosition',[0 0 18.5 12])
%print(['6subplot'], '-dtiff', '-r900')


figure(1001)

hold on
plot(t(1:6),res(1:6,1),'ob')
plot(t(1:6),res(1:6,5),'xr')
plot([3 6.5],[3 6.5],'k:')
errorbar(t(1:6),res(1:6,1),sqrt(res(1:6,3)),'b','LineStyle','none')
errorbar(t(1:6),res(1:6,5),sqrt(res(1:6,6)),'r','LineStyle','none')

legend('FFC','PFC')

xlabel('Biochemical experiments  log_{10}(number of molecules)')
ylabel('Bi-exponential fit log_{10}(number of molecules)')


figure(1011)

subplot(2,1,1)
%figure(1002)


hold on


%xlabel('Biochemical experiments')
ylabel('Estimates from bi-exponential fit')



plot(10.^(res(7:12,1))./10.^(res(1:6,1)),'*r')
plot(10.^(res(13:18,1))./10.^(res(1:6,1)),'*b')

plot([0,6.5],[1,1],'LineStyle','--','Color',[0 0 0])

%xlabel('Biochemical experiments Log_{10}(number of molecules)')
%ylabel('Estimates from bi-exponential fit Log_{10}(number of molecules)')
legend('2\mu','4\mu')



% figure(1002)
% hold on
% plot(10.^t(1:6),10.^res(1:6,1),'ok')
% plot(10.^t(1:6),10.^res(1:6,5),'xk')
% errorbar(10.^t(1:6),10.^res(1:6,1),sqrt(10.^res(1:6,3)),'LineStyle','none')
% errorbar(10.^t(1:6),10.^res(1:6,5),sqrt(10.^res(1:6,6)),'LineStyle','none')
subplot(2,1,2)


hold on


xlabel('Protein')
ylabel('Ratio')



plot(10.^(res(7:12,5))./10.^(res(1:6,5)),'*r')
plot(10.^(res(13:18,5))./10.^(res(1:6,5)),'*b')

plot([0,6.5],[1,1],'LineStyle','--','Color',[0 0 0])

legend('2\mu','4\mu')









figure(1012)
%subplot(1,2,1)

hold on
plot(t(1:6),res(1:6,1),'or')
%plot(t(1:6),res(1:6,2),'xr')
%plot([3 6.5],[3 6.5],'k')
%errorbar(t(1:6),res(1:6,2),sqrt(res(1:6,3)),'b','LineStyle','none')
%errorbar(t(1:6),res(1:6,7),sqrt(res(1:6,6)),'r','LineStyle','none')


xlabel('Biochemical experiments')
ylabel('Estimates from bi-exponential fit')




plot(t(1:6),res(1:6,2),'xb')
%plot(t(1:6),res(7:12,2),'xb')

%plot(t(1:6),res(13:18,2),'+k')
%plot(t(1:6),res(13:18,2),'xk')


%plot([2.5,6.5],[2.5,6.5],'LineStyle','--','Color',[0 0 0])
%legend('\mu','2\mu','4\mu')
% 

figure(1003)
% compre the fix noise
plot(res(1:6,2)./res(1:6,1),'or')
hold on
plot(res(1:6,7)./res(1:6,5),'xb')





figure(100)
for i=1:4
    subplot(4,1,i) 
    a=[1,5,2,7];
    b=a(i);
%plot(res(1:6,1)./res(14:19,1),'xb')
hold on
%plot(res(1:6,1)./res(8:13,1),'xk')
%plot(t(1:6)'./res(1:6,1),'*g')

%plot(res(7:12,b)./res(1:6,b),'*r')
%plot(res(13:18,b)./res(1:6,b),'*b')

plot(10.^(res(7:12,b))./10.^(res(1:6,b)),'*r')
plot(10.^(res(13:18,b))./10.^(res(1:6,b)),'*b')
legend('8-13','14-19')

end




%  


hold on
figure(7)
subplot(2,2,2)

plot(t(1:6),res(1:6,1),'or')
%plot(t(1:6),res(1:6,2),'xr')
%plot([3 6.5],[3 6.5],'k')
%errorbar(t(1:6),res(1:6,2),sqrt(res(1:6,3)),'b','LineStyle','none')
%errorbar(t(1:6),res(1:6,7),sqrt(res(1:6,6)),'r','LineStyle','none')


xlabel('Biochemical experiments')
ylabel('Estimates from bi-exponential fit')




plot(t(1:6),res(7:12,1),'xb')
%plot(t(1:6),res(7:12,2),'xb')

plot(t(1:6),res(13:18,1),'+k')
%plot(t(1:6),res(13:18,2),'xk')


plot([2.5,6.5],[2.5,6.5],'LineStyle','--','Color',[0 0 0])
%ylabel('log10(median number of molecules)')
% ylabel('Result from bi-exponential fit, log10(median number of molecules)')
% xlabel('Result from biochemical experiments, log10(median number of molecules)')

xlabel('Biochemical experiments Log_{10}(number of molecules)')
ylabel('Estimates from bi-exponential fit Log_{10}(number of molecules)')
legend('\mu','2\mu','4\mu')


subplot(2,2,3)
%figure(1002)
%subplot(1,2,1)

hold on
plot(t(1:6),res(1:6,5),'or')
%plot(t(1:6),res(1:6,2),'xr')
%plot([3 6.5],[3 6.5],'k')
%errorbar(t(1:6),res(1:6,2),sqrt(res(1:6,3)),'b','LineStyle','none')
%errorbar(t(1:6),res(1:6,7),sqrt(res(1:6,6)),'r','LineStyle','none')


xlabel('Biochemical experiments')
ylabel('Estimates from bi-exponential fit')




plot(t(1:6),res(7:12,5),'xb')
%plot(t(1:6),res(7:12,2),'xb')

plot(t(1:6),res(13:18,5),'+k')
%plot(t(1:6),res(13:18,2),'xk')


plot([2.5,6.5],[2.5,6.5],'LineStyle','--','Color',[0 0 0])
legend('\mu','2\mu','4\mu')
% 
for i=1:1
figure(777)
%subplot(4,1,1)
a=[1 2 5 7 9];



e=1;%a(i);
%subplot(5,1,i)
plot(t(1:6),res([1:6],e),'o')
hold on
plot(t(1:6),res([7:12],e),'o')
plot(t(1:6),res([13:18],e),'o')
%plot(t(1:6),t(1:6),'+')
plot([3.5 6],[3.5 6],':')

legend('\mu','2\mu','4\mu')
xlabel('Biochemical experiments')
ylabel('Gaussian Processes Estimating')
end