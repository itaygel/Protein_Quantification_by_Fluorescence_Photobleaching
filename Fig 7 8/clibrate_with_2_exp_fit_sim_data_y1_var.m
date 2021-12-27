%close all
format % long
clearvars -except  yt  peter_Swain_data dd
tic
fit_for_all=1

namber_of_experiments=21;

for j=1:8%namber_of_experiments
   
j

f=1;
nf=1;

small_rmse=0;
hige_rmse=0;

nnn(j)=round(1000*10^(j/6));
%[ yt,frame_num ] = build_flur_matrix_p_var(45,5000,nnn(j),1,false,0,750,0 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ yt1,frame_num ] = build_flur_matrix(45,3000,round(0.7*nnn(j)),1-0.85,false,0,750,0 );
[ yt2,frame_num ] = build_flur_matrix(45,3000,nnn(j)-round(0.7*nnn(j)),1-0.95,false,0,750,0 );

yt=yt1+yt2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







yt=yt';
%yt=peter_Swain_data(j).yt;
f_name='sim data';


h=yt(:,3:5)./yt(:,2:4);
p_test(j)=mean(h(:));
p_var_test(j)=var(h(:));


stop_frame=0;


    start_frame=1;
end_frame=45-stop_frame;%40;
f=1;
y0=yt(1);




%start_frame=1;
frame_num=end_frame-start_frame+1;
[cell_num,~]=size(yt);
x=zeros(cell_num,frame_num);
xx=zeros(cell_num,frame_num);
xxx=zeros(cell_num,frame_num-1);
xxxx=zeros(cell_num,frame_num-1);
ni=zeros(cell_num,frame_num-1);


if fit_for_all==1
    y1=yt(:,1);
       % [fitresult, gof] =two_exp_fit(mean(yt));
       [fitresult, gof] =two_exp_fit(1000*mean(yt./y1));
    end

for i=1:cell_num
    
    
    
    y=yt(i,start_frame:end_frame)';
   %y=yt_big(i,:)';
    %[fitresult, gof] = two_exp_fit(y);
    %y1=y(1);
    
    if fit_for_all==0
    [fitresult, gof] = two_exp_fit(y);
    %y=yt(i,start_frame:end_frame)';
   % if gof.rsquare>0.95
    end
   rmse_all(i,j)=gof.rmse;
   
   fitresult_f(i,j)=fitresult.f;
   fitresult_ab(i,j)=fitresult.a+fitresult.c;
   fitresult_abc(i,j)=fitresult.f/(fitresult.a+fitresult.c);
   
  
   if gof.adjrsquare>0.95
   %if gof.rmse< all_rmse(j)*100   %800*f  %4000
        %ey=feval(fitresult, 1:frame_num);
        ey=y1(i)*feval(fitresult, 1:frame_num)/1000 -fitresult.f;
        %ey=ey*y1/10;
        
        y=y-fitresult.f;
        
        
%         if compare_to_prevus==1
%             p=ey(2:end)./ey(1:end-1);
%             ey(2:end)=y(1:end-1).*p;
%             
%          end   
            
            p=ey(2:end)./ey(1:end-1);
            ey_ctp=ey;      
            ey_ctp(2:end)=y(1:end-1).*p; 
  
            
  
   
   %close all
       
%       figure(2000)
% hold off
%  plot(ey)
%  hold on
%  plot(y)
%  xlabel(num2str(i))
%  ylabel(num2str(gof.rmse))
%  
 
 
 
 
 
 
    
    f1(j).rsquare(i)=gof.rsquare;
    f1(j).rmse(i)=gof.rmse;
    f1(j).sse(i)=gof.sse;
    
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

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    m_fix=(1-1/length(y));
%  
%   y_noise_fix=y_noise_fix/m_fix;  
%      y_noise_fix_ctp=y_noise_fix_ctp/m_fix;  
%   y_noise=y_noise/m_fix;  
%   y_noise_ctp=y_noise_ctp/m_fix;  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % y_noise_fix=y_noise_fix';
    
    
       % ni(i,1:frame_num-1)=y(1)./(y_noise_ctp(2:end)./(y(1:end-1).*p.*(1-p)));

    ni(i,1:frame_num-1)=y0./(y_noise_fix_ctp(2:end)./(y(1:end-1).*p.*(1-p)));
    
    %%% x=ffc
    %%% xx=ffc fix
    % % %  xxx=pfc
    % % %  xxxx=pfc fix
    %%%    ni= rotenberg fotmula fix
    
   x(i,1:frame_num)=(y0/y(1))*ey.*(y(1)-ey)./y_noise;
   xx(i,1:frame_num)=(y0/y(1))*ey.*(y(1)-ey)./y_noise_fix;
   %xxx(i,1:frame_num)=ey_ctp.*(y(1)-ey_ctp)./y_noise_ctp;
   %xxxx(i,1:frame_num)=ey.*(y(1)-ey_ctp)./y_noise_fix_pfc;
     xxx(i,1:frame_num-1)=(y0./y(2:end)).*ey_ctp(2:end).*(y(1:end-1)-ey_ctp(2:end))./y_noise_ctp(2:end);
   xxxx(i,1:frame_num-1)=(y0./y(2:end)).*ey_ctp(2:end).*(y(1:end-1)-ey_ctp(2:end))./y_noise_fix_ctp(2:end); 
   
   
          hige_rmse=hige_rmse+1;

   else
       
       small_rmse=small_rmse+1;
       i
       j
       gof.adjrsquare
%        if gof.adjrsquare>0.998
%    figure()
%    hold on
%    plot(y)
%    plot(ey)
%        end
       
                 
             figure(2636)
             ey=feval(fitresult, 1:frame_num);
   hold off
   plot(y)
   hold on
   plot(ey)
   xlabel(num2str(gof.adjrsquare))
    end
end

    small_rmsel(j)= small_rmse;
    hige_rmsel(j)= hige_rmse;





figure()
hold on
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

x=x(:,2:end);
v=log10(real(x(:)));
v=real(v);
v(v<0)=[];







%subplot(3,round(namber_of_experiments/3),j) %123
%subplot(1,namber_of_experiments,j-7)
hold on
%a=histogram(real(v),'DisplayStyle','stairs','LineWidth',2,'Normalization','probability');
%a=histogram(real(v),'LineWidth',0.5,'Normalization','probability');

vv=log10(real(xx(:)));
vv=real(vv);
vv(vv<-10000000000000)=[];

hist_max(j,1) = max_hist(v);
hist_max(j,2) = log10(nnn(j));

%histogram(real(vv))


vvv=log10(real(xxx(:)));
vvv=real(vvv);
vvv(vvv<-10000000000000)=[];

%histogram(real(vvv),'DisplayStyle','stairs','LineWidth',2,'Normalization','probability')
%histogram(real(vvv),'LineWidth',0.5,'Normalization','probability')
%plot_19ex_data(peter_Swain_data(j).ex19_plot,a.Values)
vvvv=log10(real(xxxx(:)));
vvvv=real(vvvv);
vvvv(vvvv<-10000000000000)=[];


%legend('standart','fix')
%title([f_name ' ' num2str(median(v(:))) ' ' num2str(median(vv(:)))])


xlim([1 8])

%plotyy(nan,nan,peter_Swain_data(j).ex19_plot(:,1),peter_Swain_data(j).ex19_plot(:,2))
if j==1
    xlabel('Log_{10}(number of molecules)')
    ylabel('Frequency')
end
vvv(vvv>9999999999999999)=[];


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

hist_max
   %%% x=ffc  1
    %%% xx=ffc fix   2
    % % %  xxx=pfc    5
    % % %  xxxx=pfc fix   7
    %%%    ni= rotenberg fotmula fix  (vv1)   9
end


% res(:,1)./t'
% res(:,2)./t'
% res(:,5)./t'
%  
%  
%  res(:,1)-t'
% res(:,2)-t'
% res(:,5)-t'

%T = array2table(res,'VariableNames',{'ffc','ffcfix','ffcvar','ffc_fix_var','ffc_fix','ffc_var','dfw','thwsr'},'RowNames',{'sdhg','sdgds','sdg','brtg','sgs','bdt'})


%set(gcf,'PaperUnits','cent','PaperPosition',[0 0 18.5 12])
%print(['6subplot'], '-dtiff', '-r900')

%figure(1002)








figure(502)
plot(hist_max(:,2),hist_max(:,1),'ok')
hold on
plot([2 7],[2 7])
plot(hist_max(:,2),res(:,1),'xr')
save
res