%x=
format long


%p=
clearvars -except  yt_all  p_all_all


for j=1:length(yt_all)
p=p_all_all(j,:);
yt=yt_all(j,:);
p_all(1)=exp(p(2));

for i=2:44
   
    p_all(i)=p_all(i-1)*exp(p(i+1));

    
    
    
end

Ey=yt(1)*p_all;
vary=(Ey-yt(2:end)).^2;

x0=Ey.*(yt(1)-Ey)./vary;
%m(j)=mean(x0);
%ma(j,1:length(yt_all)-1)=(x0);
1:length(yt_all)
end



close all
histogram(ma(:))
mean(m)
mm = m(~isnan(m));


mm(mm>mean(mm)+2*sqrt(var(mm)))=[];
mm(mm<0)=[];
mean(mm)
hold on

histogram(mm)