function nsga2_eval_perf(pf)

N=size(pf,1);
dim=size(pf,2);
d=zeros(N,1);
m=zeros(N,1);
s=zeros(N,1);
min_f=min(pf);
max_f=max(pf);
std_f=std(pf);
sum2=0;
sum1=zeros(N,1);

for i=1:N
    for j=1:dim
        d(i)=d(i) + ((pf(i,j)- min_f(j))/(max_f(j)-min_f(j)))^2;
        s(i)=s(i) + ((pf(i,j)- min_f(j))/std_f(j))^2;
        if pf(i,j)==min_f(j)
            mu1=1;
        else if min_f(j)<pf(i,j) && pf(i,j)<max_f(j)
                mu1=(max_f(j)-pf(i,j))/(max_f(j)-min_f(j));
            else
                mu1=0;
            end
        end
        sum1(i)=sum1(i)+mu1;
        sum2=sum2+mu1;
    end
    m(i)=max(abs(pf(i,:)-min_f));
end

mu=sum1/sum2;
[min_d,indx_d]=min(d);
[min_m,indx_m]=min(m);
[min_s,indx_s]=min(s);
[max_mu,indx_mu]=max(mu);

disp(min_d)
disp(indx_d);
disp(min_m);
disp(indx_m);
disp(min_s);
disp(indx_s);
disp(max_mu);
disp(indx_mu);

end