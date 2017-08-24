n=1000;
%y=rand(n,1);
m=3;
Parament=cell(500,2);
Expectation=zeros(500,1);

SP=xlsread('J:\×ÀÃæ\Project\^GSPC.csv');
SP=SP(:,end-1);
rSP=(SP(2:end)-SP(1:end-1))./SP(1:end-1);
y=SP;

x=zeros(3,1);
x(1)=sum(y<=prctile(y,33));
x(3)=sum(y>prctile(y,66));
x(2)=length(y)-x(1)-x(3);

for em=1:500
P=rand(m);
for i=1:m
   P(i,:)=P(i,:)./sum(P(i,:)',1);
end

mu=rand(1,m);
%mu=mu./sum(mu);
sigma=abs(rand(1,m))*5;
%sigma=sigma./sum(sigma);

lmq=1;

while lmq<1000


Ptem=x/sum(x); 

   
alpha=zeros(n+1,m);
alpha(1,:)=Ptem;



for t=1:n
    for i=1:m
        sum0=0;
        for j=1:m 
            sum0=sum0+alpha(t,j)*P(j,i)*(1/(y(t+1)*((2*pi*sigma(i)^2)^0.5))*exp(-(log(y(t+1))-mu(i))^2/(2*sigma(i)^2)));
             %sum0=sum0+alpha(t,j)*P(j,i)*(1/(1*((2*pi*sigma(i)^2)^0.5))*exp(-(log(y(t+1))-log(y(t)))^2/(2*sigma(i)^2)));
        end
        alpha(t+1,i)=sum0;
    end
end

beta=zeros(n,m);
beta(n,:)=1;
for t=n-1:-1:1
   for i=1:m
       sum0=0;
      for j=1:m
      sum0=sum0+beta(t+1,j)*P(i,j)*(1/((2*pi*sigma(j)^2)^0.5)*exp(-(log(y(t+1))-mu(j))^2/(2*sigma(j)^2)));
      end
      beta(t,i)=sum0;
   end
end
Ptem2=P;
for i=1:m
    for j=1:m
        sum1=0;
        for t=2:n
           sum1=sum1+alpha(t,i)*P(i,j)*(1/((2*pi*sigma(j)^2)^0.5)*exp(-(log(y(t))-mu(j))^2/(2*sigma(j)^2)))*beta(t,j);
        end
        sum2=0;
        for j1=1:m
            for t1=2:n
                sum2=sum2+alpha(t1,i)*P(i,j1)*(1/((2*pi*sigma(j1)^2)^0.5)*exp(-(log(y(t1))-mu(j1))^2/(2*sigma(j1)^2)))*beta(t1,j1);
            end
        end
        Ptem2(i,j)=sum1/sum2;
    end
end

P=Ptem2;
mu_new=mu;
for i=1:m
    sum1=0;
    sum2=0;
    for t=1:n
        sum1=sum1+log(y(t))*alpha(t+1,i)*beta(t,i);
    end
    for t=1:n
        sum2=sum2+alpha(t+1,i)*beta(t,i);
    end
mu_new(i)=sum1/sum2;
sum3=0;
    for t=1:n
        sum3=sum3+(log(y(t))-mu_new(i))^2*alpha(t+1,i)*beta(t,i);
    end
sigma_new=(sum3/sum2)^0.5;
end

error=sum(abs((mu_new-mu).*(sigma_new-sigma)));
if error<0.1
     break
else
    mu=mu_new;
    sigma=sigma_new;
end
lmq=lmq+1;
end
tParament(1,:)=mu';
tParament(2,:)=sigma';
Parament{em,1}=tParament;
Parament{em,2}=alpha;
Expectation=sum(log((1/((2*pi.*sigma^2)^0.5).*exp(-(log(y)-mu).^2./(2.*sigma.^2)))));
end
Select_paramentor=Parament{Expectation==max(Expectation),1};
Select_alpha=Parament{Expectation==max(Expectation),2};
Pi=zeros(length(Select_alpha(:,1)),length(Select_alpha(1,:)));
for t=1:length(Select_alpha(:,1))
    Pi(t,:)=Select_alpha(t,:)/sum(Select_alpha(t,:));
end

Time=1:length(rSP);
X1=Time(Pi(:,1)>0.95);
X3=Time(Pi(:,3)>0.95);
X2=Time(Pi(:,3)<=0.95&&Pi(:,1)<=0.95);
X1(X1==1)=[];
X2(X2==1)=[];
X3(X3==1)=[];
X1_Period=[];
for i1=1:length(X1)
    X1_Period=[X1_Period,X1(i1)-1];
    X1_Period=[X1_Period,X1(i1)];
    X1_Period=[X1_Period,X1(i1)+1];
    X1_Period=[X1_Period,X1(i1)+2];
    X1_Period=[X1_Period,X1(i1)+3];
end

X2_Period=[];
for i1=1:length(X1)
    X2_Period=[X2_Period,X2(i1)-1];
    X2_Period=[X2_Period,X2(i1)];
    X2_Period=[X2_Period,X2(i1)+1];
    X2_Period=[X2_Period,X2(i1)+2];
    X2_Period=[X2_Period,X2(i1)+3];
end

X3_Period=[];
for i1=1:length(X1)
    X3_Period=[X3_Period,X3(i1)-1];
    X3_Period=[X3_Period,X3(i1)];
    X3_Period=[X3_Period,X3(i1)+1];
    X3_Period=[X3_Period,X3(i1)+2];
    X3_Period=[X3_Period,X3(i1)+3];
end

X1_Period=unique(X1_Period);
X2_Period=unique(X2_Period);
X3_Period=unique(X3_Period);



