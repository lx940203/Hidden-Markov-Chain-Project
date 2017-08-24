
%y=rand(n,1);
m=3;
Parament=cell(500,2);
Expectation=zeros(500,1);

SP=xlsread('J:\×ÀÃæ\Project\^GSPC.csv');
SP=SP(:,end-1);
rBAC2=(SP(2:end)-SP(1:end-1))./SP(1:end-1);
y=SP;
n=length(y);
x=zeros(3,1);
x(1)=sum(y<=prctile(y,33));
x(3)=sum(y>prctile(y,66));
x(2)=length(y)-x(1)-x(3);

for em=1:500
P=rand(m);
for i=1:m
   P(i,:)=P(i,:)./sum(P(i,:)',1);
end

mu=rand(1,m)*0.01+mean(log(y));
%mu=mu./sum(mu);
sigma=abs(rand(1,m))*100000;
%sigma=sigma./sum(sigma);

lmq=1;

while lmq<1000


Ptem=x/sum(x); 

   
alpha=zeros(n+1,m);
alpha(1,:)=Ptem;



for t=1:n-1
    for i=1:m
        sum0=0;
        for j=1:m 
            sum0=sum0+alpha(t,j)*P(j,i)*(1/(y(t+1)*((2*pi*sigma(i)^2)^0.5))*exp(-(log(y(t+1))-mu(i))^2/(2*sigma(i)^2)));
            %sum0=sum0+alpha(t,j)*P(j,i)*(1/(1*((2*pi*sigma(i)^2)^0.5))*exp(-(log(y(t+1))-log(y(t+1)))^2/(2*sigma(i)^2)));
        end
        alpha(t+1,i)=sum0;
    end
    alpha(t+1,:)=alpha(t+1,:)/sum(alpha(t+1,:));
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
   beta(t,:)=beta(t,:)/sum(beta(t,:));
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
    sum1(i)=0;
    sum2(i)=0;
    for t=1:n
        sum1(i)=sum1(i)+log(y(t))*alpha(t+1,i)*beta(t,i);
    end
    for t=1:n
        sum2(i)=sum2(i)+alpha(t+1,i)*beta(t,i);
    end
mu_new(i)=sum1(i)/sum2(i);
sum3(i)=0;
    for t=1:n
        sum3(i)=sum3(i)+(log(y(t))-mu_new(i))^2*alpha(t+1,i)*beta(t,i);
    end
sigma_new(i)=(sum3(i)/sum2(i))^0.5;
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
Expectation(em,1)=sum(sum((1./((2*pi.*sigma.^2).^0.5).*exp(-(log(y)-mu).^2./(2.*sigma.^2)))));
end
Select_paramentor=Parament{Expectation==max(Expectation),1};
Select_alpha=Parament{Expectation==max(Expectation),2};
Pi=zeros(length(Select_alpha(:,1)),length(Select_alpha(1,:)));
for t=1:length(Select_alpha(:,1))
    Pi(t,:)=Select_alpha(t,:)/sum(Select_alpha(t,:));
end
Pi=(rand(length(rBAC2),3)-0.5)*5+0.5;
Pi=exp(Pi)./(1+exp(Pi));
Time=1:length(rBAC2);
X1=Time(Pi(:,1)>0.95);
X3=Time(Pi(:,3)>0.95);

for p=1:length(rBAC2)
    ptem(p)=Pi(p,3)<=0.95&&Pi(p,1)<=0.95;
end
X2=Time(ptem);
X1(X1==1)=[];
X2(X2==1)=[];
X3(X3==1)=[];
X1(X1==4401)=[];
X2(X2==4401)=[];
X3(X3==4401)=[];
X1(X1==4400)=[];
X2(X2==4400)=[];
X3(X3==4400)=[];
X1(X1==4399)=[];
X2(X2==4399)=[];
X3(X3==4399)=[];
X1(X1==4398)=[];
X2(X2==4398)=[];
X3(X3==4398)=[];
X1_Period=[];
for i1=1:length(X1)
    X1_Period=[X1_Period,X1(i1)-1];
    X1_Period=[X1_Period,X1(i1)];
    X1_Period=[X1_Period,X1(i1)+1];
    X1_Period=[X1_Period,X1(i1)+2];
    X1_Period=[X1_Period,X1(i1)+3];
end

X2_Period=[];
for i1=1:length(X2)
    X2_Period=[X2_Period,X2(i1)-1];
    X2_Period=[X2_Period,X2(i1)];
    X2_Period=[X2_Period,X2(i1)+1];
    X2_Period=[X2_Period,X2(i1)+2];
    X2_Period=[X2_Period,X2(i1)+3];
end

X3_Period=[];
for i1=1:length(X3)
    X3_Period=[X3_Period,X3(i1)-1];
    X3_Period=[X3_Period,X3(i1)];
    X3_Period=[X3_Period,X3(i1)+1];
    X3_Period=[X3_Period,X3(i1)+2];
    X3_Period=[X3_Period,X3(i1)+3];
end

X1_Period=unique(X1_Period);
X2_Period=unique(X2_Period);
X3_Period=unique(X3_Period);

BAC=xlsread('J:\×ÀÃæ\Project\BAC.csv');
BAC=BAC(:,end-1);
rBAC=(BAC(2:end)-BAC(1:end-1))./BAC(1:end-1);
rBAC1=rBAC(X1_Period);
rBAC2=rBAC(X2_Period);
rBAC3=rBAC(X3_Period);

%rBAC2=(SP(2:end)-SP(1:end-1))./SP(1:end-1);

%rBAC3=rBAC3(1:end-1,end-1);
rSPs=sort(rBAC2);
rBACs=sort(rBAC1);
VIXs=sort(rBAC3);
SP_DtI=rSPs(floor(length(rBAC2)/3));
SP_ItU=rSPs(floor(length(rBAC2)/3*2));
BAC_DtI=rBACs(floor(length(rBAC1)/3));
BAC_ItU=rBACs(floor(length(rBAC1)/3*2));
VIX_DtI=VIXs(floor(length(rBAC3)/3));
VIX_ItU=VIXs(floor(length(rBAC3)/3*2));

m=3;
rBAC1=[rBAC1,zeros(length(rBAC1),1),zeros(length(rBAC1),1),zeros(length(rBAC1),1)];
rBAC2=[rBAC2,zeros(length(rBAC2),1),zeros(length(rBAC2),1),zeros(length(rBAC2),1)];
rBAC3=[rBAC3,zeros(length(rBAC3),1),zeros(length(rBAC3),1),zeros(length(rBAC3),1)];

TrMBAC=zeros(m);
TrMSP=zeros(m);
TrMVIX=zeros(m);

for i=1:length(rBAC1)
    if rBAC1(i,1)>BAC_ItU
        rBAC1(i,2)=1;
    elseif rBAC1(i,1)<BAC_DtI
        rBAC1(i,2)=-1;
    else rBAC1(i,2)=0;
    end;
end

for i=1:length(rBAC2)
    if rBAC2(i,1)>SP_ItU
        rBAC2(i,2)=1;
    elseif rBAC2(i,1)<SP_DtI
        rBAC2(i,2)=-1;
    else rBAC2(i,2)=0;
    end;
end

for i=1:length(rBAC3)
    if rBAC3(i,1)>VIX_ItU
        rBAC3(i,2)=1;
    elseif rBAC3(i,1)<VIX_DtI
        rBAC3(i,2)=-1;
    else rBAC3(i,2)=0;
    end;
end

for i=1:length(rBAC1)-1
    if rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1
        rBAC1(i,3)=11;
        TrMBAC(1,1)=TrMBAC(1,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0
        rBAC1(i,3)=12;
        TrMBAC(1,2)=TrMBAC(1,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1
        rBAC1(i,3)=13;
        TrMBAC(1,3)=TrMBAC(1,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1
        rBAC1(i,3)=21;
        TrMBAC(2,1)=TrMBAC(2,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0
        rBAC1(i,3)=22;
        TrMBAC(2,2)=TrMBAC(2,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1
        rBAC1(i,3)=23;
        TrMBAC(2,3)=TrMBAC(2,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1
        rBAC1(i,3)=31;
        TrMBAC(3,1)=TrMBAC(3,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0
        rBAC1(i,3)=32;
        TrMBAC(3,2)=TrMBAC(3,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1
        rBAC1(i,3)=33;
    TrMBAC(3,3)=TrMBAC(3,3)+1;
    end;
end

for i=1:length(rBAC2)-1
    if rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1
        rBAC2(i,3)=11;
        TrMSP(1,1)=TrMSP(1,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0
        rBAC2(i,3)=12;
        TrMSP(1,2)=TrMSP(1,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1
        rBAC2(i,3)=13;
        TrMSP(1,3)=TrMSP(1,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1
        rBAC2(i,3)=21;
        TrMSP(2,1)=TrMSP(2,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0
        rBAC2(i,3)=22;
        TrMSP(2,2)=TrMSP(2,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1
        rBAC2(i,3)=23;
        TrMSP(2,3)=TrMSP(2,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1
        rBAC2(i,3)=31;
        TrMSP(3,1)=TrMSP(3,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0
        rBAC2(i,3)=32;
        TrMSP(3,2)=TrMSP(3,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==01
        rBAC2(i,3)=33;
    TrMSP(3,3)=TrMSP(3,3)+1;
    end;
end

for i=1:length(rBAC3)-1
    if rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1
        rBAC3(i,3)=11;
        TrMVIX(1,1)=TrMVIX(1,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0
        rBAC3(i,3)=12;
        TrMVIX(1,2)=TrMVIX(1,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1
        rBAC3(i,3)=13;
        TrMVIX(1,3)=TrMVIX(1,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1
        rBAC3(i,3)=21;
        TrMVIX(2,1)=TrMVIX(2,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0
        rBAC3(i,3)=22;
        TrMVIX(2,2)=TrMVIX(2,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1
        rBAC3(i,3)=23;
        TrMVIX(2,3)=TrMVIX(2,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1
        rBAC3(i,3)=31;
        TrMVIX(3,1)=TrMVIX(3,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0
        rBAC3(i,3)=32;
        TrMVIX(3,2)=TrMVIX(3,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==01
        rBAC3(i,3)=33;
    TrMVIX(3,3)=TrMVIX(3,3)+1;
    end;
end

TrMBACp=zeros(3);
TrMSPCp=zeros(3);
TrMVIXp=zeros(3);


TrMBACp(1,:)=TrMBAC(1,:)./sum(TrMBAC(1,:));
TrMBACp(2,:)=TrMBAC(2,:)./sum(TrMBAC(2,:));
TrMBACp(3,:)=TrMBAC(3,:)./sum(TrMBAC(3,:));

TrMSPp(1,:)=TrMSP(1,:)./sum(TrMSP(1,:));
TrMSPp(2,:)=TrMSP(2,:)./sum(TrMSP(2,:));
TrMSPp(3,:)=TrMSP(3,:)./sum(TrMSP(3,:));

TrMVIXp(1,:)=TrMVIX(1,:)./sum(TrMVIX(1,:));
TrMVIXp(2,:)=TrMVIX(2,:)./sum(TrMVIX(2,:));
TrMVIXp(3,:)=TrMVIX(3,:)./sum(TrMVIX(3,:));


SigSigSP=length(rBAC2);
T_BAC2=0;
T_BAC1=0;
T_BAC3=0;

for i=1:m
    for k=1:m
      T_BAC1=T_BAC1+2*TrMBAC(i,k)*log((TrMBAC(i,k)*sum(sum(TrMBAC)))/(sum(TrMBAC(i,:))*sum(TrMBAC(:,k))));
    end
end

for i=1:m
    for k=1:m
      T_BAC2=T_BAC2+2*TrMSP(i,k)*log((TrMSP(i,k)*SigSigSP)/(sum(TrMSP(i,:))*sum(TrMSP(:,k))));
    end
end

for i=1:m
    for k=1:m
        if TrMVIX(i,k)==0
            continue
        else
      T_BAC3=T_BAC3+2*TrMVIX(i,k)*log((TrMVIX(i,k)*sum(sum(TrMVIX)))/(sum(TrMVIX(i,:))*sum(TrMVIX(:,k))));
        end
    end
end

chi2cdf(T_BAC1,4)
chi2cdf(T_BAC2,4)
chi2cdf(T_BAC3,4)

Tr2MBAC=zeros(m,m,m);
Tr2MSP=zeros(m,m,m);
Tr2MVIX=zeros(m,m,m);

T2_BAC1=0;
T2_BAC2=0;
T2_BAC3=0;

for i=1:length(rBAC1)-2
    if rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=111;
        Tr2MBAC(1,1,1)=Tr2MBAC(1,1,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0
        rBAC1(i,4)=112;
        Tr2MBAC(1,1,2)=Tr2MBAC(1,1,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1
        rBAC1(i,4)=113;
        Tr2MBAC(1,1,3)=Tr2MBAC(1,1,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=121;
        Tr2MBAC(1,2,1)=Tr2MBAC(1,2,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0
        rBAC1(i,4)=122;
        Tr2MBAC(1,2,2)=Tr2MBAC(1,2,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1
        rBAC1(i,4)=123;
        Tr2MBAC(1,2,3)=Tr2MBAC(1,2,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=131;
        Tr2MBAC(1,3,1)=Tr2MBAC(1,3,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0
        rBAC1(i,4)=132;
        Tr2MBAC(1,3,2)=Tr2MBAC(1,3,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1
        rBAC1(i,4)=133;
    Tr2MBAC(1,3,3)=Tr2MBAC(1,3,3)+1;
    
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=211;
        Tr2MBAC(2,1,1)=Tr2MBAC(2,1,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0
        rBAC1(i,4)=212;
        Tr2MBAC(2,1,2)=Tr2MBAC(2,1,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1
        rBAC1(i,4)=213;
        Tr2MBAC(2,1,3)=Tr2MBAC(2,1,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=221;
        Tr2MBAC(2,2,1)=Tr2MBAC(2,2,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0
        rBAC1(i,4)=222;
        Tr2MBAC(2,2,2)=Tr2MBAC(2,2,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1
        rBAC1(i,4)=223;
        Tr2MBAC(2,2,3)=Tr2MBAC(2,2,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=231;
        Tr2MBAC(2,3,1)=Tr2MBAC(2,3,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0
        rBAC1(i,4)=232;
        Tr2MBAC(2,3,2)=Tr2MBAC(2,3,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1
        rBAC1(i,4)=233;
    Tr2MBAC(2,3,3)=Tr2MBAC(2,3,3)+1;
    
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=311;
        Tr2MBAC(3,1,1)=Tr2MBAC(3,1,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0
        rBAC1(i,4)=312;
        Tr2MBAC(3,1,2)=Tr2MBAC(3,1,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1
        rBAC1(i,4)=313;
        Tr2MBAC(3,1,3)=Tr2MBAC(3,1,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=321;
        Tr2MBAC(3,2,1)=Tr2MBAC(3,2,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0
        rBAC1(i,4)=322;
        Tr2MBAC(3,2,2)=Tr2MBAC(3,2,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1
        rBAC1(i,4)=323;
        Tr2MBAC(3,2,3)=Tr2MBAC(3,2,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1
        rBAC1(i,4)=331;
        Tr2MBAC(3,3,1)=Tr2MBAC(3,3,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0
        rBAC1(i,4)=332;
        Tr2MBAC(3,3,2)=Tr2MBAC(3,3,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1
        rBAC1(i,4)=333;
    Tr2MBAC(3,3,3)=Tr2MBAC(3,3,3)+1;
    end;
end



for i=1:length(rBAC2)-2
    if rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=111;
        Tr2MSP(1,1,1)=Tr2MSP(1,1,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0
        rBAC2(i,4)=112;
        Tr2MSP(1,1,2)=Tr2MSP(1,1,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1
        rBAC2(i,4)=113;
        Tr2MSP(1,1,3)=Tr2MSP(1,1,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=121;
        Tr2MSP(1,2,1)=Tr2MSP(1,2,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0
        rBAC2(i,4)=122;
        Tr2MSP(1,2,2)=Tr2MSP(1,2,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1
        rBAC2(i,4)=123;
        Tr2MSP(1,2,3)=Tr2MSP(1,2,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=131;
        Tr2MSP(1,3,1)=Tr2MSP(1,3,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0
        rBAC2(i,4)=132;
        Tr2MSP(1,3,2)=Tr2MSP(1,3,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1
        rBAC2(i,4)=133;
    Tr2MSP(1,3,3)=Tr2MSP(1,3,3)+1;
    
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=211;
        Tr2MSP(2,1,1)=Tr2MSP(2,1,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0
        rBAC2(i,4)=212;
        Tr2MSP(2,1,2)=Tr2MSP(2,1,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1
        rBAC2(i,4)=213;
        Tr2MSP(2,1,3)=Tr2MSP(2,1,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=221;
        Tr2MSP(2,2,1)=Tr2MSP(2,2,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0
        rBAC2(i,4)=222;
        Tr2MSP(2,2,2)=Tr2MSP(2,2,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1
        rBAC2(i,4)=223;
        Tr2MSP(2,2,3)=Tr2MSP(2,2,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=231;
        Tr2MSP(2,3,1)=Tr2MSP(2,3,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0
        rBAC2(i,4)=232;
        Tr2MSP(2,3,2)=Tr2MSP(2,3,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1
        rBAC2(i,4)=233;
    Tr2MSP(2,3,3)=Tr2MSP(2,3,3)+1;
    
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=311;
        Tr2MSP(3,1,1)=Tr2MSP(3,1,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0
        rBAC2(i,4)=312;
        Tr2MSP(3,1,2)=Tr2MSP(3,1,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1
        rBAC2(i,4)=313;
        Tr2MSP(3,1,3)=Tr2MSP(3,1,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=321;
        Tr2MSP(3,2,1)=Tr2MSP(3,2,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0
        rBAC2(i,4)=322;
        Tr2MSP(3,2,2)=Tr2MSP(3,2,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1
        rBAC2(i,4)=323;
        Tr2MSP(3,2,3)=Tr2MSP(3,2,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1
        rBAC2(i,4)=331;
        Tr2MSP(3,3,1)=Tr2MSP(3,3,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0
        rBAC2(i,4)=332;
        Tr2MSP(3,3,2)=Tr2MSP(3,3,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1
        rBAC2(i,4)=333;
    Tr2MSP(3,3,3)=Tr2MSP(3,3,3)+1;
    end;
end




for i=1:length(rBAC3)-2
    if rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=111;
        Tr2MVIX(1,1,1)=Tr2MVIX(1,1,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0
        rBAC3(i,4)=112;
        Tr2MVIX(1,1,2)=Tr2MVIX(1,1,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1
        rBAC3(i,4)=113;
        Tr2MVIX(1,1,3)=Tr2MVIX(1,1,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=121;
        Tr2MVIX(1,2,1)=Tr2MVIX(1,2,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0
        rBAC3(i,4)=122;
        Tr2MVIX(1,2,2)=Tr2MVIX(1,2,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1
        rBAC3(i,4)=123;
        Tr2MVIX(1,2,3)=Tr2MVIX(1,2,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=131;
        Tr2MVIX(1,3,1)=Tr2MVIX(1,3,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0
        rBAC3(i,4)=132;
        Tr2MVIX(1,3,2)=Tr2MVIX(1,3,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1
        rBAC3(i,4)=133;
    Tr2MVIX(1,3,3)=Tr2MVIX(1,3,3)+1;
    
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=211;
        Tr2MVIX(2,1,1)=Tr2MVIX(2,1,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0
        rBAC3(i,4)=212;
        Tr2MVIX(2,1,2)=Tr2MVIX(2,1,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1
        rBAC3(i,4)=213;
        Tr2MVIX(2,1,3)=Tr2MVIX(2,1,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=221;
        Tr2MVIX(2,2,1)=Tr2MVIX(2,2,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0
        rBAC3(i,4)=222;
        Tr2MVIX(2,2,2)=Tr2MVIX(2,2,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1
        rBAC3(i,4)=223;
        Tr2MVIX(2,2,3)=Tr2MVIX(2,2,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=231;
        Tr2MVIX(2,3,1)=Tr2MVIX(2,3,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0
        rBAC3(i,4)=232;
        Tr2MVIX(2,3,2)=Tr2MVIX(2,3,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1
        rBAC3(i,4)=233;
    Tr2MVIX(2,3,3)=Tr2MVIX(2,3,3)+1;
    
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=311;
        Tr2MVIX(3,1,1)=Tr2MVIX(3,1,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0
        rBAC3(i,4)=312;
        Tr2MVIX(3,1,2)=Tr2MVIX(3,1,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1
        rBAC3(i,4)=313;
        Tr2MVIX(3,1,3)=Tr2MVIX(3,1,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=321;
        Tr2MVIX(3,2,1)=Tr2MVIX(3,2,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0
        rBAC3(i,4)=322;
        Tr2MVIX(3,2,2)=Tr2MVIX(3,2,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1
        rBAC3(i,4)=323;
        Tr2MVIX(3,2,3)=Tr2MVIX(3,2,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1
        rBAC3(i,4)=331;
        Tr2MVIX(3,3,1)=Tr2MVIX(3,3,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0
        rBAC3(i,4)=332;
        Tr2MVIX(3,3,2)=Tr2MVIX(3,3,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1
        rBAC3(i,4)=333;
    Tr2MVIX(3,3,3)=Tr2MVIX(3,3,3)+1;
    end;
end


Tr2MBACp(1,1,:)=Tr2MBAC(1,1,:)./sum(Tr2MBAC(1,1,:));
Tr2MBACp(1,2,:)=Tr2MBAC(1,2,:)./sum(Tr2MBAC(1,2,:));
Tr2MBACp(1,3,:)=Tr2MBAC(1,3,:)./sum(Tr2MBAC(1,3,:));

Tr2MBACp(2,1,:)=Tr2MBAC(2,1,:)./sum(Tr2MBAC(2,1,:));
Tr2MBACp(2,2,:)=Tr2MBAC(2,2,:)./sum(Tr2MBAC(2,2,:));
Tr2MBACp(2,3,:)=Tr2MBAC(2,3,:)./sum(Tr2MBAC(2,3,:));

Tr2MBACp(3,1,:)=Tr2MBAC(3,1,:)./sum(Tr2MBAC(3,1,:));
Tr2MBACp(3,2,:)=Tr2MBAC(3,2,:)./sum(Tr2MBAC(3,2,:));
Tr2MBACp(3,3,:)=Tr2MBAC(3,3,:)./sum(Tr2MBAC(3,3,:));



Tr2MSPp(1,1,:)=Tr2MSP(1,1,:)./sum(Tr2MSP(1,1,:));
Tr2MSPp(1,2,:)=Tr2MSP(1,2,:)./sum(Tr2MSP(1,2,:));
Tr2MSPp(1,3,:)=Tr2MSP(1,3,:)./sum(Tr2MSP(1,3,:));

Tr2MSPp(2,1,:)=Tr2MSP(2,1,:)./sum(Tr2MSP(2,1,:));
Tr2MSPp(2,2,:)=Tr2MSP(2,2,:)./sum(Tr2MSP(2,2,:));
Tr2MSPp(2,3,:)=Tr2MSP(2,3,:)./sum(Tr2MSP(2,3,:));

Tr2MSPp(3,1,:)=Tr2MSP(3,1,:)./sum(Tr2MSP(3,1,:));
Tr2MSPp(3,2,:)=Tr2MSP(3,2,:)./sum(Tr2MSP(3,2,:));
Tr2MSPp(3,3,:)=Tr2MSP(3,3,:)./sum(Tr2MSP(3,3,:));


Tr2MVIXp(1,1,:)=Tr2MVIX(1,1,:)./sum(Tr2MVIX(1,1,:));
Tr2MVIXp(1,2,:)=Tr2MVIX(1,2,:)./sum(Tr2MVIX(1,2,:));
Tr2MVIXp(1,3,:)=Tr2MVIX(1,3,:)./sum(Tr2MVIX(1,3,:));

Tr2MVIXp(2,1,:)=Tr2MVIX(2,1,:)./sum(Tr2MVIX(2,1,:));
Tr2MVIXp(2,2,:)=Tr2MVIX(2,2,:)./sum(Tr2MVIX(2,2,:));
Tr2MVIXp(2,3,:)=Tr2MVIX(2,3,:)./sum(Tr2MVIX(2,3,:));

Tr2MVIXp(3,1,:)=Tr2MVIX(3,1,:)./sum(Tr2MVIX(3,1,:));
Tr2MVIXp(3,2,:)=Tr2MVIX(3,2,:)./sum(Tr2MVIX(3,2,:));
Tr2MVIXp(3,3,:)=Tr2MVIX(3,3,:)./sum(Tr2MVIX(3,3,:));


T2_BAC1=0;
T2_BAC2=0;
T2_BAC3=0;

for i=1:m
    for j=1:m
        for k=1:m
      T2_BAC1=T2_BAC1-0.1*Tr2MBAC(i,j,k)*log((Tr2MBAC(i,j,k)*sum(sum(sum(Tr2MBAC)))/(sum(sum(Tr2MBAC(i,:,:)))*sum(sum(Tr2MBAC(:,:,k))))));
        end
    end
end

for i=1:m
    for j=1:m
        for k=1:m
      T2_BAC2=T2_BAC2-0.1*Tr2MSP(i,j,k)*log((Tr2MSP(i,j,k)*sum(sum(sum(Tr2MSP)))/(sum(sum(Tr2MSP(i,:,:)))*sum(sum(Tr2MSP(:,:,k))))));
        end
    end
end

for i=1:m
    for j=1:m
    for k=1:m
        if Tr2MVIX(i,j,k)==0
            continue
        else
      T2_BAC3=T2_BAC3-0.1*Tr2MVIX(i,j,k)*log(Tr2MVIX(i,j,k)*sum(sum(sum(Tr2MVIX)))/(sum(sum(Tr2MVIX(i,:,:)))*sum(sum(Tr2MVIX(:,:,k)))));
        end
    end
    end
end

chi2cdf(T2_BAC1,12)
chi2cdf(T2_BAC2,12)
chi2cdf(T2_BAC3,12)

Tr3MBAC=zeros(m,m,m,m);
Tr3MSP=zeros(m,m,m,m);
Tr3MVIX=zeros(m,m,m,m);

T3_BAC1=0;
T3_BAC2=0;
T3_BAC3=0;

for i=1:length(rBAC1)-3
    if  rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1111;
        Tr3MBAC(1,1,1,1)=Tr3MBAC(1,1,1,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1112;
        Tr3MBAC(1,1,1,2)=Tr3MBAC(1,1,1,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1113;
        Tr3MBAC(1,1,1,3)=Tr3MBAC(1,1,1,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1121;
        Tr3MBAC(1,1,2,1)=Tr3MBAC(1,1,2,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1122;
        Tr3MBAC(1,1,2,2)=Tr3MBAC(1,1,2,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1123;
        Tr3MBAC(1,1,2,3)=Tr3MBAC(1,1,2,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1131;
        Tr3MBAC(1,1,3,1)=Tr3MBAC(1,1,3,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1132;
        Tr3MBAC(1,1,3,2)=Tr3MBAC(1,1,3,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1133;
    Tr3MBAC(1,1,3,3)=Tr3MBAC(1,1,3,3)+1;
    


    elseif  rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1211;
        Tr3MBAC(1,2,1,1)=Tr3MBAC(1,2,1,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1212;
        Tr3MBAC(1,2,1,2)=Tr3MBAC(1,2,1,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1213;
        Tr3MBAC(1,2,1,3)=Tr3MBAC(1,2,1,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1221;
        Tr3MBAC(1,2,2,1)=Tr3MBAC(1,2,2,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1222;
        Tr3MBAC(1,2,2,2)=Tr3MBAC(1,2,2,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1223;
        Tr3MBAC(1,2,2,3)=Tr3MBAC(1,2,2,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1231;
        Tr3MBAC(1,2,3,1)=Tr3MBAC(1,2,3,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1232;
        Tr3MBAC(1,2,3,2)=Tr3MBAC(1,2,3,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1233;
    Tr3MBAC(1,2,3,3)=Tr3MBAC(1,2,3,3)+1;
    


    elseif  rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1311;
        Tr3MBAC(1,3,1,1)=Tr3MBAC(1,3,1,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1312;
        Tr3MBAC(1,3,1,2)=Tr3MBAC(1,3,1,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1313;
        Tr3MBAC(1,3,1,3)=Tr3MBAC(1,3,1,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1321;
        Tr3MBAC(1,3,2,1)=Tr3MBAC(1,3,2,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1322;
        Tr3MBAC(1,3,2,2)=Tr3MBAC(1,3,2,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1323;
        Tr3MBAC(1,3,2,3)=Tr3MBAC(1,3,2,3)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=1331;
        Tr3MBAC(1,3,3,1)=Tr3MBAC(1,3,3,1)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=1332;
        Tr3MBAC(1,3,3,2)=Tr3MBAC(1,3,3,2)+1;
    elseif rBAC1(i,2)==-1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=1333;
    Tr3MBAC(1,3,3,3)=Tr3MBAC(1,3,3,3)+1;



    elseif  rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2111;
        Tr3MBAC(2,1,1,1)=Tr3MBAC(2,1,1,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2112;
        Tr3MBAC(2,1,1,2)=Tr3MBAC(2,1,1,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2113;
        Tr3MBAC(2,1,1,3)=Tr3MBAC(2,1,1,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2121;
        Tr3MBAC(2,1,2,1)=Tr3MBAC(2,1,2,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2122;
        Tr3MBAC(2,1,2,2)=Tr3MBAC(2,1,2,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2123;
        Tr3MBAC(2,1,2,3)=Tr3MBAC(2,1,2,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2131;
        Tr3MBAC(2,1,3,1)=Tr3MBAC(2,1,3,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2132;
        Tr3MBAC(2,1,3,2)=Tr3MBAC(2,1,3,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2133;
    Tr3MBAC(2,1,3,3)=Tr3MBAC(2,1,3,3)+1;
    


    elseif  rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2211;
        Tr3MBAC(2,2,1,1)=Tr3MBAC(2,2,1,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2212;
        Tr3MBAC(2,2,1,2)=Tr3MBAC(2,2,1,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2213;
        Tr3MBAC(2,2,1,3)=Tr3MBAC(2,2,1,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2221;
        Tr3MBAC(2,2,2,1)=Tr3MBAC(2,2,2,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2222;
        Tr3MBAC(2,2,2,2)=Tr3MBAC(2,2,2,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2223;
        Tr3MBAC(2,2,2,3)=Tr3MBAC(2,2,2,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2231;
        Tr3MBAC(2,2,3,1)=Tr3MBAC(2,2,3,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2232;
        Tr3MBAC(2,2,3,2)=Tr3MBAC(2,2,3,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2233;
    Tr3MBAC(2,2,3,3)=Tr3MBAC(2,2,3,3)+1;
    


    elseif  rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2311;
        Tr3MBAC(2,3,1,1)=Tr3MBAC(2,3,1,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2312;
        Tr3MBAC(2,3,1,2)=Tr3MBAC(2,3,1,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2313;
        Tr3MBAC(2,3,1,3)=Tr3MBAC(2,3,1,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2321;
        Tr3MBAC(2,3,2,1)=Tr3MBAC(2,3,2,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2322;
        Tr3MBAC(2,3,2,2)=Tr3MBAC(2,3,2,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2323;
        Tr3MBAC(2,3,2,3)=Tr3MBAC(2,3,2,3)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=2331;
        Tr3MBAC(2,3,3,1)=Tr3MBAC(2,3,3,1)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=2332;
        Tr3MBAC(2,3,3,2)=Tr3MBAC(2,3,3,2)+1;
    elseif rBAC1(i,2)==0&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=2333;
    Tr3MBAC(2,3,3,3)=Tr3MBAC(2,3,3,3)+1;

    elseif  rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3111;
        Tr3MBAC(3,1,1,1)=Tr3MBAC(3,1,1,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3112;
        Tr3MBAC(3,1,1,2)=Tr3MBAC(3,1,1,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3113;
        Tr3MBAC(3,1,1,3)=Tr3MBAC(3,1,1,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3121;
        Tr3MBAC(3,1,2,1)=Tr3MBAC(3,1,2,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3122;
        Tr3MBAC(3,1,2,2)=Tr3MBAC(3,1,2,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3123;
        Tr3MBAC(3,1,2,3)=Tr3MBAC(3,1,2,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3131;
        Tr3MBAC(3,1,3,1)=Tr3MBAC(3,1,3,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3132;
        Tr3MBAC(3,1,3,2)=Tr3MBAC(3,1,3,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==-1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3133;
    Tr3MBAC(3,1,3,3)=Tr3MBAC(3,1,3,3)+1;
    


    elseif  rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3211;
        Tr3MBAC(3,2,1,1)=Tr3MBAC(3,2,1,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3212;
        Tr3MBAC(3,2,1,2)=Tr3MBAC(3,2,1,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3213;
        Tr3MBAC(3,2,1,3)=Tr3MBAC(3,2,1,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3221;
        Tr3MBAC(3,2,2,1)=Tr3MBAC(3,2,2,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3222;
        Tr3MBAC(3,2,2,2)=Tr3MBAC(3,2,2,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3223;
        Tr3MBAC(3,2,2,3)=Tr3MBAC(3,2,2,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3231;
        Tr3MBAC(3,2,3,1)=Tr3MBAC(3,2,3,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3232;
        Tr3MBAC(3,2,3,2)=Tr3MBAC(3,2,3,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==0&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3233;
    Tr3MBAC(3,2,3,3)=Tr3MBAC(3,2,3,3)+1;
    


    elseif  rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3311;
        Tr3MBAC(3,3,1,1)=Tr3MBAC(3,3,1,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3312;
        Tr3MBAC(3,3,1,2)=Tr3MBAC(3,3,1,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==-1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3313;
        Tr3MBAC(3,3,1,3)=Tr3MBAC(3,3,1,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3321;
        Tr3MBAC(3,3,2,1)=Tr3MBAC(3,3,2,1)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3322;
        Tr3MBAC(3,3,2,2)=Tr3MBAC(3,3,2,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==0&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3323;
        Tr3MBAC(3,3,2,3)=Tr3MBAC(3,3,2,3)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==-1
        rBAC1(i,4)=3331;
        Tr3MBAC(3,3,3,2)=Tr3MBAC(3,3,3,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==0
        rBAC1(i,4)=3332;
        Tr3MBAC(3,3,3,2)=Tr3MBAC(3,3,3,2)+1;
    elseif rBAC1(i,2)==1&&rBAC1(i+1,2)==1&&rBAC1(i+2,2)==1&&rBAC1(i+3,2)==1
        rBAC1(i,4)=3333;
    Tr3MBAC(3,3,3,3)=Tr3MBAC(3,3,3,3)+1;

    end;
end


for i=1:length(rBAC2)-3
    if  rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1111;
        Tr3MSP(1,1,1,1)=Tr3MSP(1,1,1,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1112;
        Tr3MSP(1,1,1,2)=Tr3MSP(1,1,1,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1113;
        Tr3MSP(1,1,1,3)=Tr3MSP(1,1,1,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1121;
        Tr3MSP(1,1,2,1)=Tr3MSP(1,1,2,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1122;
        Tr3MSP(1,1,2,2)=Tr3MSP(1,1,2,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1123;
        Tr3MSP(1,1,2,3)=Tr3MSP(1,1,2,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1131;
        Tr3MSP(1,1,3,1)=Tr3MSP(1,1,3,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1132;
        Tr3MSP(1,1,3,2)=Tr3MSP(1,1,3,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1133;
    Tr3MSP(1,1,3,3)=Tr3MSP(1,1,3,3)+1;
    


    elseif  rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1211;
        Tr3MSP(1,2,1,1)=Tr3MSP(1,2,1,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1212;
        Tr3MSP(1,2,1,2)=Tr3MSP(1,2,1,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1213;
        Tr3MSP(1,2,1,3)=Tr3MSP(1,2,1,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1221;
        Tr3MSP(1,2,2,1)=Tr3MSP(1,2,2,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1222;
        Tr3MSP(1,2,2,2)=Tr3MSP(1,2,2,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1223;
        Tr3MSP(1,2,2,3)=Tr3MSP(1,2,2,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1231;
        Tr3MSP(1,2,3,1)=Tr3MSP(1,2,3,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1232;
        Tr3MSP(1,2,3,2)=Tr3MSP(1,2,3,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1233;
    Tr3MSP(1,2,3,3)=Tr3MSP(1,2,3,3)+1;
    


    elseif  rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1311;
        Tr3MSP(1,3,1,1)=Tr3MSP(1,3,1,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1312;
        Tr3MSP(1,3,1,2)=Tr3MSP(1,3,1,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1313;
        Tr3MSP(1,3,1,3)=Tr3MSP(1,3,1,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1321;
        Tr3MSP(1,3,2,1)=Tr3MSP(1,3,2,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1322;
        Tr3MSP(1,3,2,2)=Tr3MSP(1,3,2,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1323;
        Tr3MSP(1,3,2,3)=Tr3MSP(1,3,2,3)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=1331;
        Tr3MSP(1,3,3,1)=Tr3MSP(1,3,3,1)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=1332;
        Tr3MSP(1,3,3,2)=Tr3MSP(1,3,3,2)+1;
    elseif rBAC2(i,2)==-1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=1333;
    Tr3MSP(1,3,3,3)=Tr3MSP(1,3,3,3)+1;



    elseif  rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2111;
        Tr3MSP(2,1,1,1)=Tr3MSP(2,1,1,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2112;
        Tr3MSP(2,1,1,2)=Tr3MSP(2,1,1,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2113;
        Tr3MSP(2,1,1,3)=Tr3MSP(2,1,1,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2121;
        Tr3MSP(2,1,2,1)=Tr3MSP(2,1,2,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2122;
        Tr3MSP(2,1,2,2)=Tr3MSP(2,1,2,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2123;
        Tr3MSP(2,1,2,3)=Tr3MSP(2,1,2,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2131;
        Tr3MSP(2,1,3,1)=Tr3MSP(2,1,3,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2132;
        Tr3MSP(2,1,3,2)=Tr3MSP(2,1,3,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2133;
    Tr3MSP(2,1,3,3)=Tr3MSP(2,1,3,3)+1;
    


    elseif  rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2211;
        Tr3MSP(2,2,1,1)=Tr3MSP(2,2,1,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2212;
        Tr3MSP(2,2,1,2)=Tr3MSP(2,2,1,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2213;
        Tr3MSP(2,2,1,3)=Tr3MSP(2,2,1,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2221;
        Tr3MSP(2,2,2,1)=Tr3MSP(2,2,2,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2222;
        Tr3MSP(2,2,2,2)=Tr3MSP(2,2,2,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2223;
        Tr3MSP(2,2,2,3)=Tr3MSP(2,2,2,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2231;
        Tr3MSP(2,2,3,1)=Tr3MSP(2,2,3,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2232;
        Tr3MSP(2,2,3,2)=Tr3MSP(2,2,3,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2233;
    Tr3MSP(2,2,3,3)=Tr3MSP(2,2,3,3)+1;
    


    elseif  rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2311;
        Tr3MSP(2,3,1,1)=Tr3MSP(2,3,1,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2312;
        Tr3MSP(2,3,1,2)=Tr3MSP(2,3,1,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2313;
        Tr3MSP(2,3,1,3)=Tr3MSP(2,3,1,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2321;
        Tr3MSP(2,3,2,1)=Tr3MSP(2,3,2,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2322;
        Tr3MSP(2,3,2,2)=Tr3MSP(2,3,2,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2323;
        Tr3MSP(2,3,2,3)=Tr3MSP(2,3,2,3)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=2331;
        Tr3MSP(2,3,3,1)=Tr3MSP(2,3,3,1)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=2332;
        Tr3MSP(2,3,3,2)=Tr3MSP(2,3,3,2)+1;
    elseif rBAC2(i,2)==0&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=2333;
    Tr3MSP(2,3,3,3)=Tr3MSP(2,3,3,3)+1;

    elseif  rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3111;
        Tr3MSP(3,1,1,1)=Tr3MSP(3,1,1,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3112;
        Tr3MSP(3,1,1,2)=Tr3MSP(3,1,1,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3113;
        Tr3MSP(3,1,1,3)=Tr3MSP(3,1,1,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3121;
        Tr3MSP(3,1,2,1)=Tr3MSP(3,1,2,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3122;
        Tr3MSP(3,1,2,2)=Tr3MSP(3,1,2,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3123;
        Tr3MSP(3,1,2,3)=Tr3MSP(3,1,2,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3131;
        Tr3MSP(3,1,3,1)=Tr3MSP(3,1,3,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3132;
        Tr3MSP(3,1,3,2)=Tr3MSP(3,1,3,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==-1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3133;
    Tr3MSP(3,1,3,3)=Tr3MSP(3,1,3,3)+1;
    


    elseif  rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3211;
        Tr3MSP(3,2,1,1)=Tr3MSP(3,2,1,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3212;
        Tr3MSP(3,2,1,2)=Tr3MSP(3,2,1,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3213;
        Tr3MSP(3,2,1,3)=Tr3MSP(3,2,1,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3221;
        Tr3MSP(3,2,2,1)=Tr3MSP(3,2,2,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3222;
        Tr3MSP(3,2,2,2)=Tr3MSP(3,2,2,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3223;
        Tr3MSP(3,2,2,3)=Tr3MSP(3,2,2,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3231;
        Tr3MSP(3,2,3,1)=Tr3MSP(3,2,3,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3232;
        Tr3MSP(3,2,3,2)=Tr3MSP(3,2,3,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==0&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3233;
    Tr3MSP(3,2,3,3)=Tr3MSP(3,2,3,3)+1;
    


    elseif  rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3311;
        Tr3MSP(3,3,1,1)=Tr3MSP(3,3,1,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3312;
        Tr3MSP(3,3,1,2)=Tr3MSP(3,3,1,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==-1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3313;
        Tr3MSP(3,3,1,3)=Tr3MSP(3,3,1,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3321;
        Tr3MSP(3,3,2,1)=Tr3MSP(3,3,2,1)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3322;
        Tr3MSP(3,3,2,2)=Tr3MSP(3,3,2,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==0&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3323;
        Tr3MSP(3,3,2,3)=Tr3MSP(3,3,2,3)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==-1
        rBAC2(i,4)=3331;
        Tr3MSP(3,3,3,2)=Tr3MSP(3,3,3,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==0
        rBAC2(i,4)=3332;
        Tr3MSP(3,3,3,2)=Tr3MSP(3,3,3,2)+1;
    elseif rBAC2(i,2)==1&&rBAC2(i+1,2)==1&&rBAC2(i+2,2)==1&&rBAC2(i+3,2)==1
        rBAC2(i,4)=3333;
    Tr3MSP(3,3,3,3)=Tr3MSP(3,3,3,3)+1;

    end;
end


for i=1:length(rBAC3)-3
    if  rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1111;
        Tr3MVIX(1,1,1,1)=Tr3MVIX(1,1,1,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1112;
        Tr3MVIX(1,1,1,2)=Tr3MVIX(1,1,1,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1113;
        Tr3MVIX(1,1,1,3)=Tr3MVIX(1,1,1,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1121;
        Tr3MVIX(1,1,2,1)=Tr3MVIX(1,1,2,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1122;
        Tr3MVIX(1,1,2,2)=Tr3MVIX(1,1,2,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1123;
        Tr3MVIX(1,1,2,3)=Tr3MVIX(1,1,2,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1131;
        Tr3MVIX(1,1,3,1)=Tr3MVIX(1,1,3,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1132;
        Tr3MVIX(1,1,3,2)=Tr3MVIX(1,1,3,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1133;
    Tr3MVIX(1,1,3,3)=Tr3MVIX(1,1,3,3)+1;
    


    elseif  rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1211;
        Tr3MVIX(1,2,1,1)=Tr3MVIX(1,2,1,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1212;
        Tr3MVIX(1,2,1,2)=Tr3MVIX(1,2,1,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1213;
        Tr3MVIX(1,2,1,3)=Tr3MVIX(1,2,1,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1221;
        Tr3MVIX(1,2,2,1)=Tr3MVIX(1,2,2,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1222;
        Tr3MVIX(1,2,2,2)=Tr3MVIX(1,2,2,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1223;
        Tr3MVIX(1,2,2,3)=Tr3MVIX(1,2,2,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1231;
        Tr3MVIX(1,2,3,1)=Tr3MVIX(1,2,3,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1232;
        Tr3MVIX(1,2,3,2)=Tr3MVIX(1,2,3,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1233;
    Tr3MVIX(1,2,3,3)=Tr3MVIX(1,2,3,3)+1;
    


    elseif  rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1311;
        Tr3MVIX(1,3,1,1)=Tr3MVIX(1,3,1,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1312;
        Tr3MVIX(1,3,1,2)=Tr3MVIX(1,3,1,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1313;
        Tr3MVIX(1,3,1,3)=Tr3MVIX(1,3,1,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1321;
        Tr3MVIX(1,3,2,1)=Tr3MVIX(1,3,2,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1322;
        Tr3MVIX(1,3,2,2)=Tr3MVIX(1,3,2,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1323;
        Tr3MVIX(1,3,2,3)=Tr3MVIX(1,3,2,3)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=1331;
        Tr3MVIX(1,3,3,1)=Tr3MVIX(1,3,3,1)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=1332;
        Tr3MVIX(1,3,3,2)=Tr3MVIX(1,3,3,2)+1;
    elseif rBAC3(i,2)==-1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=1333;
    Tr3MVIX(1,3,3,3)=Tr3MVIX(1,3,3,3)+1;



    elseif  rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2111;
        Tr3MVIX(2,1,1,1)=Tr3MVIX(2,1,1,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2112;
        Tr3MVIX(2,1,1,2)=Tr3MVIX(2,1,1,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2113;
        Tr3MVIX(2,1,1,3)=Tr3MVIX(2,1,1,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2121;
        Tr3MVIX(2,1,2,1)=Tr3MVIX(2,1,2,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2122;
        Tr3MVIX(2,1,2,2)=Tr3MVIX(2,1,2,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2123;
        Tr3MVIX(2,1,2,3)=Tr3MVIX(2,1,2,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2131;
        Tr3MVIX(2,1,3,1)=Tr3MVIX(2,1,3,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2132;
        Tr3MVIX(2,1,3,2)=Tr3MVIX(2,1,3,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2133;
    Tr3MVIX(2,1,3,3)=Tr3MVIX(2,1,3,3)+1;
    


    elseif  rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2211;
        Tr3MVIX(2,2,1,1)=Tr3MVIX(2,2,1,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2212;
        Tr3MVIX(2,2,1,2)=Tr3MVIX(2,2,1,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2213;
        Tr3MVIX(2,2,1,3)=Tr3MVIX(2,2,1,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2221;
        Tr3MVIX(2,2,2,1)=Tr3MVIX(2,2,2,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2222;
        Tr3MVIX(2,2,2,2)=Tr3MVIX(2,2,2,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2223;
        Tr3MVIX(2,2,2,3)=Tr3MVIX(2,2,2,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2231;
        Tr3MVIX(2,2,3,1)=Tr3MVIX(2,2,3,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2232;
        Tr3MVIX(2,2,3,2)=Tr3MVIX(2,2,3,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2233;
    Tr3MVIX(2,2,3,3)=Tr3MVIX(2,2,3,3)+1;
    


    elseif  rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2311;
        Tr3MVIX(2,3,1,1)=Tr3MVIX(2,3,1,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2312;
        Tr3MVIX(2,3,1,2)=Tr3MVIX(2,3,1,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2313;
        Tr3MVIX(2,3,1,3)=Tr3MVIX(2,3,1,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2321;
        Tr3MVIX(2,3,2,1)=Tr3MVIX(2,3,2,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2322;
        Tr3MVIX(2,3,2,2)=Tr3MVIX(2,3,2,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2323;
        Tr3MVIX(2,3,2,3)=Tr3MVIX(2,3,2,3)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=2331;
        Tr3MVIX(2,3,3,1)=Tr3MVIX(2,3,3,1)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=2332;
        Tr3MVIX(2,3,3,2)=Tr3MVIX(2,3,3,2)+1;
    elseif rBAC3(i,2)==0&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=2333;
    Tr3MVIX(2,3,3,3)=Tr3MVIX(2,3,3,3)+1;

    elseif  rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3111;
        Tr3MVIX(3,1,1,1)=Tr3MVIX(3,1,1,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3112;
        Tr3MVIX(3,1,1,2)=Tr3MVIX(3,1,1,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3113;
        Tr3MVIX(3,1,1,3)=Tr3MVIX(3,1,1,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3121;
        Tr3MVIX(3,1,2,1)=Tr3MVIX(3,1,2,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3122;
        Tr3MVIX(3,1,2,2)=Tr3MVIX(3,1,2,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3123;
        Tr3MVIX(3,1,2,3)=Tr3MVIX(3,1,2,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3131;
        Tr3MVIX(3,1,3,1)=Tr3MVIX(3,1,3,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3132;
        Tr3MVIX(3,1,3,2)=Tr3MVIX(3,1,3,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==-1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3133;
    Tr3MVIX(3,1,3,3)=Tr3MVIX(3,1,3,3)+1;
    


    elseif  rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3211;
        Tr3MVIX(3,2,1,1)=Tr3MVIX(3,2,1,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3212;
        Tr3MVIX(3,2,1,2)=Tr3MVIX(3,2,1,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3213;
        Tr3MVIX(3,2,1,3)=Tr3MVIX(3,2,1,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3221;
        Tr3MVIX(3,2,2,1)=Tr3MVIX(3,2,2,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3222;
        Tr3MVIX(3,2,2,2)=Tr3MVIX(3,2,2,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3223;
        Tr3MVIX(3,2,2,3)=Tr3MVIX(3,2,2,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3231;
        Tr3MVIX(3,2,3,1)=Tr3MVIX(3,2,3,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3232;
        Tr3MVIX(3,2,3,2)=Tr3MVIX(3,2,3,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==0&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3233;
    Tr3MVIX(3,2,3,3)=Tr3MVIX(3,2,3,3)+1;
    


    elseif  rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3311;
        Tr3MVIX(3,3,1,1)=Tr3MVIX(3,3,1,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3312;
        Tr3MVIX(3,3,1,2)=Tr3MVIX(3,3,1,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==-1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3313;
        Tr3MVIX(3,3,1,3)=Tr3MVIX(3,3,1,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3321;
        Tr3MVIX(3,3,2,1)=Tr3MVIX(3,3,2,1)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3322;
        Tr3MVIX(3,3,2,2)=Tr3MVIX(3,3,2,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==0&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3323;
        Tr3MVIX(3,3,2,3)=Tr3MVIX(3,3,2,3)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==-1
        rBAC3(i,4)=3331;
        Tr3MVIX(3,3,3,2)=Tr3MVIX(3,3,3,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==0
        rBAC3(i,4)=3332;
        Tr3MVIX(3,3,3,2)=Tr3MVIX(3,3,3,2)+1;
    elseif rBAC3(i,2)==1&&rBAC3(i+1,2)==1&&rBAC3(i+2,2)==1&&rBAC3(i+3,2)==1
        rBAC3(i,4)=3333;
    Tr3MVIX(3,3,3,3)=Tr3MVIX(3,3,3,3)+1;

    end;
end


Tr3MBACp(1,1,1,:)=Tr3MBAC(1,1,1,:)./sum(Tr3MBAC(1,1,1,:));
Tr3MBACp(1,1,2,:)=Tr3MBAC(1,1,2,:)./sum(Tr3MBAC(1,1,2,:));
Tr3MBACp(1,1,3,:)=Tr3MBAC(1,1,3,:)./sum(Tr3MBAC(1,1,3,:));

Tr3MBACp(1,2,1,:)=Tr3MBAC(1,2,1,:)./sum(Tr3MBAC(1,2,1,:));
Tr3MBACp(1,2,2,:)=Tr3MBAC(1,2,2,:)./sum(Tr3MBAC(1,2,2,:));
Tr3MBACp(1,2,3,:)=Tr3MBAC(1,2,3,:)./sum(Tr3MBAC(1,2,3,:));

Tr3MBACp(1,3,1,:)=Tr3MBAC(1,3,1,:)./sum(Tr3MBAC(1,3,1,:));
Tr3MBACp(1,3,2,:)=Tr3MBAC(1,3,2,:)./sum(Tr3MBAC(1,3,2,:));
Tr3MBACp(1,3,3,:)=Tr3MBAC(1,3,3,:)./sum(Tr3MBAC(1,3,3,:));


Tr3MBACp(2,1,1,:)=Tr3MBAC(2,1,1,:)./sum(Tr3MBAC(2,1,1,:));
Tr3MBACp(2,1,2,:)=Tr3MBAC(2,1,2,:)./sum(Tr3MBAC(2,1,2,:));
Tr3MBACp(2,1,3,:)=Tr3MBAC(2,1,3,:)./sum(Tr3MBAC(2,1,3,:));

Tr3MBACp(2,2,1,:)=Tr3MBAC(2,2,1,:)./sum(Tr3MBAC(2,2,1,:));
Tr3MBACp(2,2,2,:)=Tr3MBAC(2,2,2,:)./sum(Tr3MBAC(2,2,2,:));
Tr3MBACp(2,2,3,:)=Tr3MBAC(2,2,3,:)./sum(Tr3MBAC(2,2,3,:));

Tr3MBACp(2,3,1,:)=Tr3MBAC(2,3,1,:)./sum(Tr3MBAC(2,3,1,:));
Tr3MBACp(2,3,2,:)=Tr3MBAC(2,3,2,:)./sum(Tr3MBAC(2,3,2,:));
Tr3MBACp(2,3,3,:)=Tr3MBAC(2,3,3,:)./sum(Tr3MBAC(2,3,3,:));

Tr3MBACp(3,1,1,:)=Tr3MBAC(3,1,1,:)./sum(Tr3MBAC(3,1,1,:));
Tr3MBACp(3,1,2,:)=Tr3MBAC(3,1,2,:)./sum(Tr3MBAC(3,1,2,:));
Tr3MBACp(3,1,3,:)=Tr3MBAC(3,1,3,:)./sum(Tr3MBAC(3,1,3,:));

Tr3MBACp(3,2,1,:)=Tr3MBAC(3,2,1,:)./sum(Tr3MBAC(3,2,1,:));
Tr3MBACp(3,2,2,:)=Tr3MBAC(3,2,2,:)./sum(Tr3MBAC(3,2,2,:));
Tr3MBACp(3,2,3,:)=Tr3MBAC(3,2,3,:)./sum(Tr3MBAC(3,2,3,:));

Tr3MBACp(3,3,1,:)=Tr3MBAC(3,3,1,:)./sum(Tr3MBAC(3,3,1,:));
Tr3MBACp(3,3,2,:)=Tr3MBAC(3,3,2,:)./sum(Tr3MBAC(3,3,2,:));
Tr3MBACp(3,3,3,:)=Tr3MBAC(3,3,3,:)./sum(Tr3MBAC(3,3,3,:));


Tr3MSPp(1,1,1,:)=Tr3MSP(1,1,1,:)./sum(Tr3MSP(1,1,1,:));
Tr3MSPp(1,1,2,:)=Tr3MSP(1,1,2,:)./sum(Tr3MSP(1,1,2,:));
Tr3MSPp(1,1,3,:)=Tr3MSP(1,1,3,:)./sum(Tr3MSP(1,1,3,:));

Tr3MSPp(1,2,1,:)=Tr3MSP(1,2,1,:)./sum(Tr3MSP(1,2,1,:));
Tr3MSPp(1,2,2,:)=Tr3MSP(1,2,2,:)./sum(Tr3MSP(1,2,2,:));
Tr3MSPp(1,2,3,:)=Tr3MSP(1,2,3,:)./sum(Tr3MSP(1,2,3,:));

Tr3MSPp(1,3,1,:)=Tr3MSP(1,3,1,:)./sum(Tr3MSP(1,3,1,:));
Tr3MSPp(1,3,2,:)=Tr3MSP(1,3,2,:)./sum(Tr3MSP(1,3,2,:));
Tr3MSPp(1,3,3,:)=Tr3MSP(1,3,3,:)./sum(Tr3MSP(1,3,3,:));


Tr3MSPp(2,1,1,:)=Tr3MSP(2,1,1,:)./sum(Tr3MSP(2,1,1,:));
Tr3MSPp(2,1,2,:)=Tr3MSP(2,1,2,:)./sum(Tr3MSP(2,1,2,:));
Tr3MSPp(2,1,3,:)=Tr3MSP(2,1,3,:)./sum(Tr3MSP(2,1,3,:));

Tr3MSPp(2,2,1,:)=Tr3MSP(2,2,1,:)./sum(Tr3MSP(2,2,1,:));
Tr3MSPp(2,2,2,:)=Tr3MSP(2,2,2,:)./sum(Tr3MSP(2,2,2,:));
Tr3MSPp(2,2,3,:)=Tr3MSP(2,2,3,:)./sum(Tr3MSP(2,2,3,:));

Tr3MSPp(2,3,1,:)=Tr3MSP(2,3,1,:)./sum(Tr3MSP(2,3,1,:));
Tr3MSPp(2,3,2,:)=Tr3MSP(2,3,2,:)./sum(Tr3MSP(2,3,2,:));
Tr3MSPp(2,3,3,:)=Tr3MSP(2,3,3,:)./sum(Tr3MSP(2,3,3,:));

Tr3MSPp(3,1,1,:)=Tr3MSP(3,1,1,:)./sum(Tr3MSP(3,1,1,:));
Tr3MSPp(3,1,2,:)=Tr3MSP(3,1,2,:)./sum(Tr3MSP(3,1,2,:));
Tr3MSPp(3,1,3,:)=Tr3MSP(3,1,3,:)./sum(Tr3MSP(3,1,3,:));

Tr3MSPp(3,2,1,:)=Tr3MSP(3,2,1,:)./sum(Tr3MSP(3,2,1,:));
Tr3MSPp(3,2,2,:)=Tr3MSP(3,2,2,:)./sum(Tr3MSP(3,2,2,:));
Tr3MSPp(3,2,3,:)=Tr3MSP(3,2,3,:)./sum(Tr3MSP(3,2,3,:));

Tr3MSPp(3,3,1,:)=Tr3MSP(3,3,1,:)./sum(Tr3MSP(3,3,1,:));
Tr3MSPp(3,3,2,:)=Tr3MSP(3,3,2,:)./sum(Tr3MSP(3,3,2,:));
Tr3MSPp(3,3,3,:)=Tr3MSP(3,3,3,:)./sum(Tr3MSP(3,3,3,:));

Tr3MVIXp(1,1,1,:)=Tr3MVIX(1,1,1,:)./sum(Tr3MVIX(1,1,1,:));
Tr3MVIXp(1,1,2,:)=Tr3MVIX(1,1,2,:)./sum(Tr3MVIX(1,1,2,:));
Tr3MVIXp(1,1,3,:)=Tr3MVIX(1,1,3,:)./sum(Tr3MVIX(1,1,3,:));

Tr3MVIXp(1,2,1,:)=Tr3MVIX(1,2,1,:)./sum(Tr3MVIX(1,2,1,:));
Tr3MVIXp(1,2,2,:)=Tr3MVIX(1,2,2,:)./sum(Tr3MVIX(1,2,2,:));
Tr3MVIXp(1,2,3,:)=Tr3MVIX(1,2,3,:)./sum(Tr3MVIX(1,2,3,:));

Tr3MVIXp(1,3,1,:)=Tr3MVIX(1,3,1,:)./sum(Tr3MVIX(1,3,1,:));
Tr3MVIXp(1,3,2,:)=Tr3MVIX(1,3,2,:)./sum(Tr3MVIX(1,3,2,:));
Tr3MVIXp(1,3,3,:)=Tr3MVIX(1,3,3,:)./sum(Tr3MVIX(1,3,3,:));


Tr3MVIXp(2,1,1,:)=Tr3MVIX(2,1,1,:)./sum(Tr3MVIX(2,1,1,:));
Tr3MVIXp(2,1,2,:)=Tr3MVIX(2,1,2,:)./sum(Tr3MVIX(2,1,2,:));
Tr3MVIXp(2,1,3,:)=Tr3MVIX(2,1,3,:)./sum(Tr3MVIX(2,1,3,:));

Tr3MVIXp(2,2,1,:)=Tr3MVIX(2,2,1,:)./sum(Tr3MVIX(2,2,1,:));
Tr3MVIXp(2,2,2,:)=Tr3MVIX(2,2,2,:)./sum(Tr3MVIX(2,2,2,:));
Tr3MVIXp(2,2,3,:)=Tr3MVIX(2,2,3,:)./sum(Tr3MVIX(2,2,3,:));

Tr3MVIXp(2,3,1,:)=Tr3MVIX(2,3,1,:)./sum(Tr3MVIX(2,3,1,:));
Tr3MVIXp(2,3,2,:)=Tr3MVIX(2,3,2,:)./sum(Tr3MVIX(2,3,2,:));
Tr3MVIXp(2,3,3,:)=Tr3MVIX(2,3,3,:)./sum(Tr3MVIX(2,3,3,:));

Tr3MVIXp(3,1,1,:)=Tr3MVIX(3,1,1,:)./sum(Tr3MVIX(3,1,1,:));
Tr3MVIXp(3,1,2,:)=Tr3MVIX(3,1,2,:)./sum(Tr3MVIX(3,1,2,:));
Tr3MVIXp(3,1,3,:)=Tr3MVIX(3,1,3,:)./sum(Tr3MVIX(3,1,3,:));

Tr3MVIXp(3,2,1,:)=Tr3MVIX(3,2,1,:)./sum(Tr3MVIX(3,2,1,:));
Tr3MVIXp(3,2,2,:)=Tr3MVIX(3,2,2,:)./sum(Tr3MVIX(3,2,2,:));
Tr3MVIXp(3,2,3,:)=Tr3MVIX(3,2,3,:)./sum(Tr3MVIX(3,2,3,:));

Tr3MVIXp(3,3,1,:)=Tr3MVIX(3,3,1,:)./sum(Tr3MVIX(3,3,1,:));
Tr3MVIXp(3,3,2,:)=Tr3MVIX(3,3,2,:)./sum(Tr3MVIX(3,3,2,:));
Tr3MVIXp(3,3,3,:)=Tr3MVIX(3,3,3,:)./sum(Tr3MVIX(3,3,3,:));


T3_BAC1=0;
T3_BAC2=0;
T3_BAC3=0;

for i=1:m
    for j=1:m
        for h=1:m
        for k=1:m
      T3_BAC1=T3_BAC1-0.1*min(Tr3MBAC(i,j,h,k)*log((Tr3MBAC(i,j,h,k)*sum(sum(sum(sum(Tr3MBAC))))/(sum(sum(sum(Tr3MBAC(i,:,:,:))))*sum(sum(sum(Tr3MBAC(:,:,:,k))))))),0);
        end
        end
    end
end

for i=1:m
    for j=1:m
        for h=1:m
        for k=1:m
      T3_BAC2=T3_BAC2-0.1*min(Tr3MSP(i,j,h,k)*log((Tr3MSP(i,j,h,k)*sum(sum(sum(sum(Tr3MSP))))/(sum(sum(sum(Tr3MSP(i,:,:,:))))*sum(sum(sum(Tr3MSP(:,:,:,k))))))),0);
        end
        end
    end
end

for i=1:m
    for j=1:m
        for h=1:m
    for k=1:m
        if Tr3MVIX(i,j,k)==0
            continue
        else
      T3_BAC3=T3_BAC3-0.1*min(Tr3MVIX(i,j,h,k)*log(Tr3MVIX(i,j,h,k)*sum(sum(sum(sum(Tr3MVIX))))/(sum(sum(sum(Tr3MVIX(i,:,:,:))))*sum(sum(sum(Tr3MVIX(:,:,:,k)))))),0);
        end
    end
        end
    end
end

chi2cdf(T3_BAC1,36)
chi2cdf(T3_BAC2,36)
chi2cdf(T3_BAC3,36)


