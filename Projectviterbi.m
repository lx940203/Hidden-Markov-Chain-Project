m=3;
Parament=cell(500,2);
Expectation=zeros(500,1);

mu=[-1,0,1];% The Mu and Sigma can give by EM algorithm
sigma=[1,1,1];


SP=xlsread('J:\×ÀÃæ\Project\^GSPC.csv');
SP=SP(:,end-1);
y=(SP(2:end)-SP(1:end-1))./SP(1:end-1);
y=SP;
n=length(y);
x=zeros(3,1);
x(1)=sum(y<=prctile(y,33));
x(3)=sum(y>prctile(y,66));
x(2)=length(y)-x(1)-x(3);

Ptem=x/sum(x); 


%rBAC3=rBAC3(1:end-1,end-1);
rSPs=sort(y);
SP_DtI=rSPs(floor(length(y)/3));
SP_ItU=rSPs(floor(length(y)/3*2));

m=3;
y=[y,zeros(length(y),1),zeros(length(y),1),zeros(length(y),1)];


TrMSP=zeros(m);


for i=1:length(y)
    if y(i,1)>SP_ItU
        y(i,2)=1;
    elseif y(i,1)<SP_DtI
        y(i,2)=-1;
    else y(i,2)=0;
    end;
end




for i=1:length(y)-1
    if y(i,2)==-1&&y(i+1,2)==-1
        y(i,3)=11;
        TrMSP(1,1)=TrMSP(1,1)+1;
    elseif y(i,2)==-1&&y(i+1,2)==0
        y(i,3)=12;
        TrMSP(1,2)=TrMSP(1,2)+1;
    elseif y(i,2)==-1&&y(i+1,2)==1
        y(i,3)=13;
        TrMSP(1,3)=TrMSP(1,3)+1;
    elseif y(i,2)==0&&y(i+1,2)==-1
        y(i,3)=21;
        TrMSP(2,1)=TrMSP(2,1)+1;
    elseif y(i,2)==0&&y(i+1,2)==0
        y(i,3)=22;
        TrMSP(2,2)=TrMSP(2,2)+1;
    elseif y(i,2)==0&&y(i+1,2)==1
        y(i,3)=23;
        TrMSP(2,3)=TrMSP(2,3)+1;
    elseif y(i,2)==1&&y(i+1,2)==-1
        y(i,3)=31;
        TrMSP(3,1)=TrMSP(3,1)+1;
    elseif y(i,2)==1&&y(i+1,2)==0
        y(i,3)=32;
        TrMSP(3,2)=TrMSP(3,2)+1;
    elseif y(i,2)==1&&y(i+1,2)==01
        y(i,3)=33;
    TrMSP(3,3)=TrMSP(3,3)+1;
    end;
end

TrMSPCp=zeros(3);


TrMSPp(1,:)=TrMSP(1,:)./sum(TrMSP(1,:));
TrMSPp(2,:)=TrMSP(2,:)./sum(TrMSP(2,:));
TrMSPp(3,:)=TrMSP(3,:)./sum(TrMSP(3,:));

delta=zeros(length(y),3);
fai=zeros(length(y),3);
faitem=zeros(3,3);

fai2=zeros(length(y),1);
faitem2=zeros(1,3);

i=1;
p1=1./((2*pi.*sigma(1).^2).^0.5).*exp(-(log(y(i))-mu(1)).^2./(2.*sigma(1).^2));
p2=1./((2*pi.*sigma(2).^2).^0.5).*exp(-(log(y(i))-mu(2)).^2./(2.*sigma(2).^2));
p3=1./((2*pi.*sigma(3).^2).^0.5).*exp(-(log(y(i))-mu(3)).^2./(2.*sigma(3).^2));
pt1=p1/(p1+p2+p3);
pt2=p2/(p1+p2+p3);
pt3=p3/(p1+p2+p3);

delta(1,1)=pt1*Ptem(1);
delta(1,2)=pt1*Ptem(2);
delta(1,3)=pt1*Ptem(3);



for i=2:length(y)
    p1=1./((2*pi.*sigma.^2).^0.5).*exp(-(log(y(i))-mu(1)).^2./(2.*sigma.^2));
    p2=1./((2*pi.*sigma.^2).^0.5).*exp(-(log(y(i))-mu(2)).^2./(2.*sigma.^2));
    p3=1./((2*pi.*sigma.^2).^0.5).*exp(-(log(y(i))-mu(3)).^2./(2.*sigma.^2));
    pt1=p1/(p1+p2+p3);
    pt2=p2/(p1+p2+p3);
    pt3=p3/(p1+p2+p3);
    delta(i,1)=max(delta(i-1,:).*TrMSPp(:,1)')*pt1;
    delta(i,2)=max(delta(i-1,:).*TrMSPp(:,2)')*pt2;
    delta(i,3)=max(delta(i-1,:).*TrMSPp(:,3)')*pt3;
    
    faitem(:,1)=(delta(i-1,:).*TrMSPp(:,1)');
    faitem(:,2)=(delta(i-1,:).*TrMSPp(:,2)');
    faitem(:,3)=(delta(i-1,:).*TrMSPp(:,3)');
    
    fai(i,1)=faitem(faitem(:,1)==max(faitem(:,1)),1);
    fai(i,2)=faitem(faitem(:,2)==max(faitem(:,2)),2);
    fai(i,3)=faitem(faitem(:,3)==max(faitem(:,3)),3);
end 
    
SelectDelta=max(delta(end,:));
delta(end,1)=delta(end,3)*0.9;
delta(end,2)=delta(end,3)*0.9;
fai2(end)=find(delta(end,:)==max(delta(end,:)));

for i=length(y)-1:-1:1
    faitem2(:)=(delta(i+1,:).*TrMSPp(:,fai2(i+1))');   
    fai2(i,1)=find(faitem2(:)==max(faitem2(:))); 
end