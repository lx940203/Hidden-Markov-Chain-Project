BAC=xlsread('E:\×ÀÃæ\Project\BAC.csv');
BAC=BAC(:,end-1);
rBAC=(BAC(2:end)-BAC(1:end-1))./BAC(1:end-1);
SP=xlsread('E:\×ÀÃæ\Project\^GSPC.csv');
SP=SP(:,end-1);
rSP=(SP(2:end)-SP(1:end-1))./SP(1:end-1);
VIX=xlsread('E:\×ÀÃæ\Project\^VIX.csv',1);
VIX=VIX(1:end-1,end-1);
rSPs=sort(rSP);
rBACs=sort(rBAC);
VIXs=sort(VIX);
SP_DtI=rSPs(floor(length(rSP)/3));
SP_ItU=rSPs(floor(length(rSP)/3*2));
BAC_DtI=rBACs(floor(length(rBAC)/3));
BAC_ItU=rBACs(floor(length(rBAC)/3*2));
VIX_DtI=VIXs(floor(length(VIX)/3));
VIX_ItU=VIXs(floor(length(VIX)/3*2));

m=3;
rBAC=[rBAC,zeros(length(rBAC),1),zeros(length(rBAC),1)];
rSP=[rSP,zeros(length(rSP),1),zeros(length(rBAC),1)];
VIX=[VIX,zeros(length(VIX),1),zeros(length(rBAC),1)];

TrMBAC=zeros(m);
TrMSP=zeros(m);
TrMVIX=zeros(m);

for i=1:length(rBAC)
    if rBAC(i,1)>BAC_ItU
        rBAC(i,2)=1;
    elseif rBAC(i,1)<BAC_DtI
        rBAC(i,2)=-1;
    else rBAC(i,2)=0;
    end;
end

for i=1:length(rSP)
    if rSP(i,1)>SP_ItU
        rSP(i,2)=1;
    elseif rSP(i,1)<SP_DtI
        rSP(i,2)=-1;
    else rSP(i,2)=0;
    end;
end

for i=1:length(VIX)
    if VIX(i,1)>VIX_ItU
        VIX(i,2)=1;
    elseif VIX(i,1)<VIX_DtI
        VIX(i,2)=-1;
    else VIX(i,2)=0;
    end;
end

for i=1:length(rBAC)-1
    if rBAC(i,2)==-1&&rBAC(i+1,2)==-1
        rBAC(i,3)=11;
        TrMBAC(1,1)=TrMBAC(1,1)+1;
    elseif rBAC(i,2)==-1&&rBAC(i+1,2)==0
        rBAC(i,3)=12;
        TrMBAC(1,2)=TrMBAC(1,2)+1;
    elseif rBAC(i,2)==-1&&rBAC(i+1,2)==1
        rBAC(i,3)=13;
        TrMBAC(1,3)=TrMBAC(1,3)+1;
    elseif rBAC(i,2)==0&&rBAC(i+1,2)==-1
        rBAC(i,3)=21;
        TrMBAC(2,1)=TrMBAC(2,1)+1;
    elseif rBAC(i,2)==0&&rBAC(i+1,2)==0
        rBAC(i,3)=22;
        TrMBAC(2,2)=TrMBAC(2,2)+1;
    elseif rBAC(i,2)==0&&rBAC(i+1,2)==1
        rBAC(i,3)=23;
        TrMBAC(2,3)=TrMBAC(2,3)+1;
    elseif rBAC(i,2)==1&&rBAC(i+1,2)==-1
        rBAC(i,3)=31;
        TrMBAC(3,1)=TrMBAC(3,1)+1;
    elseif rBAC(i,2)==1&&rBAC(i+1,2)==0
        rBAC(i,3)=32;
        TrMBAC(3,2)=TrMBAC(3,2)+1;
    elseif rBAC(i,2)==1&&rBAC(i+1,2)==01
        rBAC(i,3)=33;
    TrMBAC(3,3)=TrMBAC(3,3)+1;
    end;
end

for i=1:length(rSP)-1
    if rSP(i,2)==-1&&rSP(i+1,2)==-1
        rSP(i,3)=11;
        TrMSP(1,1)=TrMSP(1,1)+1;
    elseif rSP(i,2)==-1&&rSP(i+1,2)==0
        rSP(i,3)=12;
        TrMSP(1,2)=TrMSP(1,2)+1;
    elseif rSP(i,2)==-1&&rSP(i+1,2)==1
        rSP(i,3)=13;
        TrMSP(1,3)=TrMSP(1,3)+1;
    elseif rSP(i,2)==0&&rSP(i+1,2)==-1
        rSP(i,3)=21;
        TrMSP(2,1)=TrMSP(2,1)+1;
    elseif rSP(i,2)==0&&rSP(i+1,2)==0
        rSP(i,3)=22;
        TrMSP(2,2)=TrMSP(2,2)+1;
    elseif rSP(i,2)==0&&rSP(i+1,2)==1
        rSP(i,3)=23;
        TrMSP(2,3)=TrMSP(2,3)+1;
    elseif rSP(i,2)==1&&rSP(i+1,2)==-1
        rSP(i,3)=31;
        TrMSP(3,1)=TrMSP(3,1)+1;
    elseif rSP(i,2)==1&&rSP(i+1,2)==0
        rSP(i,3)=32;
        TrMSP(3,2)=TrMSP(3,2)+1;
    elseif rSP(i,2)==1&&rSP(i+1,2)==01
        rSP(i,3)=33;
    TrMSP(3,3)=TrMSP(3,3)+1;
    end;
end

for i=1:length(VIX)-1
    if VIX(i,2)==-1&&VIX(i+1,2)==-1
        VIX(i,3)=11;
        TrMVIX(1,1)=TrMVIX(1,1)+1;
    elseif VIX(i,2)==-1&&VIX(i+1,2)==0
        VIX(i,3)=12;
        TrMVIX(1,2)=TrMVIX(1,2)+1;
    elseif VIX(i,2)==-1&&VIX(i+1,2)==1
        VIX(i,3)=13;
        TrMVIX(1,3)=TrMVIX(1,3)+1;
    elseif VIX(i,2)==0&&VIX(i+1,2)==-1
        VIX(i,3)=21;
        TrMVIX(2,1)=TrMVIX(2,1)+1;
    elseif VIX(i,2)==0&&VIX(i+1,2)==0
        VIX(i,3)=22;
        TrMVIX(2,2)=TrMVIX(2,2)+1;
    elseif VIX(i,2)==0&&VIX(i+1,2)==1
        VIX(i,3)=23;
        TrMVIX(2,3)=TrMVIX(2,3)+1;
    elseif VIX(i,2)==1&&VIX(i+1,2)==-1
        VIX(i,3)=31;
        TrMVIX(3,1)=TrMVIX(3,1)+1;
    elseif VIX(i,2)==1&&VIX(i+1,2)==0
        VIX(i,3)=32;
        TrMVIX(3,2)=TrMVIX(3,2)+1;
    elseif VIX(i,2)==1&&VIX(i+1,2)==01
        VIX(i,3)=33;
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


SigSigSP=length(rSP);
T_SP=0;
T_BAC=0;
T_VIX=0;

for i=1:m
    for k=1:m
      T_BAC=T_BAC+2*TrMBAC(i,k)*log((TrMBAC(i,k)*sum(sum(TrMBAC)))/(sum(TrMBAC(i,:))*sum(TrMBAC(:,k))));
    end
end

for i=1:m
    for k=1:m
      T_SP=T_SP+2*TrMSP(i,k)*log((TrMSP(i,k)*SigSigSP)/(sum(TrMSP(i,:))*sum(TrMSP(:,k))));
    end
end

for i=1:m
    for k=1:m
        if TrMVIX(i,k)==0
            continue
        else
      T_VIX=T_VIX+2*TrMVIX(i,k)*log((TrMVIX(i,k)*sum(sum(TrMVIX)))/(sum(TrMVIX(i,:))*sum(TrMVIX(:,k))));
        end
    end
end

chi2cdf(T_BAC,4)
chi2cdf(T_SP,4)
chi2cdf(T_VIX,4)


alpha0=[rand,rand,rand];
alpha0=alpha0./sum(alpha0);
