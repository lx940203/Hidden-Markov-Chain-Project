
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


