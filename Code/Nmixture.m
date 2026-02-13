function r=Nmixture(delta,data,B1)
ID=all(data(:,B1)<=1);
if sum(ID)==0 %Numerical 
    dis1=pdist2(data(:,B1),data(:,B1));
    rk=dis1==0;
elseif sum(ID)==length(ID) %Nominal
    dis2=pdist2(data(:,B1),data(:,B1));
    a01 = mean(delta(:,B1));
    rk=dis2<=a01;
else %hybrid
    [nL,mL]=find(ID==0);
    BL=B1(:,mL);
    dis_L = pdist2(data(:,BL),data(:,BL));
    rk_L = dis_L==0;

    [nS,mS]=find(ID==1);
    BS=B1(:,mS);
    dis_S = pdist2(data(:,BS),data(:,BS));
    a02 = mean(delta(:,BS));
    rk_S = dis_S<=a02;
    rk = min(rk_L,rk_S);
end
r = rk;
end