%%Computing reduct from hybrid data in Neighbourhood information system.
%%Neighborhood attribute quality degree is employed as the heuristic rule.
%%Please refer to the following article.
%%Yi Li, Shenhong Lei, Xinyu Su, Zhong Yuan, Ji Li, Xingqiang Tan
%%Unsupervised Hybrid Attribute Selection Based on Variable Precision
%%Neighborhood Rough Sets, Expert Systems with Applications.
%%Uploaded by Yi Li on February 13, 2026. 
%%E-mail:scmzly@scun.edu.cn

function select_feature=uar_NAQD(data,k)
%% input
% data is data matrix without decision, where rows for objects and columns for attributes. 
% Numerical attributes should be normalized into [0,1]
% Variable precision threshold parameter beta [0.5,1.0] with step size 0.1.
%% output
% a reduct--the set of selected attributes.

%% Calculate the standard deviation of the numerical data
[row, attribute]=size(data);
delta=zeros(1,attribute);
for j=1:attribute
    if min(data(:,j))==0 && max(data(:,j))==1
       delta(j)=std(data(:,j),1);
    end
end

%% Calculate the neighborhood relation matrices ssk and ssr.
for i=1:attribute
       eval(['ssk' num2str(i) '=[];']);
       dis=pdist2(data(:,i),data(:,i));
       if delta(i)>0
           r0=dis<=delta(i);
       else
           r0=dis==0;
       end
       r=r0;
       eval(['ssk' num2str(i) '=r;']);
end

for i=1:attribute
       eval(['ssr' num2str(i) '=[];']);
       dis=pdist2(data(:,i),data(:,i));
       if delta(i)>0
           r0=dis<=delta(i);
       else
           r0=dis==0;
       end
       r=r0;
       eval(['ssr' num2str(i) '=r;']);
end

%% Algorithms topic section
B=1:attribute;
x=0;
value_sig=0;
value_inc=0;
red=[];
for l=1:attribute
sig=[];
inc=[];
for l_1=1:length(B)
    r1=eval(['ssk' num2str(B(l_1))]);
    for l_2=1:attribute
        r_SIN=eval(['ssr' num2str(l_2)]);
        [r_SIN_temp,~,r_SIN_ic]=unique(r_SIN,'rows');
        m=zeros(row);
        ks=zeros(row);
        for l_3=1:size(r_SIN_temp,1)
            sxx_Min=min(r1,repmat(r_SIN_temp(l_3,:),row,1));
            for l_4=1:row
                value = sum(sxx_Min(l_4,:),2)/sum(r1(l_4,:),2);
                if value>=k
                    m(l_3,l_4)=1;
                    ks(l_3,l_4)=value;
                end
            end
        end
        importance_SIN=sum(max(m,[],1));
        importance_INC=sum(max(ks,[],1));
        if sum(max(ks,[],1))==0
            inc(l_1,l_2)=0;
        else
            inc(l_1,l_2)=importance_INC/importance_SIN;
        end
        sig(l_1,l_2)=importance_SIN/row; 
    end
end
[m0,n0]=size(sig);
incl = inc - repmat(value_inc,m0,1);
incl(find(incl==0))=[0.0001];
sigl = sig - repmat(value_sig,m0,1);
Aqd=incl.*sigl;
[x1,n1]=max(mean(Aqd,2));
x11=mean(sig(n1,:));
x=[x;x11];
len=length(x);
value_sig=sig(n1,:);
value_inc=inc(n1,:);
if abs(x(len)-x(len-1))>0.0001
    red=[red B(n1)];
    B=setdiff(B,B(n1));
    for j2=1:length(B)
         B1=[red B(j2)];
         rk = Nmixture(delta,data,B1);
         eval(['ssk' num2str(B(j2)) '=rk;']);
    end
else
    break
end

end

if length(red)==attribute

   select_feature=red(1:end-1);
else
   select_feature=red;
end
end



