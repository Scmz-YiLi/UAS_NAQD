clear all
clc
load Example.mat
trandata=data;
trandata(:,1:3)=normalize(trandata(:,1:3),'range');
k=0.6;
select_feature=uar_NAQD(trandata,k)

