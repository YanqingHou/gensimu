function [xsv,ysv,zsv,elsv,no_sv,rows]=choosesats(xs,ys,zs,plh,xyz,cutoff,num)
% 返回值：rows，所有no_sv中选num的组合数
%       no_sv 当前可见卫星个数
%       elsv  所有组合的卫星高度角矩阵，每一列为一个组合，共有rows列
%       xsv，ysv，zsv  所有组合的卫星x y z位置，每一列为一个组合，共有rows列，计rows个组合
el       = cpelev1(xs,ys,zs,plh,xyz);
idx1     = find (el > cutoff);
no_sv    = length(idx1);

if num>no_sv
    xsv=0;
    ysv=0;
    zsv=0;
    elsv=0;
    rows=0;
    return;
end
C=nchoosek(idx1,num);%每一行是一个选择
rows=size(C,1);
xsv=zeros(num,rows);
ysv=xsv; zsv=xsv; elsv=xsv;
for i=1:rows
    idx=C(i,:)';
    xsv(:,i)      = xs(idx);
    ysv(:,i)       = ys(idx);
    zsv(:,i)      = zs(idx);
    elsv(:,i)       = el(idx);
end
% refsat   = find(el(idx1)==max(el(idx1)));
