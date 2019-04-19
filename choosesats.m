function [xsv,ysv,zsv,elsv,no_sv,rows]=choosesats(xs,ys,zs,plh,xyz,cutoff,num)
% ����ֵ��rows������no_sv��ѡnum�������
%       no_sv ��ǰ�ɼ����Ǹ���
%       elsv  ������ϵ����Ǹ߶ȽǾ���ÿһ��Ϊһ����ϣ�����rows��
%       xsv��ysv��zsv  ������ϵ�����x y zλ�ã�ÿһ��Ϊһ����ϣ�����rows�У���rows�����
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
C=nchoosek(idx1,num);%ÿһ����һ��ѡ��
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
