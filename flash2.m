h = figure();				% ����ͼ�δ���
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% �ر���صľ�����ʾ����Ϊ�����˷ǹ����ӿڣ�
jFrame = get(h,'JavaFrame');	% ��ȡ�ײ� Java �ṹ��ؾ��
pause(0.1);					% �� Win 10��Matlab 2017b �����²���ͣ�ٻᱨ Java �ײ���󡣸��˸�����Ҫ���Խ���ʵ����֤
set(jFrame,'Maximized',1);	%���������Ϊ�棨0 Ϊ�٣�
pause(0.1);					% ����ʵ���з��������ͣ�٣����ڿ����������仯������ȡ�Ĵ��ڴ�С����ԭ���ĳߴ硣���˸�����Ҫ���Խ���ʵ����֤
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% ����ؾ�������

s = [Xc(1),Xc(end),Yc(1),Yc(end),min(min(min(Q1))) - 0.1,max(max(max(Q1))) + 0.1];
[xc,yc] = meshgrid(Xc,Yc);
stopL = 0;
STOPL = 0;
stopR = 0;
STOPR = 0;
for i = 1:N
    if Xc(i) > 0.3 && stopL == 0
        stopL = 1;
        STOPL = i;
    end
    if Xc(i) > 0.7 && stopR == 0
        stopR = 1;
        STOPR = i;
    end
end
TT = 100;
t0 = T(end)/TT;
%for k = 1:10
for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    mesh(yc,xc,Q1(:,:,j));
    axis(s)
    colormap(cool)
    %contour(yc,xc,Q1(:,:,j),40);
    title(T(j))
    pause(0.0001);
end
%end
%[~,j] = min(abs(T - 3));
%mesh(yc,xc,QF(:,:,j));
%axis(s);
% [xc,yc] = meshgrid(Xc(STOPL:STOPR),Yc(STOPL:STOPR));
% [~,j] = min(abs(T - 4));
% colormap(cool)
% contour(yc,xc,Q1(STOPL:STOPR,STOPL:STOPR,j),30);