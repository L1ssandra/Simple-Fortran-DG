h = figure();				% 创建图形窗口
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% 关闭相关的警告提示（因为调用了非公开接口）
jFrame = get(h,'JavaFrame');	% 获取底层 Java 结构相关句柄
pause(0.1);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
pause(0.1);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% 打开相关警告设置

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