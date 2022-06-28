% 显示解随时间变化规律,从Fortran写好的文件中读取数据
% flashFortran.m
%xa = 0; xb = 2*pi; ya = 0; yb = 2*pi; % Tang Vortex
%xa = -10; xb = 10; ya = -10; yb = 10; % 平移
%xa = 0;xb = 1;ya = 0;yb = 1;
xa = 0;xb = 10;ya = 0;yb = 10; % 等熵涡
TIME = 2000;
Q1 = str2num(fileread('Q2.txt'));
Q1 = Q1(:,1);
T = fileread('T.txt');
T = str2num(T);

STOP = length(T);

N = round(sqrt(length(Q1)/STOP));
% 网格
h = (xb - xa)/N;
h1 = h/2;
X = zeros(N,1);
Y = zeros(N,1);
for i = 1:N + 1
    X(i) = xa + (i - 1)*h;
    Y(i) = ya + (i - 1)*h;
end

% 整格点
Xc = (X(1:end - 1) + X(2:end))/2;
Yc = (Y(1:end - 1) + Y(2:end))/2;




Q1 = reshape(Q1,N,N,STOP);
stop = 0;

%QF = QF(:,:,(1:STOP));
T = T(1:STOP);



h = figure();
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(h,'JavaFrame');
pause(0.1);
set(jFrame,'Maximized',1);
pause(0.1);
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');


[xc,yc] = meshgrid(Xc,Yc);
xlabel('X');
ylabel('Y');
grid on;
%s = [X(1),X(end),Y(1),Y(end),min(min(min(Q(:,:,:,4)))) - 0.1,max(max(max(Q(:,:,:,4)))) + 0.1];
s = [Xc(1),Xc(end),Yc(1),Yc(end),min(min(min(Q1))) - 0.1,max(max(max(Q1))) + 0.1];
axis(s);
TT = 100;
t0 = T(end)/TT;
for i = 1:101
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - tt));
    mesh(yc,xc,Q1(:,:,j));
    %contour(yc,xc,Q1(:,:,j,1));
    title(T(j));
    colormap(cool)
    axis(s);
    pause(0.0001);
end
mesh(yc,xc,Q1(:,:,j))
axis(s)
%plot(T,DIV);
% [~,j] = min(abs(T - 0.5));
% contour(yc,xc,QF(:,:,j,1),25);
% [~,j] = min(abs(T - 2));
% contour(Xc,Yc,QF(:,:,j,1),25);
% [~,j] = min(abs(T - 3));
% contour(Xc,Yc,QF(:,:,j,1),25);
% [~,j] = min(abs(T - 4));
% contour(Xc,Yc,QF(:,:,j,1),25);