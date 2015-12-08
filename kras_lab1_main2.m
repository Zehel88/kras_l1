function varargout = kras_lab1_main2(varargin)
% KRAS_LAB1_MAIN2 MATLAB code for kras_lab1_main2.fig
%      KRAS_LAB1_MAIN2, by itself, creates a new KRAS_LAB1_MAIN2 or raises the existing
%      singleton*.
%
%      H = KRAS_LAB1_MAIN2 returns the handle to a new KRAS_LAB1_MAIN2 or the handle to
%      the existing singleton*.
%
%      KRAS_LAB1_MAIN2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KRAS_LAB1_MAIN2.M with the given input arguments.
%
%      KRAS_LAB1_MAIN2('Property','Value',...) creates a new KRAS_LAB1_MAIN2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kras_lab1_main2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kras_lab1_main2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kras_lab1_main2

% Last Modified by GUIDE v2.5 08-Dec-2015 16:10:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kras_lab1_main2_OpeningFcn, ...
                   'gui_OutputFcn',  @kras_lab1_main2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before kras_lab1_main2 is made visible.
function kras_lab1_main2_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
clc
% ====================================================================
syms p1 p2
% Заданная функция
% f=6*p1^2-5*p1*p2+2*p2^2+3*p1^3;
f=7*p1^2+4*p1*p2+p2^2-p2^3;
%% Нахождение критических точек заданной функции
% Частные производные
dfdp1=diff(f,p1);
dfdp2=diff(f,p2);
% Поиск нулей
Sres=solve(dfdp1,dfdp2);
xp1=double(Sres.p1);
xp2=double(Sres.p2);
% Выбор подходящих решений ( не 0 / не комплексные )
%% Построение линии уровня 
% Находим значение функции в критической точке
in_f=inline(f);
f_z=in_f(xp1(2),xp2(2));
% Cтроим линию уровня соответствующую найденному значению
x=-0.5:0.001:0.5;
y=-0.5:0.001:0.5;
[X,Y]=meshgrid(x,y);
Z=7.*X.^2+4.*X.*Y+Y.^2-Y.^3;
% Z=X.*Y.*-5.0+X.^2.*6.0+X.^3.*3.0+Y.^2.*2.0;

axes(handles.axes1)
contour(X,Y,Z,[f_z f_z]);grid on;hold on;

% Вписываем окружность, как многообразие с простым атласом
% CR=0.2;
CR=0.27;
plotCircle(0,0,CR);xlabel(['Многообразие и окружность M_0=S^1(0,',num2str(CR),')']);


%% Находим коэффициент пропорциональности
syms a x1 x2 y b z
% Ф-я в локальных координатах карты 1
fa=(collect(subs(f,[p1 p2],[y*a a*sqrt(CR^2-y^2)]),a))-f_z;


[c,p]=coeffs(fa,a);
np=find(sym2poly(sum(p))==1);


% Нахождение коэф. пропорц. для 1 карты
n=50;


j=1;
for i=-CR:(2*CR/n):CR
%   Замена символьных переменных 
    w(np)=double(subs(c,y,i));
%   Нахождение корней многочлена
    r=(roots(w));
%   Нахождение корней >0
    r=r(find(r>0));
%   Нахождение корней <1
%     r=r(find(r<1));

%   Нахождение действительных корней  
a(j)=min(r);
%     try
%     a(j)=r(find(imag(r)==0));
%     catch
%     j=j-1;
%     end
%   буферные переменные
    k1(j)=i;
    j=j+1;
end
% Построение гарфика коэф. проп 1 карты
axes(handles.axes6)
plot(k1,a);grid on;xlabel({'Коэф. пропорциональности ','для карты (U_{1},\phi_{1})'});

%% Ф-я в локальных координатах карты 2
fb=subs(f,[p1 p2],[b.*sqrt(CR^2-z.^2) b.*z])-f_z;

[c,p]=coeffs(fb,b);
np=find(sym2poly(sum(p))==1);

% Нахождение коэф. пропорц. для 2 карты
j=1;
for i=CR:-(2*CR/n):-CR
%  Замена символьных переменных 
    w(np)=double(subs(c,z,i));

%   Нахождение корней многочлена
    r=(roots(w));
%   Нахождение корней >0
    r=r(find(r>0));
%   Нахождение корней <1
%     r=r(find(r<1));
%   Нахождение действительных корней  
b(j)=min(r);
%     try
%     b(j)=r(find(imag(r)==0));
%     catch
%     j=j-1;
%     end
%   буферные переменные
    k2(j)=i;
    j=j+1;
end
% Построение гарфика коэф. проп 2 карты
axes(handles.axes7)
plot(k2,b);grid on;xlabel({'Коэф. пропорциональности ','для карты (U_{2},\phi_{2})'});

%% Построение функции в координатном пространстве {G1(y),y}
axes(handles.axes8)
% fy=sqrt(CR^2-y.^2);

for i=1:numel(k1)
% G1(i)=a1(i)^4*fy(i)^3*sqrt(CR^2-fy(i)^2)+a1(i)^3*(CR^2-fy(i)^2)*fy(i)-a1(i)^5*(CR^2-fy(i)^2)^(5/2);
p1=a(i)*k1(i);
p2=a(i)*sqrt(CR^2-k1(i)^2);
G1(i)=7*p1^2+4*p1*p2+p2^2-p2^3;
end

plot(k1,G1);grid on;
xlabel({'График функции F(p)',' в координатах карты (U_{1},\phi_{1})'});
hold on

%% Построение функции в координатном пространстве {G2(z),z}

% fz=sqrt(CR^2-z.^2);
for i=1:numel(k2)
% G2(i)=a2(i)^4*fz(i)*(CR^2-fz(i)^2)^(3/2)+a2(i)^3*fz(i)^2*sqrt(CR^2-fz(i)^2)-a2(i)^5*fz(i)^5;
p1=b(i)*(sqrt(CR^2-k2(i)^2));
p2=b(i)*k2(i);
G2(i)=7*p1^2+4*p1*p2+p2^2-p2^3;
end

plot(k2,G2);grid on;
xlabel({'График функции F(p)' ,'в координатах карты  (U_{2},\phi_{2})'})
%% 

axes(handles.axes2)
y=CR.*cos(pi:-pi/n:0);
numel(a)
numel(y)
p_1=a.*y;
p_2=a.*(sqrt(CR^2-y.^2));
G_1=7.*p_1.^2+4.*p_1.*p_2+p_2.^2-p_2.^3;

plot(y,G_1);
grid on

axes(handles.axes3)
z=CR.*sin(pi/2:-pi/n:-pi/2);

p__1=b.*(sqrt(CR^2-z.^2));
p__2=b.*z;
figure(1)
plot(p__1)
figure(2)
plot(p__2)
G_2=7.*p__1.^2+4.*p__1.*p__2+p__2.^2-p__2.^3;

plot(z,G_2);
grid on

DaTa.a=a;
DaTa.b=b;
save('DaTa.mat','DaTa');

%% 
% Функция для задания 3
% x_range=-1:0.05:1;
% y_range=-1:0.05:1;
% [p31,p32] = meshgrid(x_range,y_range);

% f3=p31.^3.*p32+p32.^2.*p31-p32.^5;
% % f3=6.*p31.^2-5.*p31.*p32+2.*p32.^2+3.*p31.^3;
% axes(handles.axes5)
% mesh(p31,p32,f3,'UIContextMenu',f2pCM2);
% xlabel('График функции F(p)                     ')

% DaTa.X=X;
% DaTa.Y=Y;
% DaTa.Z=Z;

% save('DaTa.mat','DaTa');


% UIWAIT makes kras_lab1_main2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function plotCircle (xc, yc, R)
plot(xc + R * cos(0:0.001:2*pi), yc + R * sin(0:0.001:2*pi));


function cmHandle = f2pCM2(n_con)
   cmHandle = uicontextmenu;
   uimenu(cmHandle,'Label','Построить в отдельном окне','Callback',@newfigureplot2);
   
function newfigureplot2(x1,x2)

h=gco;
figure;
if (isempty(h.ZData)==0)
mesh(h.XData,h.YData,h.ZData);
else
    load('f2p.mat');
    load('Pref.mat');
    contour(h.XData,h.YData,f2p.func2plot,Pref{3,2});hold on;
    quiver(h.XData,h.YData,h.UData,h.VData);hold on;
    gr = gradient(f2p.func2plot,.5,.5);

[zero_grad_ind_x,zero_grad_ind_y]=find(gr==0);
x_range=Pref{4,2};y_range=Pref{5,2};

     for i=1:numel(zero_grad_ind_y)
         plot(x_range(zero_grad_ind_y(i)),y_range(zero_grad_ind_x(i)),'*');
         hold on;
     end
    hold off;
end

% --- Outputs from this function are returned to the command line.
function varargout = kras_lab1_main2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function NFP_Callback(hObject, eventdata, handles)
% hObject    handle to NFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles,X)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% -------------------Theory------------------------
function Teo_Callback(hObject, eventdata, handles)
% Открываем док MuPad
open('Kras_lab1_teo.doc');
