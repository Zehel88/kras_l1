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

% Last Modified by GUIDE v2.5 10-Dec-2015 21:54:12

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
f=6*p1^2-5*p1*p2+2*p2^2+3*p1^3;
% f=7*p1^2+4*p1*p2+p2^2-p2^3;
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
x=-1:0.05:1;
y=-1:0.05:1;
[X,Y]=meshgrid(x,y);
Z=X.*Y.*-5.0+X.^2.*6.0+X.^3.*3.0+Y.^2.*2.0;

axes(handles.axes1)
contour(X,Y,Z,[f_z f_z]);grid on;hold on;

% Описываем окружностью, как многообразием с простым атласом
CR=0.2;
plotCircle(0,0,CR);
xlabel({'x_1,p_1',['Многообразие и окружность M_0=S^1(0,',num2str(CR),')']});
ylabel('x_2,p_2');
DaTa.X=X; DaTa.Y=Y; DaTa.Z=Z; DaTa.f_z=f_z; 
DaTa.CR=CR;

%% Находим коэффициент пропорциональности
% Кол-во точек
n=50;
% ОО параметра t
t=pi:-pi/n:0;
x1=CR*cos(t);
x2=CR*sin(t);

% Нахождение коэф. пропорц. для 1 карты
for i=1:length(t)
% Значение функции для конкретной точки
    A=6.*x1(i).^2-5.*x1(i).*x2(i)+2.*x2(i).^2;
    B=3.*x1(i).^3;

if B<=10^-5 || A<=10^-5
    B=0;
end
% Функция в конкретной точке (пол. вид)
    fAB=[B A 0 -f_z];
% Нахождение корней полинома
    r=roots(fAB);
% Определение коэф. пропорциональности
for j=1:length(r)
    if  imag(r(j))==0 && r(j)>0
       a(i)=r(j);
    end
 end
end

y=x1;

% Построение гарфика коэф. проп 1 карты
axes(handles.axes6)
plot(y,a);grid on;
xlabel({'y','Коэф. пропорциональности ','для карты (U_{1},\phi_{1})'});
ylabel('\alpha(y)');
%% Построение функции в координатном пространстве {G1(y),y}
axes(handles.axes8)

p1=a.*y;
p2=a.*sqrt(CR^2-y.^2);
G1=6.*p1.^2-5.*p1.*p2+2.*p2.^2+p1.^3;

plot(y,G1);grid on;
xlabel({'y','График функции F(p)',' в координатах карты (U_{1},\phi_{1})'});
ylabel('G_1(y)');



%% Ф-я в локальных координатах карты 2
% ОО параметра t
t=-pi/2:pi/n:pi/2;
x1=CR*cos(t);
x2=CR*sin(t);
% Нахождение коэф. пропорц. для 1 карты
for i=1:length(t)
% Значение функции для конкретной точки
    A=6.*x1(i).^2-5.*x1(i).*x2(i)+2.*x2(i).^2;
    B=3.*x1(i).^3;

if B<=10^-5 || A<=10^-5
    B=0;
end
% Функция в конкретной точке (пол. вид)
    fAB=[B A 0 -f_z];
% Нахождение корней полинома
    r=roots(fAB);
% Определение коэф. пропорциональности
for j=1:length(r)
    if  imag(r(j))==0 && r(j)>0
       b(i)=r(j);
    end
 end
end
z=x2;
% Построение гарфика коэф. проп 2 карты
axes(handles.axes7)
plot(z,b);grid on;xlabel({'z','Коэф. пропорциональности ','для карты (U_{2},\phi_{2})'});
ylabel('\beta(z)');


%% Построение функции в координатном пространстве {G2(z),z}

p1=b.*sqrt(CR^2-z.^2);
p2=b.*z;
G2=6.*p1.^2-5.*p1.*p2+2.*p2.^2+p1.^3;
axes(handles.axes2)
plot(z,G2);grid on;
xlabel({'z','График функции F(p)' ,'в координатах карты  (U_{2},\phi_{2})'})
ylabel('G_2(z)');
%% Построение графиков функций f(p) в координатных пространствах {G2(z),z},{G1(y),y}
axes(handles.axes3)
plot(z,G2);grid on;hold on;
plot(y,G1);
xlabel({'y,z','Графики f(p)' ,'в коорд. просранствах  (G2(z),z),(G1(y),y)'})
ylabel('G_1(y),G_2(z)');

%% Сохранение необходимых данных 
DaTa.a=a; DaTa.b=b;
DaTa.y=y; DaTa.z=z;
DaTa.G1=G1; DaTa.G2=G2;
save('DaTa.mat','DaTa');



function plotCircle (xc, yc, R)
plot(xc + R * cos(0:0.001:2*pi), yc + R * sin(0:0.001:2*pi));


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

%% Построение на отдельных фаормах 
% ----------------Puc 1-----------------------------------
function Untitled_6_Callback(hObject, eventdata, handles,X)
load('DaTa.mat');
figure(1)
contour(DaTa.X,DaTa.Y,DaTa.Z,[DaTa.f_z DaTa.f_z]);grid on;hold on;
plotCircle(0,0,DaTa.CR);
xlabel({'x_1,p_1',['Многообразие и окружность M_0=S^1(0,',num2str(DaTa.CR),')']});
ylabel('x_2,p_2');

% ----------------------Puc 2-----------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
load('DaTa.mat');
figure(2)
plot(DaTa.y,DaTa.a);grid on;
xlabel({'y','Коэф. пропорциональности ','для карты (U_{1},\phi_{1})'});
ylabel('\alpha(y)');

% ---------------------Puc 3-------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
load('DaTa.mat');
figure(3)
plot(DaTa.z,DaTa.b);grid on;xlabel({'z','Коэф. пропорциональности ','для карты (U_{2},\phi_{2})'});
ylabel('\beta(z)');

% -------------------Theory------------------------
function Teo_Callback(hObject, eventdata, handles)
% Открываем док MuPad
open('Kras_lab1_teo.doc');


% ----------------------Puc 4----------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
load('DaTa.mat');
figure(4)
plot(DaTa.y,DaTa.G1);grid on;
xlabel({'y','График функции F(p)',' в координатах карты (U_{1},\phi_{1})'});
ylabel('G_1(y)');

% ----------------------Puc 5----------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
load('DaTa.mat');
figure(5)
plot(DaTa.z,DaTa.G2);grid on;
xlabel({'z','График функции F(p)',' в координатах карты (U_{2},\phi_{2})'});
ylabel('G_2(z)');


% -----------------------Puc 6-------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
load('DaTa.mat');
figure(6)
pl1=plot(DaTa.z,DaTa.G2);grid on;hold on;
plot(DaTa.y,DaTa.G1);
xlabel({'y,z','Графики f(p)' ,'в коорд. просранствах  (G2(z),z),(G1(y),y)'})
ylabel('G_1(y),G_2(z)');


function Untitled_12_Callback(hObject, eventdata, handles)
%% Part 2 
syms x1 x2 x3 t
eq1=x1^2+x1*x2*x3;
eq2=x2^2*sin(3*x1)+x3;
eq3=x1*cos(x3)-x3^3;
% Нахождение равновесной точки
soeq=fsolve(@sysfun,[-5,-5,-5])
% Якобиан системы
J=[diff(eq1,x1) diff(eq1,x2) diff(eq1,x3);
   diff(eq2,x1) diff(eq2,x2) diff(eq2,x3);
   diff(eq3,x1) diff(eq3,x2) diff(eq3,x3)]

for i=1:3
a{i,1}=(int(subs(J(i,1),[x1,x2,x3],[(soeq(1)+t*(x1-soeq(1))) (soeq(2)+t*(x2-soeq(2))) (soeq(3)+t*(x3-soeq(3)))]),t,0,1));
a{i,2}=(int(subs(J(i,1),[x1,x2,x3],[(soeq(1)+t*(x1-soeq(1))) (soeq(2)+t*(x2-soeq(2))) (soeq(3)+t*(x3-soeq(3)))]),t,0,1));
a{i,3}=(int(subs(J(i,1),[x1,x2,x3],[(soeq(1)+t*(x1-soeq(1))) (soeq(2)+t*(x2-soeq(2))) (soeq(3)+t*(x3-soeq(3)))]),t,0,1));
end
save('a.mat','a');


% Данная система нел.ур
function res=sysfun(x)
res=[x(1)^2+x(1)*x(2)*x(3);
x(2)^2*sin(3*x(1))+x(3);
x(1)*cos(x(3))-x(3)^3];


