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

% Last Modified by GUIDE v2.5 07-Dec-2015 11:24:10

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

% ====================================================================
syms p1 p2
% Заданная функция
f=6*p1^2-5*p1*p2+2*p2^2+3*p1^3;
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
% Вписываем окружность, как многообразие с простым атласом
CR=0.2;
plotCircle(0,0, CR);xlabel(['Многообразие и окружность M_0=S^1(0,',num2str(CR),')']);
%% Находим коэффициент пропорциональности
syms a1 x1 x2 y1 a2 z1
% Ф-я в локальных координатах карты 1
fa1=subs(f,[p1 p2],[a1.*y1 a1.*sqrt(CR^2-y1.^2)]);
% Нахождение коэф. пропорц. для 1 карты
j=1;
for i=-CR:(2*CR/50):CR
    fi1=subs(fa1-CR,y1,i);
    r=roots(sym2poly(fi1));
    t1=find(r>0);
    a1(j)=min(r(t1));
    k1(j)=i;
    j=j+1; 
end
% Ф-я в локальных координатах карты 2
fa2=subs(f,[p1 p2],[a2.*sqrt(CR^2-z1.^2) a2.*z1]);
% Нахождение коэф. пропорц. для 2 карты
j=1;
for i=-CR:(2*CR/50):CR
    fi2=subs(fa2-CR,z1,i);
    r=roots(sym2poly(fi2));
    t1=find(r>0);
    a2(j)=min(r(t1));
    k2(j)=i;
    j=j+1;   
end
% Построение гарфика коэф. проп 1 карты
axes(handles.axes6)
plot(k1,a1);grid on;xlabel({'Коэф. пропорциональности ','для карты (U_{1},\phi_{1})'});
% Построение гарфика коэф. проп 2 карты
axes(handles.axes7)
plot(k2,a2);grid on;xlabel({'Коэф. пропорциональности ','для карты (U_{2},\phi_{2})'});

% Построение функции в координатном пространстве {G1(y),y}
axes(handles.axes2)
fy=-1:2/50:1;
% fy=sqrt(CR^2-y.^2);
for i=1:numel(fy)
G1(i)=a1(i)^4*fy(i)^3*sqrt(CR^2-fy(i)^2)+a1(i)^3*(CR^2-fy(i)^2)*fy(i)-a1(i)^5*(CR^2-fy(i)^2)^(5/2);
end
% G1=6*a^2.*y.^2-5*a^2.*y.*sqrt(CR^2-y.^2)+2*a^2.*(CR^2-y.^2)+3*a^3.*y.^3;
plot(fy,G1);grid on;
xlabel({'График функции F(p)',' в координатах карты (U_{1},\phi_{1})'})

% Построение функции в координатном пространстве {G2(z),z}
axes(handles.axes3)
fz=-1:2/50:1;
% fz=sqrt(CR^2-z.^2);
for i=1:numel(fz)
G2(i)=a2(i)^4*fz(i)*(CR^2-fz(i)^2)^(3/2)+a2(i)^3*fz(i)^2*sqrt(CR^2-fz(i)^2)-a2(i)^5*fz(i)^5;
end
% G2=6*a^2.*(CR^2-z.^2)-5*a^2.*(CR^2-z.^2).*z+2*a^2.*z.^2+3*a^3.*(CR^2-z.^2).^(3/2);
plot(fz,G2);grid on;
xlabel({'График функции F(p)' ,'в координатах карты  (U_{2},\phi_{2})'})

% Функция для задания 3
x_range=-1:0.05:1;
y_range=-1:0.05:1;
[p31,p32] = meshgrid(x_range,y_range);

f3=p31.^3.*p32+p32.^2.*p31-p32.^5;
% f3=6.*p31.^2-5.*p31.*p32+2.*p32.^2+3.*p31.^3;
axes(handles.axes5)
mesh(p31,p32,f3,'UIContextMenu',f2pCM2);
xlabel('График функции F(p)                     ')


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


% -----------------------INFOBTN------------------------------
function infobtn_ClickedCallback(hObject, eventdata, handles)

