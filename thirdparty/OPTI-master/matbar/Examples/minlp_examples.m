%%  A collection of mixed-integer nonlinear global optimization problems
% provided with the MATLAB/BARON interface from http://www.minlp.com/
% The collection was put together by J. Currie and N. Sahinidis
%#ok<*NASGU,*ASGLU,*NOPTS>

%% nvs01
clc
%Objective
fun = @(x) 0.04712385*(900 + x(1)^2)^0.5*x(2);
%Nonlinear Constraints
nlcon = @(x) [ 420.169404664517*(900 + x(1)^2)^0.5 - x(3)*x(1)*x(2);
               -x(3);
               (2960.87631843 + 18505.4769901875*x(2)^2)/(7200 + x(1)^2) - x(3)];           
cl = [0;-100;0];
cu = [0;Inf;Inf];
%Bounds
lb = [0;0;0];
ub = [200;200;100];
%Integer Variables
xtype = 'IIC';

%Solve
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu,xtype,[],baronset('filekp',0)) 

%% nvs02
clc
%Objective
fun = @(x) 9.99999999999999e-5*(5.3578547*x(3)^2 + 0.8356891*x(1)*x(5) + 37.293239*x(1)) + 5.9207859;
%Nonlinear Constraints
nlcon = @(x) [ -(0.0056858*x(2)*x(5) + 0.0006262*x(1)*x(4) - 0.0022053*x(3)*x(5)) + x(6);
               -(0.0071317*x(2)*x(5) + 0.0029955*x(1)*x(2) + 0.0021813*x(3)^2) + x(7);
               -(0.0047026*x(3)*x(5) + 0.0012547*x(1)*x(3) + 0.0019085*x(3)*x(4)) + x(8)];
cl = [85.334407;80.51249;9.300961];
cu = [85.334407;80.51249;9.300961];
%Bounds
lb = [0;0;0;0;0;0;90;20];
ub = [200;200;200;200;200;92;110;25];
%Integer Variables
xtype = 'IIIIICCC';
%Starting point
x0 = [NaN;NaN;NaN;NaN;NaN;NaN;90;20];

%Solve
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu,xtype,x0) 

%% nvs03
clc
%Objective
fun = @(x) (-8 + x(1))^2 + (-2 + x(2))^2;
%Linear Constraints
A = [-1/3 -1];
rl = -4.5;
ru = Inf;
%Nonlinear Constraints
nlcon = @(x) -0.1*x(1)^2 + x(2);
cl = 0;
cu = Inf;
%Bounds
lb = [0;0];
ub = [200;200];
%Integer Variables
xtype = 'II';

%Solve
[x,fval,ef,info] = baron(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xtype) 

%% nvs04
clc
%Objective
fun = @(x) 100*(0.5 - (0.6 + x(1))^2 + x(2))^2 + (0.4 - x(1))^2;
%Bounds
lb = [0;0];
ub = [200;200];
%Integer Variables
xtype = 'II';

%Solve
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,[],[],[],xtype) 

%% nvs05
clc
%Objective
fun = @(x) 1.10471*x(3)^2*x(4) + 0.04811*x(1)*x(2)*(14 + x(4));
%Linear Constraints
A = [0 1 -1 0 0 0 0 0];
rl = 0;
ru = Inf;
%Nonlinear Constraints
nlcon = @(x) [ -4243.28147100424/(x(3)*x(4)) + x(5);
               -(0.25*x(4)^2 + (0.5*x(1) + 0.5*x(3))^2)^0.5 + x(7);
               -(59405.9405940594 + 2121.64073550212*x(4))*x(7)/(x(3)*x(4)*(0.0833333333333333*x(4)^2 + (0.5*x(1) + 0.5*x(3))^2)) + x(6);
               -0.5*x(4)/x(7) + x(8);
               -(x(5)^2 + 2*x(5)*x(6)*x(8) + x(6)^2)^0.5;
               -504000/(x(1)^2*x(2));
               0.0204744897959184*(1e15*x(2)^3*x(1)*x(1)*x(2)^3)^0.5*(1 - 0.0282346219657891*x(1));
               -0.21952/(x(1)^3*x(2))];
cl = [0;0;0;0;-13600;-30000;6000;-0.25];
cu = [0;0;0;0;Inf;Inf;Inf;Inf];
%Bounds
lb = [1;1;0.01;0.01;-Inf;-Inf;-Inf;-Inf];
ub = [200;200;200;200;Inf;Inf;Inf;Inf];
%Integer Variables
xtype = 'IICCCCCC';
%Starting Point
x0 = [1;1;1;1;1;1;2;1];

%Solve
[x,fval,ef,info] = baron(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xtype,x0) 

%% Read Example Directly from .bar File
clc
clear
[fun,lb,ub,nlcon,cl,cu,xtype,x0] = baronRead('nvs01.bar');
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu,xtype,x0) 
