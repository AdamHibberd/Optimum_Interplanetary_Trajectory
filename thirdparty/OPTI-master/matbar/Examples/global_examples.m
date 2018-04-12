%% A collection of nonlinear global optimization problems
% provided with the MATLAB/BARON interface from http://www.minlp.com/
% The collection was put together by J. Currie and N. Sahinidis
%#ok<*NASGU,*ASGLU,*NOPTS>

%% e01
clc
%Objective
fun = @(x) -x(1) - x(2);
%Nonlinear Constraints
nlcon = @(x) prod(x);
cl = -Inf;
cu = 4;
%Bounds
lb = [0;0];
ub = [6;4];

%Solve
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu) 

%% e02
clc
%Objective
fun = @(x) x(3);
%Nonlinear Constraints
nlcon = @(x) [ 30*x(1) - 6*x(1)*x(1) - x(3);
               20*x(2) - 12*x(2)*x(2) - x(3);
               0.5*(x(1) + x(2))^2 - x(3) ];
cl = [-250;-300;-150];
cu = [-250;-300;-150];
%Bounds
lb = [0;0;0];
ub = [9.422;5.9023;267.417085245];

%Solve
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu)

%% e03
clc
%Objective
fun = @(x) -0.063*x(4)*x(7) + 5.04*x(1) + 0.035*x(2) + 10*x(3) + 3.36*x(5);
%Linear Constraints
A = [1 0 0 -1.22 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0.222;
     0 0 0 0 0 0 3 0 0 -1];
rl = [0;35.82;133];
ru = [0;35.82;133];                        
%Nonlinear Constraints
nlcon = @(x) [ 0.038*x(8)^2 - 1.098*x(8) - 0.325*x(6) + x(7);
               x(4)*x(9)*x(6) + 1000*x(3)*x(6) - 98000*x(3);
               -x(1)*x(8) + x(2) + x(5);
               0.13167*x(8)*x(1) + 1.12*x(1) - 0.00667*x(8)^2*x(1) - x(4)];
cl = [57.425;0;0;0];
cu = [57.425;0;0;Inf];
%Bounds
lb = [1;1;0;1;0;85;90;3;1.2;145];
ub = [2000;16000;120;5000;2000;93;95;12;4;162];
%Starting point
x0 = [1;1;NaN;1;NaN;85;90;3;1.2;145];

%Solve
[x,fval,ef,info] = baron(fun,A,rl,ru,lb,ub,nlcon,cl,cu,[],x0)

%% e04
clc
%Objective
fun = @(x) 400*x(1)^0.9 + 22*(-14.7 + x(2))^1.2 + x(3) + 1000;
%Nonlinear Constraints
nlcon = @(x) [ x(3)*x(1) + 144*x(4);
               -exp(11.86 - 3950/(460 + x(4))) + x(2) ];
cl = [11520;0];
cu = [Inf;0];
%Bounds
lb = [0;14.7;0;-459.67];
ub = [15.1;94.2;5371;80];
%Starting point
x0 = [NaN;14.7;NaN;NaN];

%Solve
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,baronset('filekp',1,'cplexlibname','cplex1260.dll'))

%% e05
clc
%Objective
fun = @(x) sum(x(1:3));
%Linear Constraints
A = [0 0 -4e3 0 -1e5];
rl = -5e7;
ru = -5e7;
%Nonlinear Constraints
nlcon = @(x) [ 1e5*x(4) - 120*x(1)*(300-x(4));
               1e5*x(5) - 80*x(2)*(400-x(5)) - 1e5*x(4)];
cl = [1e7;0];
cu = [1e7;0];
%Bounds
lb = [0;0;0;100;100];
ub = [15834;36250;10000;300;400];
%Starting point
x0 = [NaN;NaN;NaN;100;100];

%Solve
[x,fval,ef,info] = baron(fun,A,rl,ru,lb,ub,nlcon,cl,cu,[],x0)

%% Read Example Directly from .bar File
clc
clear
[fun,lb,ub,nlcon,cl,cu,xint,x0] = baronRead('e01.bar')
[x,fval,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu,xint,x0) 
