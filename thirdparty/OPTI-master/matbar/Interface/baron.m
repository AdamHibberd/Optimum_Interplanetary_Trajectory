%BARON Solve a NLP/MINLP using the global MINLP solver BARON
%   
%   Calling Form:
%       [x,fval,exitflag,info,allsol] = baron(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xtype,x0,opts)
%
%   x = baron(fun,A,rl,ru) solves the constrained nonlinear system using 
%   the Global MINLP solver BARON. fun must be a function 
%   specified using BARON rules (see below). A, rl and ru specify the
%   linear constraints of the form rl <= Ax <= ru. To specify a one sided
%   constraint use infinity on the empty side, and rl = ru for equality
%   constraints.
%
%   x = baron(fun,A,rl,ru,lb,ub) solves subject to the decision
%   variable bounds, lb and ub. Use -inf for unbounded lb, and inf for
%   unbounded ub.
%
%   x = baron(fun,...,ub,nlcon,cl,cu) solves subject to the nonlinear
%   constraints specified by nlcon, cl and cu of the form 
%   cl <= nlcon(x) <= cu. nlcon is a column-wise vector function and must 
%   be a function specified using BARON rules. To specify a one
%   sided constraint use infinity on the empty side, and cl = cu for an
%   equality constraint.
%
%   x = baron(fun,...,cu,xtype) solves subject to the integer and binary
%   constraints specified. xtype is a char array where 'C' is continuous,
%   'I' is integer and 'B' is binary.
%
%   x = baron(fun,...,xtype,x0) allows the user to specify an initial
%   solution guess via x0. This is an optional field as BARON can generate
%   its own starting point. To specify a partial starting point, fill
%   unknown values in the x0 vector with NaNs.
%
%   x = baron(fun,...,x0,opts) allows the user to specify options
%   specific to BARON via opts. Use baronset() to generate this options
%   structure.
%
%   [x,fval,exitflag,info] = baron(fun,...,opts) also returns the objective
%   value at the solution, the solver exitflag and an information structure
%   on the solver progress.
%
%   [x,fval,exitflag,info,allsol] = baron(...) also returns a structure
%   containing all requested solutions and objective values when numsol >
%   1.
%
%
%   Matlab functions supplied to this method are processed and passed to
%   BARON as a text model and therefore must meet certain requirements. An 
%   example function is shown below:
%
%       obj = @(x) 3*x(1)^2 - 2*x(2)^3 + 3*x(1)*x(2);
%
%       OR
%
%       function fx = objective(x)
%          fx = 3*x(1)^2 - 2*x(2)^3 + 3*x(1)*x(2);
%       end
%
%   The function must be:
%       - A function of one input argument (the decision variable vector)
%       - Contain scalar, vector, or matrix operations only (max 2D supported)
%       - Limited to a subset of functions (exp, log, log10, dot, sqrt, norm, sum, prod)
%       - Contain no conditional statements or relational operators with respect to the decision variables
%
%   For nonlinear constraints the function must be a column-wise vector:
%
%       con = @(x) [-(x(1)*x(2));
%                   3*x(1) + 2*x(2)];
%
%       OR
%
%       function cx = nlconstraints(x)
%          cx(1) = -(x(1)*x(2));
%          cx(2) = 3*x(1) + 2*x(2);
%       end
%
%
%   Note on Temporary Installation:
%   If you are working on a computer where you do not have administrator
%   privileges and wish to use BARON you can copy the p-code (.p) files and
%   the BARON executable in the BARON/Interface/ folder to the directory
%   where you wish to use BARON from.
%
%   See also baronset sbaron
