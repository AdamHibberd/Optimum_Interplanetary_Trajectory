%SBARON Solve a NLP/MINLP using the Global MINLP Solver BARON (Alternative Interface)
%
%   [x,fval,exitflag,info,allsol = sbaron('field1',value1,'field2',value2,...)
%
%   sbaron() provides an alternative interface to the BARON solver where
%   parameters can be specified using a [field,value] pairing. For a
%   detailed description of each field, see the function baron.m.
%
%   Fields may be:
%
%   Field              Value                             Data Type
%   'fun'       Objective Function                  Function Handle
%   'lb'        Lower Bounds                        Double Column Vector
%   'ub'        Upper Bounds                        Double Column Vector
%   'nlcon'     Nonlinear Constraints Function      Function Handle
%   'cl'        Nonlinear Constraints Lowed Bound   Double Column Vector
%   'cu'        Nonlinear Constraints Upper Bound   Double Column Vector
%   'A'         Linear Constraints A Matrix         Double Matrix
%   'rl'        Linear Constraints Lower Bound      Double Column Vector
%   'ru'        Linear Constraints Upper Bound      Double Column Vector
%   'xtype'     Variable Integrality                Character Row Vector
%   'x0'        Initial Solution Guess              Double Matrix / Vector
%   'opts'      Solver Options (from baronset)      Structure
%
%   
%   Example usage of this interface is shown below:
%
%           [x,fval] = sbaron('fun',objFun,'lb',lb,'ub',ub)
%
%
%   See also baron baronset
