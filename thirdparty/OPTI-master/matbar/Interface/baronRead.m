%BARONREAD  Read a BARON .bar model into MATLAB
%
%   THIS FUNCTION WAS WRITTEN FOR INTERNAL USE ONLY AND THEREFORE DOES NOT
%   ACCEPT THE FULL BARON SYNTAX. IT IS PROVIDED AS A COURTESY ONLY, IT WILL
%   NOT READ EVERY BARON PROBLEM CORRECTLY (but hopefully most). USE AT
%   YOUR OWN RISK.
%
%   [fun,lb,ub] = baronRead(filename) reads a BARON .bar model into Matlab. 
%   filename is the name of a BARON .bar file on the Matlab path, or an 
%   absolute path to the .bar. If the extension '.bar' is  not used it 
%   will be automatically appended. Returned will be an anonymous function
%   fun with the objective, and the decision variable bounds lb and ub.
%
%   [fun,lb,ub,nlcon,cl,cu] = baronRead(filename) also returns a vector
%   anonymous function of the constraints (both nonlinear and linear),
%   together with lower and upper constraint bounds.
%
%   [fun,lb,ub,nlcon,cl,cu,xtype] = baronRead(filename) returns a character
%   string of the integrality of the decision variables.
%
%   [fun,...,xtype,x0] = baronRead(filename) returns the starting point of
%   the model.
%
%   [fun,...,x0,opts] = baronRead(filename) returns a structure with four
%   fields:
%           opts.probname:  Filename excluding extension
%           opts.sense:     Objective sense ('min' or 'max')
%           opts.eqtype:    Equation types (0, 1, 2) as per baron()
%           opts.brvarpr:   Branching variable priority as per baron()
%
%   [fun,...,opts,maps] = baronRead(filename) returns a structure with two
%   fields:
%           maps.vars:      A container.Map object with variable names and indices
%           maps.eqs:       A container.Map object with equation names and indices
%
%   baronRead(filename,writeM) will write a MATLAB mat file of the
%   model read when writeM = 1, or a MATLAB m file if writeM = 2. If the 
%   model was called ex1.bar, it will be saved as ex1.mat (or .m) in the 
%   current directory.
%
%   baronRead(filename,writeM,dispBar) will create a progress bar so the
%   user can track the reading of large files. The default is dispBar = 0.
%
%   NOTE: This function DOES NOT parse x^y^z to MATLAB format!
%
% N.B.: this function does not appear to work on WIN if an underscore is in the name of the input file
