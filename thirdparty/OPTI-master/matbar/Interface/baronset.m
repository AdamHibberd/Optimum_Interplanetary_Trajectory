%BARONSET  Create or alter the options for optimization with BARON
%
% options = baronset('param1',value1,'param2',value2,...) creates a BARON
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the BARON
% default.
%
% options = baronset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (writes over) the parameters
% specified by 'param' and 'value'.
%
% options = baronset() creates an options structure with all fields set to
% BARON defaults.
%
% baronset() prints a list of all possible fields and their function.
