function fminlbfgs__doc
    
%% |fminlbfgs|
%
% Finds a local minimum of a function of several variables.

%% Usage
%
%  [x,fval,exitflag,output,grad] = fminlbfgs(fun,x0,options)

%% Description
%
% This optimizer is developed for image registration methods with large
% amounts of unknown variables.
%
% Optimization methods supported:
% * Quasi Newton Broyden–Fletcher–Goldfarb–Shanno (BFGS)
% * Limited memory BFGS (L-BFGS)
% * Steepest Gradient Descent optimization.

%% Input arguments
%
% *|fun|* 
%
% Function handle or string which is minimized, returning an error 
% value and  optionally the error gradient.
%
% Note that the speed of this optimizer can be improved by also 
% providing the gradient at X. Write the |fun| function as follows

function [f, g] = FUN(X)
    f % = ...          function value calculation at X
    if (nargout > 1)
        g % = ...      gradient calculation at X
    end
end

%%
% *|x0|* 
%
% Initial values of unknowns can be a scalar, vector or matrix (optional)

%%
% *|options|*
%
% Structure with optimizer options, made by a struct or |optimset|. 
% Note that |optimset| does not support all input options.

%% Output arguments
%
% *|x|*
%
% The found location (values) which minimize the function.

%%
% *|fval|*
%
% The function value at the solution found.

%%
% *|exitflag|*
%
% Gives value, which explain why the minimizer stopped. Possible values 
% of |exitflag|, and the corresponding exit conditions are
%
%  		1: Change in the objective function value was less than the 
%          specified tolerance TolFun.
%  		2: Change in x was smaller than the specified tolerance TolX.
%  		3: Magnitude of gradient smaller than the specified tolerance.
%  		4: Boundary fminimum reached.
%  		0: Number of iterations exceeded options.MaxIter or number of 
%          function evaluations exceeded options.FunEvals.
%      -1: Algorithm was terminated by the output function.
%  	   -2: Line search cannot find an acceptable point along the current 
%          search.

%%
% *|output|*
%
% Structure with all important ouput values and parameters

%% 
% *|grad|*
%
% the gradient at this location.
%
%% Options
%
% *|'GoalsExactAchieve'|*
%
% If set to 0, a line search method is used which uses a few function 
% calls to do a good line search. When set to 1 a normal line search 
% method with Wolfe conditions is used (default).

%%
% *|'GradConstr'|*
%
% Set this variable to true if gradient calls are cpu-expensive 
% (default). If false more gradient calls are used and less function 
% calls.

%%
% *|'HessUpdate'|*
%
% If set to 'bfgs', Broyden–Fletcher–Goldfarb–Shanno optimization is 
% used (default), when the number of unknowns is larger then 3000 the 
% function will switch to Limited memory BFGS, or if you set it to 
% 'lbfgs'. When set to 'steepdesc', steepest decent optimization is 
% used.

%%
% *|'StoreN'|* 
%
% Number of itterations used to approximate the Hessian, in L-BFGS, 
% 20 is default. A lower value may work better with non smooth functions, 
% because than the Hessian is only valid for a specific position. A 
% higher value is recommend with quadratic equations.

%%
% *|'GradObj'|*
%
% Set to 'on' if gradient available otherwise finited difference is used.

%%
% *|'Display'|*
%
% Level of display. 'off' displays no output; 'plot' displays all linesearch 
% results in figures. 'iter' displays output at each iteration; 'final' 
% displays just the final output; 'notify' displays output only if the 
% function does not converge.

%%
% *|'TolX'|* 
%
% Termination tolerance on x, default 1e-6.

%%
% *|'TolFun'|* 
%
% Termination tolerance on the function value, default 1e-6.

%%
% *|'MaxIter'|* 
%
% Maximum number of iterations allowed, default 400.

%%
% *|'MaxFunEvals'|* 
%
% Maximum number of function evaluations allowed, default 100 times the 
% amount of unknowns.

%%
% *|'DiffMaxChange'|* 
%
% Maximum stepsize used for finite difference gradients.

%%
% *|'DiffMinChange'|*
%
% Minimum stepsize used for finite difference gradients.

%%
% *|'OutputFcn'|* 
%
% User-defined function that an optimization function calls at each 
% iteration.

%%
% *|'rho'|* 
%
% Wolfe condition on gradient (c1 on wikipedia), default 0.01.

%%
% *|'sigma'|* 
% 
% Wolfe condition on gradient (c2 on wikipedia), default 0.9.

%%
% *|'tau1'|* 
%
% Bracket expansion if stepsize becomes larger, default 3.

%%
% *|'tau2'|*
%
% Left bracket reduction used in section phase,	default 0.1.

%%
% *|'tau3'|*
%
% Right bracket reduction used in section phase, default 0.5.


%% Examples

options = optimset('GradObj','on');
X = fminlbfgs(@myfun,2,options)

%%
% where myfun is a MATLAB function such as:
function [f,g] = myfun(x)
    f = sin(x) + 3;
    if ( nargout > 1 )
        g = cos(x); end
end

%% See also 
%
% <matlab:doc('optimset') optimset>, <matlab:doc('fminsearch') fminsearch>,
% <matlab:doc('fminbnd') fminbnd>, <matlab:doc('fmincon') fmincon>, <matlab:doc('fminunc') fminunc>, <matlab:doc('@') @>, <matlab:doc('inline') inline>.

end
