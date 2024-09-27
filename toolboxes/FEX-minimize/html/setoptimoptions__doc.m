%% SETOPTIMOPTIONS                   
%
% Optimization options for minimize

%% Usage
%
% |options = setoptimoptions('param1',value1, 'param2',value2, ...)|
%
% The relevant parameter value pairs are defined below.

%% Used options
%
% Options used by |minimize()| are: 
%
% *Same as in <matlab:doc('optimset') optimset>*: 
%
% |TolFun|, |TolX|, |OutputFcn|, |PlotFcn|, |MaxIter|, |MaxFunEvals|, |Display|.

%%
% *Specific to <matlab:doc('minimize') minimize>*: 
%
%       - 'TolCon'       
%           Tolerance used on any constraint. This means that minimize()
%           considers the constraints are only violated when they exceed
%           this amount of violation; otherwise, the constraints are
%           considered met. Defaults to 1e-8.
%           
%       - 'GradObj'     
%           Specifies whether the objcetive function returns gradient
%           information as its second output argument. Valid values are
%           'on' and 'off' (the default). In case this option is set to 
%           'off', gradient information is computed via finite 
%           differences.
%     
%       - 'GradConstr'  
%           Specifies whether the non-linear constraint function returns 
%           Jacobian information as its third and fourth output arguments. 
%           Valid values are 'on' and 'off' (the default). In case this
%           option is 'off', Jacobian information is computed via finite
%           differences.
%
%       - 'FinDiffType'  
%           Type of finite differences to use. Valid values are 'forward'
%           (the default), 'backward', and 'central'. Central differences
%           provide the best accuracy, but require twice as many function
%           evaluations.
%
%       - 'DiffMaxChange' 
%           Maximum change in the objective/constraint function to use when
%           computing gradient/Jacobian information with finite
%           differences. Defaults to 1e-1.
%
%       - 'DiffMinChange'  
%           Minimum change in the objective/constraint function to use when
%           computing gradient/Jacobian information with finite
%           differences. Defaults to 1e-8.
%
%       - 'AlwaysHonorConstraints'
%           By default, minimize() will assume the objective (and  
%           constraint) function(s) can be evaluated at ANY point in 
%           RN-space; the initial estimate does not have to lie in the 
%           feasible region, and intermediate solutions are also allowed 
%           to step outside this area. this is equal to setting this option 
%           to 'none'. 
%
%           If the objective function cannot be evaluated outside the
%           feasible region, set this argument to 'bounds' (bound 
%           constraints will never be broken) or 'all' (also the linear 
%           constraints will never be broken). Note that the non-linear 
%           constraints will remain satisfied within options.TolCon. 
%       
%           When using 'Bounds' or 'All', the initial estimate [x0]
%           MUST be feasible. If it is not feasible, an error is produced
%           before the objective function is ever evaluated. 
%
%       - 'Algorithm'
%           By default, this is set to MATLAB's own derivative-free 
%           Nelder-Mead algorithm, implemented in fminsearch. minimize() 
%           supports another algotithm, fminlbfgs: a limited-memory,
%           Broyden/Fletcher/Goldfarb/Shanno optimizer. Use this algorithm
%           when your objective function has many free variables, e.g., [x]
%           is large.
%
%       - 'popsize'
%           Used by the global optimization routine. This is the number of
%           randomized initial values to use, and thus the number of times
%           to repeat the call to minimize(). Defaults to 20× the number of
%           elements in [x0].

%%
% Specific to <matlab:doc('fminlbfgs') fminlbfgs>:
%
%       - 'GoalsExactAchieve'
%           If set to 0, a line search method is used which uses a few
%           function calls to do a good line search. When set to 1 a normal
%           line search method with Wolfe conditions is used (default).
%
%       - 'HessUpdate'
%           If set to 'bfgs', Broyden/Fletcher/Goldfarb/Shanno
%           optimization is used (default), when the number of unknowns is 
%           larger than 3000 the function will switch to Limited-memory BFGS,
%           or if you set it to 'lbfgs'. When set to 'steepdesc', steepest
%           decent optimization is used.
%
%       - 'StoreN' 
%           Number of itterations used to approximate the Hessian,
%           in L-BFGS, 20 is default. A lower value may work better with
%           non smooth functions, because than the Hessian is only valid for
%           a specific position. A higher value is recommend with quadratic 
%           equations.
%
%       - 'rho'   
%           Wolfe condition on gradient (c1 on wikipedia), default 0.01.
%
%       - 'sigma'  
%           Wolfe condition on gradient (c2 on wikipedia), default 0.9.
%
%       - 'tau1' 
%           Bracket expansion if stepsize becomes larger, default 3.
%
%       - 'tau2'
%           Left bracket reduction used in section phase, default 0.1.
%
%       - 'tau3'   
%           Right bracket reduction used in section phase, default 0.5.

%% See also 
%
% <matlab:doc('optimset') optimset>, <matlab:doc('optimget') optimget>.
