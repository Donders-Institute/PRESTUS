[![View minimize on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/24298-minimize)

[![Donate to Rody](https://i.stack.imgur.com/bneea.png)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url)

# FEX-minimize

( NOTE: adding the main folder and its sub-folders to the MATLAB search path will enable you to view the extended documentation in the MATLAB help browser. )
MINIMIZE is an improvement upon the functions FMINSEARCHBND and FMINSEARCHCON written by John d'Errico (also available on the file exchange). It solves the optimization problem
min f(x)

s.t. 

lb <= x <= ub
A * x < b
Aeq * x = beq
c(x) <= 0
ceq(x) = 0

using a coordinate transformation for the bound constraints, and penalty functions for the other constraints. The penalty functions used are pseudo-adaptive, in that they are designed to penalize heavily yet prevent overflow from ever happening.

The main differences between MINIMIZE and FMINSEARCHCON are

 - rudimentary support for global optimization problems
 - it handles (non)linear equality constraints 
 - strictness is more controllable 
 - support for FMINLBFGS 

While FMINSEARCHCON does not permit ANY function evaluation outside the feasible domain, MINIMIZE can be either allowed (default) or disallowed ('AlwaysHonorConstraints' option) to do so. 

Its behavior is similar to that of FMINCON (optimization toolbox), which makes it useful for those who do not have the optimization toolbox, but only have one-off simple problems to solve, or are considering buying the toolbox and want to practice a bit with FMINCON's interface. 

Note that MINIMIZE is by no means intended to be a full replacement of FMINCON, since the algorithms used by FMINCON are simply better. However, it does have a few notable advantages.

For relatively small problems where it is hard or impossible to come up with good initial estimates, MINIMIZE can be used effectively as a global optimization routine, providing a simple means to find good initial estimates. 

It is also particularly useful in cases where the objective function is costly to compute, and hard or impossible to differentiate analytically. In such cases, FMINCON is forced to compute the derivatives numerically, which usually takes > 60% of the computation time if you have a sizeable problem. Since FMINSEARCH is the engine for MINIMIZE , no derivatives are required, which might make it more efficient than using FMINCON. 

With the addition of FMINLBFGS (included in this publication), it is also useful for extremely large problems (over 3000 variables). Hessian information, even with BFGS solvers, can consume large amounts of memory. This problem is solved by using a Limited-memory version of the BFGS routines, implemented nicely in FMINLBFGS.


Usage:

sol = MINIMIZE(func, x0) 
sol = MINIMIZE(func, x0, lb,ub)
sol = MINIMIZE(func, x0, A,b) 
sol = MINIMIZE(func, x0, A,b, Aeq,beq) 
sol = MINIMIZE(func, x0, A,b, Aeq,beq, lb,ub) 
sol = MINIMIZE(func, x0, A,b, Aeq,beq, lb,ub, nonlcon, options) 

[sol, fval] = OPTIMIZE(func, ...)
[sol, fval, exitflag] = MINIMIZE(func, ...)
[sol, fval, exitflag, output] = MINIMIZE(func, ...)


( NOTE: adding the main folder and its sub-folders to the MATLAB search path will enable you to view the extended documentation in the MATLAB help browser. )


If you find this work useful, please consider [a donation](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url).
