function BUGS 
    % BUG: strict does not work; reported on FEX

clc


Aeq = ones(1,3); 
beq = 1; 

options = optimset;
options.Algorithm = 'fminlbfgs';
%options.HessUpdate = 'bfgs';
options.Display = 'off';
options.TolX = 1e-12;
options.TolFun = 1e-12;
options.MaxFunEvals = 1e4;
options.MaxIter = 1500;

options.AlwaysHonorConstraints = 'Bounds';

% FIXME: ignore nonlcon() if this is true?
options.ConstraintsInObjectiveFunction = true;

a0 = [0.3 0.4 0.3]'; 
%a0 = [];

lb = [0 0 0]';
ub = [1 1 1]';

[sol0,fval0,exit0,output0] = optimize(@L1,a0, [],[], Aeq,beq, lb,ub, [])
[sol1,fval1,exit1,output1] = optimize(@L1,a0, [],[], Aeq,beq, lb,ub, [],  options)
output1.message

% FIXME: fminlbfgs will consume all function evaluations without doing
% anything
%options.AlwaysHonorConstraints = 'All';
[sol2,fval2,exit2,output2] = optimize(@L1,a0, [],[], Aeq,beq, lb,ub, @nonlcon,  options)
output2.message

% just to check 
options.Algorithm = 'active-set';
[sol3,fval3,exit3,output3] = fmincon(@L1,a0, [],[], Aeq,beq, lb,ub, @nonlcon,  options)
output3.message

end

function [y, c, ceq] =  L1(x)
    y = sum(x.^x);
    c = [];
    ceq = sum(x(:))-1;
end

function [c,ceq] = nonlcon(x)
    c = [];
    ceq = sum(x(:))-1;
end