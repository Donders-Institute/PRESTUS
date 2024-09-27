function test
    
    clc
    
    
    
    options = setoptimoptions(...
        %{
        'algorithm', 'fminsearch',...        
        %} 
    'algorithm', 'fminlbfgs',...       
    'goalsexactachieve', 0,...
    'hessupdate', 'bfgs',...
    'Largescale', 'on',...
    'storeN', 20,...
        'display', 'iter',...        
        'TolCon', 1e-10,...        
        'TolX' , 1e-13,...
        'TolFun' ,1e-13,...
        'Gradobj', 'on',...
        'GradConstr', 'on',...
        'MaxIter', 1e2,...
        'MaxFunEvals', 5e4,...
        'popsize', 10);
   


    [sol, fval, exitflag, output] = ...
        minimize(@himmelblau, [10;-1], [-1 1],[6], [1 1],[sqrt(2)], ...
        5*[-1; -1], 5*[1; 1], @nonlcon, options);
sol, fval

([1 1]*sol)/sqrt(2)
[c, ceq] = nonlcon(sol)  

disp ' '
disp ' '
disp ' '

options = optimset(...
        %{
        'algorithm', 'fminsearch',...        
        %}  
'algorithm', 'active-set',...       
        'display', 'off',...        
        'TolCon', 1e-10,...        
        'TolX' , 1e-10,...
        'TolFun' ,1e-10,...
        'Gradobj', 'on',...
        'GradConstr', 'on',...
        'MaxIter', 1e4,...
        'MaxFunEvals', 5e4);  
        
    
     [sol2, fval2, ~, output] = ...
        fmincon(@himmelblau, [10;-1], [-1 1],[6], [1 1],[sqrt(2)], ...
        5*[-1; -1], 5*[1; 1], @nonlcon, options)
%     
output

% ([1 1]*sol2)/sqrt(2)
   
  [c, ceq] = nonlcon(sol2)  
   


   
    
%     [sol, fval, exitflag, output] = ...
%         minimize(@banana, [1; 1], [1 1],[3], [1 1],[2], ...
%         5*[-1 -1], 5*[1 1], @nonlcon, options)
    
    
end

function [obj, grad_obj] =  banana(x)
  
    
    obj = (1-x(1)).^2 + 100*(x(2)-x(1).^2).^2;
    grad_obj = [
        -2*(1-x(1)) - 400*(x(2)-x(1).^2).*x(1)
        200*(x(2)-x(1).^2)
    ]';
end


function [obj, grad_obj] =  himmelblau(x)

obj = (x(1).^2 + x(2) - 11).^2 + (x(1) + x(2).^2 - 7).^2;
grad_obj = [
    4*x(1)*(x(1).^2 + x(2) - 11) + 2*(x(1) + x(2).^2 - 7)
    2*(x(1).^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2).^2 - 7)
    ];

end

function [c, ceq, grad_c, grad_ceq] = nonlcon(x)
    ceq   = x(:).'*x(:) - 1;
    c = [];
    
    grad_c = [];
    grad_ceq = 2*x;
end





