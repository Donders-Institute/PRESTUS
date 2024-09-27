%% Constrained function minimization with |fminsearch| and |fminlbfgs|

function minimize_demo
      
    %% Unconstrained optimization    
    %
    % first, define a test function:
    clc, rosen = @(x) ( (1-x(1))^2 + 105*(x(2)-x(1)^2)^2 ) /1e4;

    %%
    % This is the classical Rosenbrock function, which has a global minimum
    % at $f(x) = f([1, 1]) = 0$. The function is relatively hard to minimize,
    % because that minimum is located in a long narrow ``valley'':
    k = 0; range = -5:0.1:5;
    z = zeros(numel(range));    
    for ii = range
        m = 0; k = k + 1;        
        for jj = range
            m = m + 1;
            z(k, m) = rosen([ii, jj]);            
        end
    end  
    [y, x] = meshgrid(range, range);
    
    S = surf(x, y, z, 'linestyle', 'none');
    view(-213, 38)
    axis tight    
            
    shading interp
    material metal
    lighting gouraud
    colormap('hot')
    
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.8)     
        
    %%
    % Optimizing the fully unconstrained problem with |minimize| indeed 
    % finds the global minimum:
    solution = minimize(rosen, [3 3])    
        
    
    %% Optimization with bound constraints
    %
    % Imposing a lower bound on the variables gives 
    [solution, fval] = minimize(rosen, [3 3], [],[], [],[], [2 2])
    
    %%
    % in the figure, this looks like
    zz = z;   zz(x > 2 & y > 2) = NaN;  
    ZZ = z;   ZZ(x < 2 & y < 2) = NaN;
    
    figure('renderer', 'opengl');
    hold on
    
    S(1) = surf(x, y, zz,... 
                'linestyle', 'none',...
                'FaceAlpha', 0.2);
            
    S(2) = surf(x, y, ZZ,...
                'linestyle', 'none');
    
    plot3(solution(1), solution(2), fval+0.5, ...
          'gx', ...
          'MarkerSize', 20,...
          'linewidth' , 5)    

    xlabel('X(1)')
    ylabel('X(2)')  
    
    view(-196, 38)
    grid on
    axis tight      
    
    shading interp
    material metal
    lighting gouraud
    colormap('hot')
    
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.8); 
        
    %%
    % Similarly, imposing an upper bound yields
    solution = minimize(rosen, [3 3], [],[], [],[], [],[0.5 0.5])
    
    zz = z;   zz(x < 0.5 & y < 0.5) = NaN;  
    ZZ = z;   ZZ(x > 0.5 & y > 0.5) = NaN;
    
    figure('renderer', 'opengl');
    hold on
    
    S(1) = surf(x, y, zz, ...
                'linestyle', 'none',...
                'FaceAlpha', 0.2);
            
    S(2) = surf(x, y, ZZ, ...
                'linestyle', 'none'); 
    
    plot3(solution(1), solution(2), fval+0.5, 'gx', ...
          'MarkerSize', 20,...
          'LineWidth', 5);        
    
    xlabel('X(1)')
    ylabel('X(2)')    
    view(201, 38)
    grid on
    axis tight      
    
    shading interp
    material metal
    lighting gouraud
    colormap('hot')
    
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.8); 
 
    %%
    % Minimize with $x_2$ fixed at 3. In this case, |minimize| simply
    % removes the variable before |fminsearch| sees it, essentially 
    % reducing the dimensionality of the problem. This is particularly
    % useful when the number of dimensions _N_ becomes large.
    minimize(rosen, [3 3], [],[], [],[], [-inf 3], [inf 3])
    
    
    %% Linear constraints
    %
    % You can use linear inequality or equality constraints. For
    % example, with the constraints
    %
    % $A*x \leq b$,              
    % $A_{eq}*x == b_{eq}$,    
    %
    % with A = [+2 +1], b = -2 and Aeq = [+1 -1], beq = -2. 
    %
    % |minimize()| finds the following result:    
    [solution, fval] = minimize(rosen, [3;3], [2 1],-2, [1 -1],-2)
        
    %% 
    % These constraints look like the following:
    xinds = 2*x+y <= -2;
    zz = z;   zz( xinds ) = inf;  
    ZZ = z;   ZZ(~xinds ) = inf;
    
    Ax = z;   Axinds = abs(2*x+y + 2) < 1e-3;        
    [x1, sortinds] = sort(x(Axinds)); 
    Ax = Ax(Axinds); Ax = Ax(sortinds);
    y1 = y(Axinds); y1 = y1(sortinds);
    
    Aeq = z;   Aeqinds = abs(x-y + 2) < 1e-3;        
    [x2, sortinds] = sort(x(Aeqinds)); 
    Aeq = Aeq(Aeqinds); Aeq = Aeq(sortinds);
    y2 = y(Aeqinds); y2 = y2(sortinds);
    
    figure('renderer', 'opengl');
    hold on
    
    l1 = line([x1(1:end-1)';x1(2:end)'],...
              [y1(1:end-1)';y1(2:end)'],...
              [Ax(1:end-1)';Ax(2:end)']);    
    l2 = line([x2(1:end-1)';x2(2:end)'],...
              [y2(1:end-1).';y2(2:end)'],...
              [Aeq(1:end-1).';Aeq(2:end)']);    
    
    S(1) = surf(x, y, zz, 'linestyle', 'none', 'FaceAlpha', 0.2);
    S(2) = surf(x, y, ZZ, 'linestyle', 'none');      
    
    l3 = plot3(solution(1)+0.4, solution(2)+0.8, fval+0.5, 'gx',...
        'MarkerSize', 20,...
        'LineWidth', 5);
    
    set(l1, 'color', 'b', 'linewidth', 2);
    set(l2, 'color', 'k', 'linewidth', 2);         
    
    view(150, 30)
    grid on
    axis tight      
    
    xlabel('X(1)', 'interpreter', 'LaTeX');
    ylabel('X(2)', 'interpreter', 'LaTeX');
    
    k = legend([l1(1); l2(1); l3],'inequality $$A\mathbf{x} \leq -2$$', ...
        'equality $$A_{eq}\mathbf{x} = -2$$', 'Solution');
    
    set(k, 'interpreter', 'LaTeX', 'location', 'NorthWest');   
    
    shading interp
    material metal
    lighting phong
    colormap('autumn')
    
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.8);
 
    %% Non-linear constraints
    %
    % Also general nonlinear constraints can be used. A simple example:
    %
    % nonlinear inequality: 
    %
    % $$\sqrt{x_1^2 + x_2^2} \leq 2$$
    %
    % nonlinear equality  : 
    %
    % $$0.2x_1^2 + 0.4x_2^3 = 1$$
        
    options = setoptimoptions(...
        'TolFun', 1e-6, ...
        'TolX'  , 1e-6, ...
        'MaxFunEvals', inf,...
        'MaxIter', 1e4);
    
    [sol, fval,...
     ~,...
     output] = minimize(rosen,...
                        [-3; 3],...
                        [],[],...
                        [],[],...
                        [],[],...
                        @nonlcon,...
                        options);
    
    %%
    % Note that |nonlcon| is a subfunction, listed below. 
    %
    % These constraints look like the following:
    zz = z;   zz(sqrt(x.^2 + y.^2) <= 2)   = inf;  
    ZZ = z;   ZZ(sqrt(x.^2 + y.^2) >= 2.2) = inf;
    zZ = z;   zZ(x.^2 + y.^3 >= 1.0 + 0.1) = inf; 
              zZ(x.^2 + y.^3 <= 1.0 - 0.1) = inf;
      
    xX = x(isfinite(zZ));  [xX, inds] = sort(xX);
    yY = y(isfinite(zZ));  yY = yY(inds);
    zZ = zZ(isfinite(zZ)); zZ = zZ(inds);
              
    figure('renderer', 'opengl');
    hold on
    
    S(1) = surf(x, y, zz, 'linestyle', 'none', 'FaceAlpha', 0.2);
    S(2) = surf(x, y, ZZ, 'linestyle', 'none');
    
    L = line([xX(1:end-1)';xX(2:end)'],[yY(1:end-1)';yY(2:end)'],[zZ(1:end-1)';zZ(2:end)']);        
    l3 = plot3(sol(1)+0.4, sol(2)+0.5, fval+1, 'gx', 'MarkerSize', 20, 'linewidth', 5);
    
    set(L, 'linewidth', 2, 'color', 'b');
    view(150, 50)
    grid on
    axis tight      
    
    k = legend([S(2); L(1); l3],'non-linear inequality $$c(x) < 0$$', ...
        'non-linear equality $$c_{eq}(x) = 0$$', 'Solution');
    set(k, 'interpreter', 'LaTeX', 'location', 'NorthWest');   
    
    shading interp
    material metal
    lighting phong
    colormap('autumn')
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.8); 
    
   
    %%
    % Note that the |output| structure contains a field |constrviolation|:
    output
    
    %%
    % The contents of which shows that all constraints have been satisfied:
    output.constrviolation
    output.constrviolation.nonlin_eq{:}
    output.constrviolation.nonlin_ineq{:}

    %% Global optimization
    %
    % This is the 2D sine-envelope-sine function. It has a single global
    % minimum at [0,0], where the function assumes a value of 0. As you can
    % imagine, it is hard to find this minimum when the initial estimates
    % is not very close to the minimum already:
    
    sinenvsin = @(x) 3*sum( (sin(sqrt(x(:).'*x(:))).^2 - 0.5)./(1 + 0.001*x(:).'*x(:)).^2 + 0.5, 1);
    
    figure('renderer', 'opengl');
    hold on
    
    k = 0; range = -10:0.1:10;
    z = zeros(numel(range));    
    for ii = range
        m = 0; k = k + 1;        
        for jj = range
            m = m + 1;
            z(k,m) = sinenvsin([ii jj]);            
        end
    end  
    [y, x] = meshgrid(range, range);
    
    S = surf(x, y, z, 'linestyle', 'none');
    
    axis equal
    view(-148,24)    
    
    shading interp
    material shiny
    lighting phong 
    colormap('autumn')
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.6); 
    
    %%
    % |minimize()| provides rudimentary support for this type of problem.
    % Omitting the initial value x0 will re-start |minimize()| several times 
    % at randomly chosen initial values in the interval [lb ub]: 
    
    options = setoptimoptions(...
        'popsize', 1e2, 'maxfunevals', 1e4,  'maxiter', 1e2);
    
    [sol,fval] = minimize(sinenvsin, ...
                          [],...
                          [],[],...
                          [],[],...
                          -[5 5], +[5 5],...
                          [],...
                          options)
    
    %%
    % Naturally, these types of problems may also have constraints:
    
    [solution,fval] = minimize(sinenvsin, [], [2 1],-2, [1 -1],-2,...
        -[5;5], +[5;5], [], options);    
    
    xinds = 2*x+y <= -2;
    zz = z;   zz( xinds ) = inf;  
    ZZ = z;   ZZ(~xinds ) = inf;
    
    Ax = z;   Axinds = abs(2*x+y + 2) < 1e-3;        
    [x1, sortinds] = sort(x(Axinds)); 
    Ax = Ax(Axinds); Ax = Ax(sortinds);
    y1 = y(Axinds); y1 = y1(sortinds);
    
    Aeq = z;   Aeqinds = abs(x-y + 2) < 1e-3;        
    [x2, sortinds] = sort(x(Aeqinds)); 
    Aeq = Aeq(Aeqinds); Aeq = Aeq(sortinds);
    y2 = y(Aeqinds); y2 = y2(sortinds);
    
    figure('renderer', 'opengl');
    hold on
    
    S(1) = surf(x, y, zz, 'linestyle', 'none', 'FaceAlpha', 0.2);
    S(2) = surf(x, y, ZZ, 'linestyle', 'none');
    
    l1 = line([x1(1:end-1)';x1(2:end)'],...
              [y1(1:end-1)';y1(2:end)'],...
              [Ax(1:end-1)';Ax(2:end)']);    
    l2 = line([x2(1:end-1)';x2(2:end)'],...
              [y2(1:end-1).';y2(2:end)'],...
              [Aeq(1:end-1).';Aeq(2:end)']);    
    l3 = plot3(solution(1)+0.5, solution(2)+2.8, fval+2,...
               'gx', ...
               'MarkerSize', 20,...
               'linewidth', 5);
           
    xlabel('X(1)', 'interpreter', 'LaTeX'); 
    ylabel('X(2)', 'interpreter', 'LaTeX');
    
    set(l1, 'color', 'r', 'linewidth', 2);
    set(l2, 'color', 'k', 'linewidth', 2);     
    
    view(150, 30)
    grid on
    axis tight          
    
    k = legend([l1(1); l2(1); l3],'inequality $$A\mathbf{x} \leq -2$$', ...
        'equality $$A_{eq}\mathbf{x} = -2$$', 'Solution');    
    set(k, 'interpreter', 'LaTeX', 'location', 'NorthWest'); 
    
    view(170,80);  
    
    shading interp
    material shiny
    lighting phong
    colormap('autumn')
    light('style', 'local', 'position', [-3 0 5]);
    set(S, 'ambientstrength', 0.6); 
    
    %% Different algorithm: |fminlbfgs|
    %
    % |minimize()| also supports |fminlbfgs|, a limited-memory, 
    % Broyden/Fletcher/Goldfarb/Shanno optimizer, implemented by Dirk-Jan
    % Kroon. Unlike |fminsearch()|, this algorithm uses gradient
    % information to improve the overall optimization speed: 
    
    options = setoptimoptions(...
        'Algorithm', 'fminsearch');    
    [solution, fval, ~, output] = minimize(...
        rosen, [-3 -18], [],[], [],[], [],[], [], options);  
    
    solution, fval
    output.funcCount
    
    options = setoptimoptions(...
        'Algorithm', 'fminlbfgs');    
    [solution, fval, ~, output] = minimize(...
        rosen, [-3 -18], [],[], [],[], [],[], [], options);    
    
    solution, fval
    output.funcCount
    
    %%
    % As can be seen, |fminlbfgs()| indeed uses less funtion evaluations to
    % come to a comparable answer. 
    %
    % The great advantage of this minimizer over |fminsearch()| is that
    % |fminlbfgs()| can handle very large problems efficiently. For
    % problems of higher dimensionality, |fminsearch()| has poor 
    % convergence properties compared to |fminlbfgs()|. 
   
    %% Supplying gradients to |fminlbfgs|
    %
    % |fminlbfgs()| needs gradient information of both the objective and
    % non-linear constraint functions. |minimize()| estimate this
    % information numerically via finite differences, but this can be
    % costly for larger problems. Therefore, |minimize()| can also accept
    % gradient information computed by the objective funcion:
    
    options = setoptimoptions(...
        'TolX', 1e-8,...
        'TolFun', 1e-8,...
        'FinDiffType', 'central',...
        'MaxFunEvals', inf,...
        'MaxIter', 1e4,...
        'GoalsExactAchieve', 0,...
        'Algorithm', 'fminlbfgs',...  % This option specifies that your
        'GradObj'  , 'on');           % objective function also provides
                                      % gradient information as its second
                                      % output argument
    
    [solution, fval, ~, output] = minimize(...
        @rosen_with_gradient, [-3 -3], [],[], [],[], [],[], @nonlcon, options);
    
    solution, fval
    output.ObjfuncCount
    output.ConstrfuncCount

    %%
    % In case of non-linearly constrained problems, also Jacobian
    % information of the non-linear constraint function can be provided:
        
    options = setoptimoptions(...
        'TolX', 1e-10,...
        'TolFun', 1e-10,...  
        'MaxFunEvals', inf,...
        'MaxIter', 1e4,...
        'GoalsExactAchieve', 0,...
        'Algorithm' , 'fminlbfgs',... % This option specifies that your 
        'GradObj'   , 'on',...        % non-linear constraint function also 
        'GradConstr', 'on');          % provides Jacobian information as its 
                                      % third and fourth output arguments
    
    [solution, fval, ~, output] = minimize(...
        @rosen_with_gradient, [-3 -3], [],[], [],[], [],[], @nonlcon_with_Jacobian, options);
    
    solution, fval
    output.ObjfuncCount
    output.ConstrfuncCount
            
    
    %% See also
    %
    % <matlab:doc('setoptimoptions') setoptimoptions>,
    % <matlab:doc('fminlbfgs') fminlbfgs>.

    
end

%% 

function [fVal, gVal] = rosen_with_gradient(x)    
    fVal = ( (1-x(1))^2 + 105*(x(2)-x(1)^2)^2 ) /1e4;
    gVal = [ -2*(1-x(1)) - 4*105*x(1)*(x(2)-x(1)^2)
             2*105*(x(2)-x(1)^2)
             ]/1e4;
end

function [c, ceq] = nonlcon(x)
    c   = x(1)^2 + x(2)^2 - 2;
    ceq = 0.2*x(1)^2 + 0.4*x(2)^3 - 1;
end

function [c, ceq, c_Jac, ceq_Jac] = nonlcon_with_Jacobian(x)
    c   = x(1)^2 + x(2)^2 - 2;
    ceq = 0.2*x(1)^2 + 0.4*x(2)^3 - 1;
    
    c_Jac   = 2*x;
    ceq_Jac = [0.4*x(1); 1.2*x(2)^2];
end

