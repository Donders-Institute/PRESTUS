%% Changelog for MINIMIZE()

%% *v9* (2016/October/25)
% * *FIXED*: Constraints in objective function were not taken into account
%   in the finalize() nested function 
%

%% *v8* (Jul 07, 2014)
% * *FIXED*: loop range issue in setoptimoptions
% * *FIXED*: simple linear constraint in one of examples was incorrectly 
%   transposed
%
% Thanks to pag (<a https://github.com/pag>) for reporting these issues.
%

%% *v7* (Mar 13, 2014)
%
% * *NEW*: massive documentation update
% * *CHANGED*: completed demo
% * *NEW*: renamed |optimize()| -> |minimize()|
% * *CHANGED*: separated |fminlbfgs| again into separate file; less
%   documentation to do that way.
%

%% *v6* (Feb 26, 2014, unpublished))
%
% * *NEW*: custom |setoptimoptions()| function, akin to 
%   |optimset()|, but suited for custom options.
% * *CHANGED*: function signature: 
%
% # |'algorithm'|  -> absorbed in |setoptimoptions()|
% # |'strictness'| -> renamed to |'AlwaysHonorConstraints'| to better mimic 
%    |fmincon()|. Also absorbed in |setoptimoptions()|.
%
% * Final tests on derivatives (for |fminlbfgs|); ready for publication
% * *CHANGED*: removed internal Nelder-Mead method, as that confuses most 
%   users. The implementation will be posted separately on the FEX
% * *CHANGED*: contact info (4 years tends to change these things...)
% * Changelog was starting to take up the majority of the function body; 
%   separated into a  file of its own
% * Major documentation update; |'doc optimze'| is supported now (XML/HTML), 
%   along with much less verbose documentation in |'help optimize'|
% * *FIXED*: countless little (non-critical) bugs 
% 
%% *v5* (Nov 2, 2010, unpublished)
%
% * Made output functions work correctly in the global routine. 
% * *FIXED*: fixed variables were inserted incorrectly in the transformation 
%   function, when using |NelderMead()|, resulting in an error.
% * If at any one iteration in |NelderMead| all function values and/or x-values 
%   are non-finite, the algorithm terminates. This is different from how 
%   |fminsearch()| handles such a case; |fminsearch| keeps on shrinking the 
%   simplex, which is pretty useless most of the time.
% * Made the overall structure more readible
% * Added option |'popsize'| for the global routine
% * *FIXED*: problem with |outputFcn|, in that |optimvalues.fval| showed the
%   PENALIZED function value, which should of course be the UNpenalized 
%   one
% * *FIXED*: problem with globalized version destroying the 
%   |ConstraintsInObjectiveFunction| field in the options structure 
%   (apparently, |options = optimset(options, 'parameter', value)|
%   removes all non-standard fields...)
%  
%% *v4* (Oct 12, 2009, unpublished)
%
% * *NEW*: |optimize()| now also supports |fminlbfgs()|, a limited-memory, 
%   BFGS quasi-newton optimizer, also available on the file exchange. 
%   This is a major change, because it requires more options and, more 
%   importantly, derivatives for the penalized constraint functions.
% * *CHANGED*: penalties are now scaled with the objective function value
% * *FIXED*: text formatting in |output.message| was all messed up for 
%   constrained problems
% * *FIXED*: |ConstrFuncCount| wasn't updated; it always resulted in 
%   |optimize()| reporting only 2 evaluations of the constraint function
% * *FIXED*: constrained global searches had a lot of problems with finalizing 
%   their solutions; these were mostly related to the proper formatting of 
%   |output|.
% * *FIXED*: another small bug; |lb| and |ub| were swapped during 
%   initialization, resulting in an error when |lb| was given but |ub| 
%   was omitted. 
% * *CHANGED*: the default algorithm back to |fminsearch|. This is of course 
%   the least buggy, and most well-known to everyone.
% * *FIXED*: corrected small mistake in the penalty function for linear 
%   equality constraints; an |abs()| was omitted, so that the sum of the 
%   constraint violation could result in lower penalties than actually 
%   deserved.
%
% 
%% *v3* (Aug 6, 2009)
%
% * *CHANGED*: removed dependency on the |optimset()| from the optimization 
%   toolbox (TolCon is not part of the diet-|optimset| most people have)
% * Included basic global optimization routine. If |x0| is omitted, but 
%   |ub| and |lb| are given, a number of initial values will be optimized, 
%   which are randomly generated in the boundaries |lb| and |ub|.
% * added one |exitflag| (-3). This is the exitflag for when the globalized 
%   algorithm has found nothing else than |inf| or |NaN| function values.
% 
% 
%% *v2*(Aug 5, 2009)
%
% * *NEW*: included |NelderMead()| algorithm, a slightly more robust and 
%   internally efficient method than |fminsearch|. Of course, |fminsearch| can 
%   still be selected through (yet another) argument.
% * one-dimensional functions could not be optimized because of a bad 
%   reference to elements of |lb| and |ub|; problem fixed. 
% * I forgot that |strictness| and |options| can also be empty; included 
%   an |isempty()| check. 
% * *FIXED*: problem with the |outputFcn|; this function was still evaluated 
%   under the assumption that |x0| is always a vector. Inserted a reshaping 
%   operation. 
% 
% 
%% *v1* (Aug 1, 2009)
%
% * |x0| can now be a matrix, just as in |fminsearch|.
% * *FIXED*: minor bug with the |strictness| setting: with multiple nonlinear 
%   constraints, an |any()| was necessary, which was not used. 
%  
%% *v0* (May 24, 2009)
% * First release.

