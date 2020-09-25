

%
% Part performing post-simulation analysis when needed
%

[solution_post] = runPostSolveTasks(problem,solution,options,data);    % Output solutions
options.plot=4;
genSolutionPlots(options,solution_post)