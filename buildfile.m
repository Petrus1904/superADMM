function plan = buildfile
import matlab.buildtool.tasks.MexTask

% Create a plan with no tasks
plan = buildplan;

% Add a task to build a MEX file
plan("mex") = MexTask("explore.c","output",Options=["-R2018a" "-silent"]);