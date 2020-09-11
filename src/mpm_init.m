% mpm init function
%==========================================================================
function [mdl, output] = mpm_init (mdl, time)

% Init active nodes
mdl.node.active = [];

% Init particle initial position
mdl.part.initial_position = mdl.part.position;

% Init node-particle maps
mdl.node.part = cell(1, mdl.node.n);
mdl.part.node = cell(1, mdl.part.n);

% Init shape function data
mdl.elem.shape_fcn = cell(1, mdl.part.n);
mdl.elem.dshape_fcn = cell(1, mdl.part.n);

% Init deformation gradient
mdl.part.deformation_grad = cell(1, mdl.part.n);
mdl.part.deformation_grad(:) = {eye(mdl.dim)};

% Init time increment
mdl.time.incr = 0;

% Init report time instants
mdl.report.curr = 0;
mdl.report.time = linspace(0, mdl.time.simulation, mdl.report.n+1);

% Update model for initial config
mdl = mpm_update (mdl, time);

% Init output data
output = struct('time', [], 'model', []);

% Report initial configuration
[mdl, output] = mpm_report (mdl, output, time);
