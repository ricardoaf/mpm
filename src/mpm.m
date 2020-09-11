% mpm core function
%==========================================================================
function output = mpm (input)

% create model
mdl = input;

% init model and output
time = 0;
[mdl, output] = mpm_init (mdl, time);

% Loop on time
tic;
while time < mdl.time.simulation
    
    % update model
    mdl = mpm_update (mdl, time);
    
    % solve mpm explicit step
    mdl = mpm_step (mdl);
    
    % update time
    time = time + mdl.time.incr;
    
    % check for report
    [mdl, output] = mpm_report (mdl, output, time);
    
end
toc