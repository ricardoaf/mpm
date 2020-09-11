% mpm report
%==========================================================================
function [mdl, output] = mpm_report (mdl, output, time)

% Check if report step
if time >= mdl.report.time(mdl.report.curr+1)
    
    % Update current report
    mdl.report.curr = mdl.report.curr + 1;
    
    % Print progress
    if mdl.report.show_progress
        fprintf('time: %8.4f  [%7.3f%%]\n', ...
            time, time/mdl.time.simulation*100);
    end
    
    % Plot particle kinematics
    if mdl.report.show_animation
        mpm_view(mdl);
    end
    
    % Increment output
    output.time = [output.time time+mdl.time.incr];
    output.model = [output.model mdl];
end
