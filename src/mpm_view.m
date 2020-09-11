% mpm view
%==========================================================================
function mpm_view (mdl)

% 2D only **
if mdl.dim~=2, return; end

global MPM_GRID_PLOT MPM_PARTICLE_PLOT
if mdl.report.curr==1
    
    figure('Color', [1 1 1]); hold on; axis equal
    x = mdl.node.coord(:,1); y = mdl.node.coord(:,2);
    MPM_GRID_PLOT = patch(x(mdl.elem.node)', y(mdl.elem.node)', 'white');
    MPM_PARTICLE_PLOT = plot(mdl.part.position(:,1), ...
        mdl.part.position(:,2), 'ro');
    hold off; drawnow;
    
else
    set(MPM_PARTICLE_PLOT, 'XData', mdl.part.position(:,1));
    set(MPM_PARTICLE_PLOT, 'YData', mdl.part.position(:,2));
    drawnow;

end
