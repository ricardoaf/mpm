% mpm view
%==========================================================================
function mpm_view (mdl)
global MPM_PARTICLE_PLOT

if mdl.dim==1
    
    if mdl.report.curr==1
        figure('Color', [1 1 1]); hold on; axis equal
        x = mdl.node.coord(:,1); x = [x x]; y = x;
        h = mdl.grid.elem_size; y(:,1) = -h/2; y(:,2) = h/2;
        patch(x', y', 'white');
        MPM_PARTICLE_PLOT = plot(mdl.part.position(:,1), ...
            mdl.part.position(:,1)*0, 'ro');
        hold off; drawnow;
    else
        set(MPM_PARTICLE_PLOT, 'XData', mdl.part.position(:,1));
        drawnow;
    end
    
elseif mdl.dim==2
    
    if mdl.report.curr==1
        figure('Color', [1 1 1]); hold on; axis equal
        x = mdl.node.coord(:,1); y = mdl.node.coord(:,2);
        patch(x(mdl.elem.node)', y(mdl.elem.node)', 'white');
        MPM_PARTICLE_PLOT = plot(mdl.part.position(:,1), ...
            mdl.part.position(:,2), 'ro');
        hold off; drawnow;
    else
        set(MPM_PARTICLE_PLOT, 'XData', mdl.part.position(:,1));
        set(MPM_PARTICLE_PLOT, 'YData', mdl.part.position(:,2));
        drawnow;
    end
    
elseif mdl.dim==3
    
    if mdl.report.curr==1
        figure('Color', [1 1 1]); hold on; axis equal
        
        x_min = mdl.grid.min_coord;
        x_max = mdl.grid.max_coord;
        h = mdl.grid.elem_size;
        
        for j = x_min(2):h:x_max(2), for k = x_min(3):h:x_max(3)
                plot3([x_min(1) x_max(1)], [j j], [k k], 'k-'); end, end
        for i = x_min(1):h:x_max(1), for k = x_min(3):h:x_max(3)
                plot3([i i], [x_min(2) x_max(2)], [k k], 'k-'); end, end
        for i = x_min(1):h:x_max(1), for j = x_min(2):h:x_max(2)
                plot3([i i], [j j], [x_min(3) x_max(3)], 'k-'); end, end

        MPM_PARTICLE_PLOT = plot3(mdl.part.position(:,1), ...
            mdl.part.position(:,2), mdl.part.position(:,3), 'ro');
        hold off; drawnow;
    else
        set(MPM_PARTICLE_PLOT, 'XData', mdl.part.position(:,1));
        set(MPM_PARTICLE_PLOT, 'YData', mdl.part.position(:,2));
        set(MPM_PARTICLE_PLOT, 'ZData', mdl.part.position(:,3));
        drawnow;
    end
end
