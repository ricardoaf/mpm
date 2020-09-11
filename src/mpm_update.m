% mpm update function
%==========================================================================
function mdl = mpm_update (mdl, time)

% Init active node flags
is_active_node = zeros(1,mdl.node.n);

% Init node-part map data
mdl.node.part = cell(1, mdl.node.n);
mdl.part.node = cell(1, mdl.part.n);

% Init shape function and its derivatives
mdl.elem.shape_fcn = cell(1, mdl.part.n);
mdl.elem.shape_fcn(:) = {zeros(mdl.node.n,1)};
mdl.elem.dshape_fcn = cell(1, mdl.part.n);
mdl.elem.dshape_fcn(:) = {sparse(mdl.node.n, mdl.dim)};

% Init max speed sum
max_cpvp = 0;

% Loop on particles
for p = 1:mdl.part.n
    
    % Get current element and particle parametric coords
    [part_elem, param_coord] = get_elem(p, mdl);
    
    % Get current element nodes and set them as active nodes
    part_conn = mdl.elem.node(part_elem,:);
    is_active_node(part_conn) = 1;
    
    % Update node-part map
    for I = part_conn
        mdl.node.part{I} = [mdl.node.part{I} p];
    end
    mdl.part.node{p} = part_conn;
    
    % Calc shape function values and its derivatives
    shape_function = ['shape_function' num2str(mdl.dim)];
    [N, dN] = feval(shape_function, param_coord, mdl.grid.elem_size);
    mdl.elem.shape_fcn{p}(part_conn) = N';
    mdl.elem.dshape_fcn{p}(part_conn,:) = dN;
    
    % Check for max speed sum update
    [cp, vp] = particle_speed (mdl, p);
    if (cp+vp)>max_cpvp, max_cpvp=cp+vp; end
    
end
% Define active nodes
mdl.node.active = find(is_active_node);

% Update time increment
dtcr = mdl.grid.elem_size/max_cpvp;
mdl.time.incr = min(mdl.time.dt_fraction * dtcr, mdl.time.simulation-time);

% mpm get part_elem function
%==========================================================================
function [elem, param_coord] = get_elem (part, mdl)

% Get particle position and grid minimum coords and elemsize
pos = mdl.part.position(part,:);
pos_min = mdl.grid.min_coord;
d = mdl.grid.elem_size;

% Calculate particle location indices
index = floor((pos-pos_min)./d) + 1;
index(:,mdl.dim+1:3) = 1;

% Calc elem and its particle parametric coords
elem = zeros(size(part));
param_coord = zeros(length(part), mdl.dim);
n_cell = ones(1,3); n_cell(1:mdl.dim) = mdl.grid.n_cell;
for p = 1:length(part)
    elem(p) = sub2ind(n_cell,index(p,1),index(p,2),index(p,3));
    % particle parametric coords on elem
    ecoord = mdl.node.coord(mdl.elem.node(elem(p),:),:);
    param_coord(p,:) = 2*(pos-min(ecoord))./(max(ecoord)-min(ecoord))-1;
end

% calc particle speed
%==========================================================================
function [cp, vp] = particle_speed (mdl, p)

% Get elastic parameters
mat = mdl.mat; E = mat.young_modulus; nu = mat.poisson_ratio;

% Calc element wave and particle speed
cp = sqrt(E*(1-nu)/(1+nu)/(1-2*nu)/mat.density);
vp = norm(mdl.part.velocity(p,:));

% mpm shape functions: L2, Q4, B8
%==========================================================================
function [N, dN] = shape_function1 (param_coord, elem_size)
r = param_coord(1);
N = 0.5*[(1-r); (1+r)];
dN = 0.5*[-1; +1]*2/elem_size;

function [N, dN] = shape_function2 (param_coord, elem_size)
r = param_coord(1); s = param_coord(2);
N = 0.25*[(1-r)*(1-s) (1+r)*(1-s) (1+r)*(1+s) (1-r)*(1+s)]';
dN = 0.25*[-1+s -1+r; 1-s -1-r; 1+s 1+r; -1-s 1-r]*2/elem_size;

function [N, dN] = shape_function3 (param_coord, elem_size)
r = param_coord(1); s = param_coord(2); t = param_coord(3);
N = 1/8*[-r*s*t+r*s+r*t+s*t-r-s-t+1, r*s*t-r*s-r*t+s*t+r-s-t+1, ...
    -r*s*t+r*s-r*t-s*t+r+s-t+1, r*s*t-r*s+r*t-s*t-r+s-t+1, ...
    r*s*t+r*s-r*t-s*t-r-s+t+1, -r*s*t-r*s+r*t-s*t+r-s+t+1, ...
    r*s*t+r*s+r*t+s*t+r+s+t+1, -r*s*t-r*s-r*t+s*t-r+s+t+1];
dN = 2/elem_size * 1/8 * ...
    [-s*t+s+t-1, -r*t+r+t-1, -r*s+r+s-1; s*t-s-t+1, r*t-r+t-1, r*s-r+s-1;...
    -s*t+s-t+1, -r*t+r-t+1, -r*s-r-s-1; s*t-s+t-1, r*t-r-t+1, r*s+r-s-1;...
    s*t+s-t-1, r*t+r-t-1, r*s-r-s+1; -s*t-s+t+1, -r*t-r-t-1, -r*s+r-s+1;...
    s*t+s+t+1, r*t+r+t+1, r*s+r+s+1; -s*t-s-t-1, -r*t-r+t+1, -r*s-r+s+1];