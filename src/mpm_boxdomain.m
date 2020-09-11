% generate mpm model using box domain
%==========================================================================
function mdl = mpm_boxdomain (len, elem_size, gap, partperelem, supp, ...
    young_modulus, poisson_ratio, density, gravity_acc, simulation_time,...
    dt_fraction, n_report, show_progress, show_animation)

% default values
if ~exist('dt_fraction', 'var'), dt_fraction=0.01; end
if ~exist('n_report', 'var'), n_report=200; end
if ~exist('show_progress', 'var'), show_progress=true; end
if ~exist('show_animation', 'var'), show_animation=true; end

% number of dimensions
%--------------------------------------------------------------------------
dim = length(len);

% boundbox
%--------------------------------------------------------------------------
grid = {};
grid.elem_size = elem_size;
grid.min_coord = gap(1,:);
grid.max_coord = len + gap(2,:);

% cells
%--------------------------------------------------------------------------
grid.n_cell = ceil((grid.max_coord-grid.min_coord)./elem_size);
grid.max_coord = grid.min_coord + grid.n_cell*elem_size;

% node coords
%--------------------------------------------------------------------------
node = {}; node.n = prod(grid.n_cell+1);
X = cell(1,3); X(:) = {0};
for i=1:dim
    X{i} = grid.min_coord(i):elem_size:grid.max_coord(i);
end
[X{2},X{1},X{3}] = meshgrid(X{2},X{1},X{3});
node.coord = [X{1}(:) X{2}(:) X{3}(:)];
node.coord = node.coord(:,1:dim);

% elem
%--------------------------------------------------------------------------
elem = {};
elem.n = prod(grid.n_cell);
elem.n_node = 2^dim;
elem.node = zeros(elem.n, elem.n_node);

% connectivities
%--------------------------------------------------------------------------
n_cell = zeros(1,3); n_cell(1:dim) = grid.n_cell;
node_matrix = reshape(1:node.n, n_cell+1);
n_node = n_cell + 1; n_cell(n_cell==0) = 1;
idx_p = 1:elem.n_node;
for m = 1:dim-1, idx_p([4*m-1 4*m]) = idx_p([4*m 4*m-1]); end
for e = 1:elem.n
    [i,j,k] = ind2sub(n_cell, e);
    idx_i = [i i+1]; idx_i(idx_i>n_node(1)) = [];
    idx_j = [j j+1]; idx_j(idx_j>n_node(2)) = [];
    idx_k = [k k+1]; idx_k(idx_k>n_node(3)) = [];
    e_node = node_matrix(idx_i,idx_j,idx_k);
    elem.node(e,:) = e_node(idx_p);
end

% supports
%--------------------------------------------------------------------------
node.support = sparse(node.n, dim);
for i = 1:dim
    node.support(support_node(node.coord, supp{i}), i) = 1;
end

% material
%--------------------------------------------------------------------------
mat = {};
mat.young_modulus = young_modulus;
mat.poisson_ratio = poisson_ratio;
mat.density = density;

% time data
%--------------------------------------------------------------------------
time = {};
time.dt_fraction = dt_fraction;
time.simulation = simulation_time;

% report data
%--------------------------------------------------------------------------
report = {};
report.n = n_report;
report.show_progress = show_progress;
report.show_animation = show_animation;

% particles
%--------------------------------------------------------------------------
part = {};
part.n = elem.n*partperelem^dim;
part.position = zeros(part.n,dim);

% use gauss quadrate for particle positions
gaussquad = {0, 1/sqrt(3)*[-1 1], sqrt(3/5)*[-1 0 1], ...
    sort([sqrt(3/7+2/7*sqrt(6/5))*[-1 1] sqrt(3/7-2/7*sqrt(6/5))*[-1 1]])};

% generate parametric positions
R = cell(1,3); R(:)={0};
for i = 1:dim, R{i} = gaussquad{partperelem}; end
[R{1},R{2},R{3}]=meshgrid(R{1},R{2},R{3});
r_coord=[R{1}(:) R{2}(:) R{3}(:)];
lbd = (r_coord(:,1:dim)+1)/2;

% generate particles and calc their positions at active elements
pid = 0; nea = 0;
for e = 1:elem.n
    e_coord = node.coord(elem.node(e,:),:);
    active_elem = true;
    for i = 1:dim
        if min(e_coord(:,i))<0 || max(e_coord(:,i))>len(i)
            active_elem = false; break;
        end
    end
    if active_elem
        nea = nea + 1;
        np = size(e_coord,1);
        part.position(pid+1:pid+np,:) = min(e_coord) + ...
            lbd.*(max(e_coord)-min(e_coord));
        pid = pid + np;
    end
end
part.n = pid;
part.position(part.n+1:end,:) = [];

% init other particle fields
part.velocity = zeros(part.n, dim);
part.density = ones(part.n, 1)*density;
part.mass = ones(part.n, 1)*density*(elem_size/partperelem)^2;
part.body_force = zeros(part.n, dim);
if dim>1, gravity_idx = 2; else, gravity_idx = 1; end 
part.body_force(:,gravity_idx) = -density*abs(gravity_acc);
part.stress = cell(1, part.n);
part.stress(:) = {zeros(dim)};

% define struct model
%--------------------------------------------------------------------------
mdl = {};
mdl.dim = dim;
mdl.grid = grid;
mdl.node = node;
mdl.elem = elem;
mdl.mat = mat;
mdl.part = part;
mdl.time = time;
mdl.report = report;

% define support nodes
%==========================================================================
function node = support_node (node_coord, supp_coord)

% init support node array
node = [];
node_coord(:,size(node_coord,2)+1:3) = 0;

% loop on support coords
for s = supp_coord', s_node = [];
    x = nan(1,3); x(1:length(s)) = s;
    
    % identify support nodes based on its coords
    if isnan(x(1)) && isnan(x(2)) && ~isnan(x(3))
        s_node = find(node_coord(:,3)==x(3));
    elseif isnan(x(1)) && ~isnan(x(2)) && isnan(x(3))
        s_node = find(node_coord(:,2)==x(2));
    elseif ~isnan(x(1)) && isnan(x(2)) && isnan(x(3))
        s_node = find(node_coord(:,1)==x(1));
    elseif isnan(x(1)) && ~isnan(x(2)) && ~isnan(x(3))
        s_node = find(node_coord(:,2)==x(2) && node_coord(:,3)==x(3));
    elseif ~isnan(x(1)) && isnan(x(2)) && ~isnan(x(3))
        s_node = find(node_coord(:,1)==x(1) && node_coord(:,3)==x(3));
    elseif ~isnan(x(1)) && ~isnan(x(2)) && isnan(x(3))
        s_node = find(node_coord(:,1)==x(1) && node_coord(:,2)==x(2));
    elseif ~isnan(x(1)) && ~isnan(x(2)) && ~isnan(x(3))
        s_node = find(node_coord(:,1)==x(1) && node_coord(:,2)==x(2) && ...
            node_coord(:,3)==x(3));
    end
    
    % append identified nodes into support node array
    node = [node; s_node];
end
