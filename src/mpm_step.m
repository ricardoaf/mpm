% mpm step
%==========================================================================
function mdl = mpm_step (mdl, stress_update)

% USL as default stress update scheme
if nargin<2 || ~any(strcmpi({'USL', 'USF', 'MUSL'}, stress_update))
    stress_update = 'USL';
end

% 1. Calc nodal mass and momentum
%--------------------------------------------------------------------------
node_mass = zeros(mdl.node.n,1);
node_momentum = zeros(mdl.node.n, mdl.dim);
for I = mdl.node.active
    mI = 0;
    for p = mdl.node.part{I}
        mp = mdl.part.mass(p);
        NIp = mdl.elem.shape_fcn{p}(I);
        mI = mI + mp*NIp;
    end
    node_mass(I) = mI;
    for i = 1:mdl.dim
        piI = 0;
        % 2. Impose essential boundary conditions
        %------------------------------------------------------------------
        if mdl.node.support(I,i)==0
            for p = mdl.node.part{I}
                mp = mdl.part.mass(p);
                NIp = mdl.elem.shape_fcn{p}(I);
                vip = mdl.part.velocity(p,i);
                piI = piI + mp*vip*NIp;
            end
        end
        node_momentum(I,i) = piI;
    end
end

% 3. For USF only, update particle density and stress
%--------------------------------------------------------------------------
if strcmpi(stress_update, 'USF')
    mdl = stress_state (mdl, node_mass, node_momentum);
end

% 4. Calc nodal internal and external forces
%--------------------------------------------------------------------------
node_internalforce = zeros(mdl.node.n, mdl.dim);
node_externalforce = zeros(mdl.node.n, mdl.dim);
for I = mdl.node.active
    for i = 1:mdl.dim
        fint_iI = 0;
        fext_iI = 0;
        if mdl.node.support(I,i)==0
            for p = mdl.node.part{I}
                mp = mdl.part.mass(p);
                rhop = mdl.part.density(p);
                for j = 1:mdl.dim
                    NIpj = mdl.elem.dshape_fcn{p}(I,j);
                    sijp = mdl.part.stress{p}(i,j);
                    fint_iI = fint_iI - NIpj*sijp*mp/rhop;
                end
                NIp = mdl.elem.shape_fcn{p}(I);
                bip = mdl.part.body_force(p,i);
                fext_iI = fext_iI + NIp*bip*mp/rhop;
            end
        end
        node_internalforce(I,i) = fint_iI;
        node_externalforce(I,i) = fext_iI;
    end
end
node_force = node_internalforce + node_externalforce;

% 5. Integrate nodal momentum
%--------------------------------------------------------------------------
dt = mdl.time.incr;
node_momentum = node_momentum + node_force * dt;

% 6. Update particle velocity and position
%--------------------------------------------------------------------------
for p = 1:mdl.part.n
    for i = 1:mdl.dim
        delta_vip = 0;
        delta_xip = 0;
        for I = mdl.part.node{p}
            fiI = node_force(I,i);
            piI = node_momentum(I,i);
            NIp = mdl.elem.shape_fcn{p}(I);
            mI = node_mass(I);
            delta_vip = delta_vip + fiI*NIp/mI*dt;
            delta_xip = delta_xip + piI*NIp/mI*dt;
        end
        mdl.part.velocity(p,i) = mdl.part.velocity(p,i) + delta_vip;
        mdl.part.position(p,i) = mdl.part.position(p,i) + delta_xip;
    end
end

% 7. For MUSL only, recalc nodal momentum
%--------------------------------------------------------------------------
node_momentum = zeros(mdl.node.n, mdl.dim);
for I = mdl.node.active
    for i = 1:mdl.dim
        piI = 0;
        % Impose essential boundary conditions
        %------------------------------------------------------------------
        if mdl.node.support(I,i)==0
            for p = mdl.node.part{I}
                mp = mdl.part.mass(p);
                NIp = mdl.elem.shape_fcn{p}(I);
                vip = mdl.part.velocity(p,i);
                piI = piI + mp*vip*NIp;
            end
        end
        node_momentum(I,i) = piI;
    end
end

% 8. For MUSL or USL only, update particle density and stress
%--------------------------------------------------------------------------
if any(strcmpi({'USL', 'MUSL'}, stress_update))
    mdl = stress_state (mdl, node_mass, node_momentum);
end

% Update stress state
%==========================================================================
function mdl = stress_state (mdl, node_mass, node_momentum)

% a. Calc nodal velocity
%--------------------------------------------------------------------------
node_velocity = node_momentum ./ repmat(node_mass,[1 mdl.dim]);

% b. Calc particle velocity gradient, deformation rate and spin tensor
%--------------------------------------------------------------------------
velocity_grad = cell(1, mdl.part.n);
deformation_rate = cell(1, mdl.part.n);
spin_tensor = cell(1, mdl.part.n);
for p = 1:mdl.part.n
    for i = 1:mdl.dim
        for j = 1:mdl.dim
            Lij = 0;
            for I = mdl.part.node{p}
                NIpj = mdl.elem.dshape_fcn{p}(I,j);
                viI = node_velocity(I,i);
                Lij = Lij + NIpj*viI;
            end
            velocity_grad{p}(i,j) = Lij;
        end
    end
    deformation_rate{p} = 0.5*(velocity_grad{p} + velocity_grad{p}');
    spin_tensor{p} = 0.5*(velocity_grad{p} - velocity_grad{p}');
end

% c. Update particle density
%--------------------------------------------------------------------------
dt = mdl.time.incr;
for p = 1:mdl.part.n
    deij = deformation_rate{p} * dt;
    mdl.part.density(p) = mdl.part.density(p)/(1 + trace(deij));
end

% d. Update particle stress state
%--------------------------------------------------------------------------
for p = 1:mdl.part.n
    % Update deformation gradient
    mdl.part.deformation_grad{p} = (eye(mdl.dim) + ...
        velocity_grad{p}*dt) * mdl.part.deformation_grad{p};
    J = det(mdl.part.deformation_grad{p});
    % Elastic properties
    E = mdl.mat.young_modulus; nu = mdl.mat.poisson_ratio;
    mu = E/2/(1+nu); lambda = E*nu/(1+nu)/(1-2*nu);
    % Calc Jaumann rate
    jaumann_kichhoff_rate = 2*mu*deformation_rate{p} + ...
        lambda*eye(mdl.dim)*trace(deformation_rate{p});
    jaumann_rate = 1/J*jaumann_kichhoff_rate;
    % Calc rotated Cauchy stress
    last_cauchy = mdl.part.stress{p};
    rotated_cauchy = last_cauchy + (last_cauchy * spin_tensor{p}' + ...
        spin_tensor{p} * last_cauchy) * dt;
    % Update Cauchy stress
    mdl.part.stress{p} = rotated_cauchy + jaumann_rate * dt;
end
