% mpm column under self-weight example (tielen 4.2)
%==========================================================================
function mpm_column_selfweight
addpath(strrep(pwd,'examples','src'));

% input data
L = [1 25]; h = 1; gap = [-h -h; h h]; ppe = 2; supp = {[],[NaN 0]};
E = 5e+4; nu = 0; rho = 1; g = 9.81; tsim = 0.9; dtf = 0.3; nrep = 200;

% define model
mdl = mpm_boxdomain(L, h, gap, ppe, supp, E, nu, rho, g, tsim, dtf, nrep);
mdl.show_animation = false;

% run model and get reported time instants
out = mpm(mdl);
time = out.time;

% calc mpm cm velocity
mpm_vel_cm = zeros(size(time));
for i = 1:length(out.time)
    mdl = out.model(i);
    mpm_vel_cm(i) = mdl.part.mass'*mdl.part.velocity(:,2);
end
mpm_vel_cm = mpm_vel_cm ./ sum(mdl.part.mass);

% calc ref cm velocity
ref_vel_cm = ref_velocity_cm (L(2), E, rho, g, 0, 5, time);

% compare results
figure('Color','w'); hold on; title('MPM self-weight column');
plot(time, ref_vel_cm, 'k-', 'DisplayName', 'Reference');
plot(time, mpm_vel_cm, 'ro', 'DisplayName', 'MPM');
xlabel('Time'); ylabel('CM Velocity Y'); legend;

%==========================================================================
function veloc = ref_velocity_cm (H, E, rho, g, p0, N, t)
veloc = zeros(size(t));
c = -16*sqrt(E/rho)/E/(pi^3);
for n = 0:N-1
    veloc = veloc + c*sin(.5*pi*(1+2*n)*sqrt(E/rho)/H*t).*((-1)^(n+1).* ...
        pi*n*p0+H*g*rho-0.5*p0*(-1)^n*pi)./((1+2*n)^3);
end
