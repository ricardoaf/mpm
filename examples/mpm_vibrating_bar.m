% mpm vibrating bar example
%==========================================================================
function mpm_vibrating_bar
addpath(strrep(pwd,'examples','src'));

% input data
L = [25 1]; h = 1; gap = [-h -h; h h]; ppe = 2;
E = 100; nu = 0; rho = 1; g = 0; tsim = 18; dtf = 0.1; nrep = 200;

% choose example
if 1 % tielen 4.1
    supp = {[0 NaN; L(1) NaN],[]};
    beta = pi/L(1); fac = 2;
else % bardenhagen (2002) 4.2
    supp = {[0 NaN],[]}; fac = 1;
    n = 1; beta = 0.5*(2*n-1)*pi/L(1);
end

% define model
mdl = mpm_boxdomain(L, h, gap, ppe, supp, E, nu, rho, g, tsim, dtf, nrep);

% set prescribed velocity
v0 = 0.1; veloc = @(x) v0*sin(beta*x);
for p = 1:mdl.part.n
    mdl.part.velocity(p,1) = veloc(mdl.part.position(p,1));
end

% run model and get reported time instants
out = mpm(mdl);
time = out.time;

% calc mpm cm velocity
mpm_vel_cm = zeros(size(time));
for i = 1:length(out.time)
    mdl = out.model(i);
    mpm_vel_cm(i) = mdl.part.mass'*mdl.part.velocity(:,1);
end
mpm_vel_cm = mpm_vel_cm ./ sum(mdl.part.mass);

% calc ref cm velocity
ref_vel_cm = fac*v0/(beta*L(1))*cos(beta*sqrt(E/rho)*time);

% compare results
figure('Color','w'); hold on;
plot(time, ref_vel_cm, 'k-', 'DisplayName', 'Reference');
plot(time, mpm_vel_cm, 'ro', 'DisplayName', 'MPM');
xlabel('Time'); ylabel('CM Velocity X'); title('MPM vibrating bar');legend;
