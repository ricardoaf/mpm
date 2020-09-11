% mpm free fall example
%==========================================================================
function mpm_freefall
addpath(strrep(pwd,'examples','src'));

% input data
L = [1 1]; h = 1; gap = [-h -20*h; h h]; ppe = 2; supp = {[],[]};
E = 1e+6; nu = 0; rho = 2000; g = 9.81; tsim = 1.95; dtf = 0.1; nrep = 200;

% define model
mdl = mpm_boxdomain(L, h, gap, ppe, supp, E, nu, rho, g, tsim, dtf, nrep);
mdl.report.show_progress = false;

% run model and get reported time instants
out = mpm(mdl);
time = out.time;

% calc mpm and ref cm displacement
mpm_disp_cm = zeros(size(time));
ref_disp_cm = zeros(size(time));
for i = 1:length(out.time)
    mdl = out.model(i);
    mpm_disp_cm(i) = mdl.part.mass'*(mdl.part.position(:,2) ...
        - mdl.part.initial_position(:,2));
    ref_disp_cm(i) = -0.5*g*time(i)^2;
end
mpm_disp_cm = mpm_disp_cm ./ sum(mdl.part.mass);

% compare results
figure('Color','w'); hold on;
plot(time, ref_disp_cm, 'k-', 'DisplayName', 'Reference');
plot(time, mpm_disp_cm, 'ro', 'DisplayName', 'MPM');
xlabel('Time'); ylabel('Y'); title('MPM free fall'); legend;
