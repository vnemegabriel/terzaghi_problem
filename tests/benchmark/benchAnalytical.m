function benchAnalytical()
% benchAnalytical  Micro-benchmark for analyticalTerzaghi.

here    = fileparts(mfilename('fullpath'));
project = fileparts(fileparts(here));
addpath(fullfile(project, 'src', 'postprocess'))

p.lambda = 4e7; p.mu = 4e7; p.alpha = 0.4; p.Mbiot = 4.17e7;
p.kperm = 1e-9; p.muf = 1; p.sigma0 = 10000; p.L = 10;

for nx = [50 200 1000]
    for nt = [1 10 100]
        x = linspace(0, p.L, nx)';
        t = linspace(1, 1000, nt);
        tic
        for r = 1:5
            [pa, ua] = analyticalTerzaghi(x, t, p); %#ok<ASGLU>
        end
        elapsed = toc / 5;
        fprintf('  Nx=%4d  Nt=%3d  call=%.4f s  (%.1f kpoints/s)\n', ...
            nx, nt, elapsed, nx*nt / elapsed / 1000)
    end
end
end
