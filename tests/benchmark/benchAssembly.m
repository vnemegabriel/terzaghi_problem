function benchAssembly()
% benchAssembly  Micro-benchmark for buildGlobalMatrices on a large mesh.
%   Reports total assembly time and per-element cost.

here    = fileparts(mfilename('fullpath'));
project = fileparts(fileparts(here));
addpath(project)
addpath(fullfile(project, 'src', 'elements'))
addpath(fullfile(project, 'src', 'assembly'))
addpath(fullfile(project, 'src', 'terzaghi_funs'))
addpath(fullfile(project, 'src', 'utils'))
addpath(fullfile(project, 'meshes'))

p.lambda = 4e7; p.mu = 4e7; p.alpha = 0.4; p.Mbiot = 4.17e7;
p.kperm = 1e-9; p.muf = 1;

for eleType = {'Q4', 'Q8'}
    fprintf('\n--- %s ---\n', eleType{1})
    for ny = [40 80 160]
        mesh = generateColumnMesh(1, 10, 4, ny, eleType{1});
        tic
        [K, H, Q, S] = buildGlobalMatrices(mesh, [], p, 2); %#ok<ASGLU>
        t = toc;
        nElem = numel(mesh.elements);
        fprintf('  ny=%3d  nElem=%4d  nNod=%5d  total=%.3f s  per-elem=%.2f us\n', ...
            ny, nElem, size(mesh.nodes,1), t, t/nElem*1e6)
    end
end
end
