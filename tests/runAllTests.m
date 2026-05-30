function results = runAllTests()
% runAllTests  Run every test in the tests/ folder and print a summary.
%
%   results = runAllTests();
%
% Each *.m file beginning with 'test' is discovered automatically and
% executed via MATLAB's function-based test framework (functiontests).

here    = fileparts(mfilename('fullpath'));
project = fileparts(here);

addpath(project)
addpath(fullfile(project, 'src', 'elements'))
addpath(fullfile(project, 'src', 'assembly'))
addpath(fullfile(project, 'src', 'solvers'))
addpath(fullfile(project, 'src', 'postprocess'))
addpath(fullfile(project, 'src', 'terzaghi_funs'))
addpath(fullfile(project, 'src', 'utils'))
addpath(fullfile(project, 'meshes'))
addpath(here)

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.TestRunProgressPlugin

suite  = TestSuite.fromFolder(here);
runner = TestRunner.withTextOutput('OutputDetail', 1);

fprintf('\n========================================\n')
fprintf(' TP3 Poroelasticity - Test Suite\n')
fprintf('========================================\n')
fprintf(' Tests discovered: %d\n', numel(suite))
fprintf('----------------------------------------\n\n')

results = runner.run(suite);

% --- summary table ---
nPassed   = sum([results.Passed]);
nFailed   = sum([results.Failed]);
nIncomp   = sum([results.Incomplete]);
totalTime = sum([results.Duration]);

fprintf('\n========================================\n')
fprintf(' Summary\n')
fprintf('========================================\n')
fprintf('   Passed     : %d\n', nPassed)
fprintf('   Failed     : %d\n', nFailed)
fprintf('   Incomplete : %d\n', nIncomp)
fprintf('   Total time : %.2f s\n', totalTime)
fprintf('========================================\n\n')

if nFailed > 0
    fprintf('FAILED TESTS:\n');
    for i = 1:numel(results)
        if results(i).Failed
            fprintf('   - %s\n', results(i).Name)
        end
    end
end
end
