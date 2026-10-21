function baseline_dir = make_baseline(source)

% function to save a snapshot of the pipeline code, used as the OLD version (reference) by "run_regression"
% (revised October 2026)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% source:               (optional) git reference (tag, branch or commit) to take the code from
%                       when it is not provided, the current files of the repository are copied (including uncommitted changes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   OUTPUT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% baseline_dir:         folder with the snapshot (<work_root>/code_baseline); it also contains "baseline_info.txt"
% WARNING: the chain of intermediate results (<work_root>/chain) is deleted, because it was generated with the previous baseline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = regression_config();
baseline_dir = fullfile(cfg.work_root, 'code_baseline');
if isfolder(baseline_dir), rmdir(baseline_dir, 's'); end                    % replace the previous baseline
mkdir(baseline_dir);
chain_dir = fullfile(cfg.work_root, 'chain');
if isfolder(chain_dir), rmdir(chain_dir, 's'); end                          % the intermediate results must be regenerated with the new baseline

if nargin < 1 || strlength(string(source)) == 0                             % snapshot of the current files (only the top folder, where the pipeline functions are)
    copyfile(fullfile(cfg.repo_dir, '*.m'), baseline_dir);
    [~, commit] = system(sprintf('git -C "%s" rev-parse HEAD', cfg.repo_dir));
    [~, changes] = system(sprintf('git -C "%s" status --short -- "*.m"', cfg.repo_dir));
    info = sprintf('Snapshot of the working folder (commit %s), taken on %s\nUncommitted changes:\n%s', strtrim(commit), string(datetime('now')), changes);
else                                                                        % code from a git reference
    status = system(sprintf('git -C "%s" archive "%s" -- "*.m" | tar -x -C "%s"', cfg.repo_dir, source, baseline_dir));
    if status ~= 0
        error('The git reference "%s" could not be exported', source);
    end
    info = sprintf('Code from the git reference %s, exported on %s\n', source, string(datetime('now')));
end
writelines(string(info), fullfile(baseline_dir, 'baseline_info.txt'));
fprintf('%s\nBaseline saved in %s\n', info, baseline_dir);

end
