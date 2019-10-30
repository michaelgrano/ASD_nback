
clearvars;
load('data/nback_remOutliers.mat')

n_tp = length(aut(1).taskIrf_miss);

for dtype = {'hit', 'FA', 'miss'}

    % no distractor (not a)
    X = [];
    Y = [];

    for sub = 1:length(aut)
        if (aut(sub).hitRate == 1) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub for miss type 1')
            continue
        end
        if (aut(sub).faRate == 0) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub for FA type 1')
            continue
        end
        if any(isnan(aut(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub for %s type 1', dtype{1})
            continue
        end
        X = [X; aut(sub).(sprintf('taskIrf_%s', dtype{1}))];
        Y = [Y; 1];
    end
    for sub = 1:length(con)
        if (con(sub).hitRate == 1) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub for miss type 0')
            continue
        end
        if (con(sub).faRate == 0) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub for FA type 0')
            continue
        end
        if any(isnan(con(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub for %s type 0', dtype{1})
            continue
        end
        X = [X; con(sub).(sprintf('taskIrf_%s', dtype{1}))];
        Y = [Y; 0];
    end
    data = struct('X',X,'Y',Y);
    save(sprintf('data/XY_%s.mat', dtype{1}), 'data')

    % distractors (a)
    X = [];
    Y = [];

    for sub = 1:length(aut_a)
        if (aut_a(sub).hitRate == 1) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub for distractor miss type 1')
            continue
        end
        if (aut_a(sub).faRate == 0) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub for distractor FA type 1')
            continue
        end
        if any(isnan(aut_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub for distractor %s type 1', dtype{1})
            continue
        end
        X = [X; aut_a(sub).(sprintf('taskIrf_%s', dtype{1}))];
        Y = [Y; 1];
    end
    for sub = 1:length(con_a)
        if (con_a(sub).hitRate == 1) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub for distractor miss type 0')
            continue
        end
        if (con_a(sub).faRate == 0) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub for distractor FA type 0')
            continue
        end
        if any(isnan(con_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub for distractor %s type 0', dtype{1})
            continue
        end
        X = [X; con_a(sub).(sprintf('taskIrf_%s', dtype{1}))];
        Y = [Y; 0];
    end
    data = struct('X',X,'Y',Y);
    save(sprintf('data/XY_a_%s.mat', dtype{1}), 'data')
end

%% take difference
for dtype = {'hit', 'FA', 'miss'}

    % no distractor (not a)
    X = [];
    Y = [];

    for sub = 1:length(aut)
        if (aut(sub).hitRate == 1) || (aut_a(sub).hitRate == 1) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub for miss type 1')
            continue
        end
        if (aut(sub).faRate == 0) || (aut_a(sub).faRate == 0) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub for FA  type 1')
            continue
        end
        if any(isnan(aut(sub).(sprintf('taskIrf_%s', dtype{1})))) || any(isnan(aut_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub for %s type 1', dtype{1})
            continue
        end
        dat = mad(aut(sub).(sprintf('taskIrf_%s', dtype{1})),1);
        dat_a = mad(aut_a(sub).(sprintf('taskIrf_%s', dtype{1})),1);
%         X = [X; (dat-dat_a)./(dat+dat_a)];
        X = [X; dat_a-dat];
        Y = [Y; 1];
    end
    for sub = 1:length(con)
        if (con(sub).hitRate == 1) || (con_a(sub).hitRate == 1) && strcmp(dtype{1}, 'miss')
            fprintf('\n bad sub for miss type 0')
            continue
        end
        if (con(sub).faRate == 0) || (con_a(sub).faRate == 0) && strcmp(dtype{1}, 'FA')
            fprintf('\n bad sub for FA  type 0')
            continue
        end
        if any(isnan(con(sub).(sprintf('taskIrf_%s', dtype{1})))) || any(isnan(con_a(sub).(sprintf('taskIrf_%s', dtype{1}))))
            fprintf('\n bad sub for %s type 0', dtype{1})
            continue
        end
        dat = mad(con(sub).(sprintf('taskIrf_%s', dtype{1})),1);
        dat_a = mad(con_a(sub).(sprintf('taskIrf_%s', dtype{1})),1);
        X = [X; dat_a-dat];
%         X = [X; (dat-dat_a)./(dat+dat_a)];
        Y = [Y; 0];
    end
    data = struct('X',X,'Y',Y);
    save(sprintf('data/XY_distractors-diff_MAD_%s.mat', dtype{1}), 'data')
end
