function schedule = build_alternating_schedule(model)

    schedule = {};
    hasAdditive = ~isempty(model.get_idx(model).additive);

    % ============================================================
    % BLOCK 1
    % ============================================================
    block1 = {};
    if model.use_I1J
        block1 = [block1, {'h'}];
    end
    if model.use_I2J
        block1 = [block1, {'i'}];
    end
    if model.use_I1I2
        block1 = [block1, {'k'}];
    end
    if hasAdditive
        schedule{end+1} = [{'additive'}, block1];
    else
        schedule{end+1} = block1;
    end

    % ============================================================
    % BLOCK 2
    % ============================================================
    block2 = {};
    if model.use_I1J
        block2 = [block2, {'g'}];
    end
    if model.use_I2J
        block2 = [block2, {'j'}];
    end
    if model.use_I1I2
        block2 = [block2, {'l'}];
    end
    if hasAdditive
        schedule{end+1} = [{'additive'}, block2];
    else
        schedule{end+1} = block2;
    end

    for b = 1:numel(schedule)
        fprintf('Block %d: ', b);
        block = schedule{b};  % current block
        if isempty(schedule{b})
            fprintf('(empty)\n');
            continue;
        end
        % print all terms in the block
        fprintf('%s ', block{:});
        fprintf('\n');
    end
end