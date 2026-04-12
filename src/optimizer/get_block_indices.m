function blockIdx = get_block_indices(model, blockNames)
    idx = model.get_idx(model);
    blockIdx = [];

    for k = 1:numel(blockNames)
        name = lower(blockNames{k});
        switch name
            case 'additive'
                blockIdx = [blockIdx, idx.additive];
            case 'h'
                blockIdx = [blockIdx, idx.h];
            case 'g'
                blockIdx = [blockIdx, idx.g];
            case 'j'
                blockIdx = [blockIdx, idx.j];
            case 'i'
                blockIdx = [blockIdx, idx.i];
            case 'k'
                blockIdx = [blockIdx, idx.k];
            case 'l'
                blockIdx = [blockIdx, idx.l];
            otherwise
                error('Unknown block name: %s', blockNames{k});
        end
    end

    blockIdx = unique(blockIdx, 'stable');
end