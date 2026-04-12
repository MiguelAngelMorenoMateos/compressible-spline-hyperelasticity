function [Aineq_blk, bineq_blk, Aeq_blk, beq_blk, lb_blk, ub_blk, x0_blk] = ...
    slice_linear_constraints(Aineq, bineq, Aeq, beq, lb, ub, x0, blockIdx)

    Aineq_blk = Aineq(:, blockIdx);
    bineq_blk = bineq;

    Aeq_blk = Aeq(:, blockIdx);
    beq_blk = beq;

    if ~isempty(Aineq_blk)
        mask = any(abs(Aineq_blk) > 1e-14, 2);
        Aineq_blk = Aineq_blk(mask,:);
        bineq_blk = bineq_blk(mask);
    end

    if ~isempty(Aeq_blk)
        mask = any(abs(Aeq_blk) > 1e-14, 2);
        Aeq_blk = Aeq_blk(mask,:);
        beq_blk = beq_blk(mask);
    end

    lb_blk = lb(blockIdx);
    ub_blk = ub(blockIdx);
    x0_blk = x0(blockIdx);
end