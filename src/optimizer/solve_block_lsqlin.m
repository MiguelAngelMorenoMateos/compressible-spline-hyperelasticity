function [model, out] = solve_block_lsqlin(model, prot, data, Aineq, bineq, Aeq, beq, x0, opts_lin, blockNames)

    [A, y, blockIdx] = build_linear_system_block(model, prot, data, blockNames);

    lb = model.lb();
    ub = model.ub();

    [Aineq_blk, bineq_blk, Aeq_blk, beq_blk, lb_blk, ub_blk, x0_blk] = ...
        slice_linear_constraints(Aineq, bineq, Aeq, beq, lb, ub, x0, blockIdx);

    [xblk, resnorm, residual, exitflag, output, lambda] = ...
        lsqlin(A, y, Aineq_blk, bineq_blk, Aeq_blk, beq_blk, lb_blk, ub_blk, x0_blk, opts_lin);

    model = update_named_blocks(model, blockNames, xblk);

    out = struct();
    out.xblk = xblk;
    out.resnorm = resnorm;
    out.residual = residual;
    out.exitflag = exitflag;
    out.output = output;
    out.lambda = lambda;
end