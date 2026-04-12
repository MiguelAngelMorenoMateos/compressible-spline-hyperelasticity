function model = update_named_blocks(model, blockNames, xblk)
    idx = model.get_idx(model);
    pos = 0;

    for k = 1:numel(blockNames)
        name = lower(blockNames{k});

        switch name
            case 'additive'
                if model.use_I1
                    n = numel(idx.I1);
                    model.I1bar_vals = xblk(pos + (1:n));
                    pos = pos + n;
                end

                if model.use_I2
                    n = numel(idx.I2);
                    model.I2bar_vals = xblk(pos + (1:n));
                    pos = pos + n;
                end

                if model.use_J
                    n = numel(idx.J);
                    model.J_vals = xblk(pos + (1:n));
                    pos = pos + n;
                end

            case 'h'
                n = numel(idx.h);
                model.h_I1_vals = xblk(pos + (1:n));
                pos = pos + n;

            case 'g'
                n = numel(idx.g);
                model.g_J_vals = xblk(pos + (1:n));
                pos = pos + n;

            case 'j'
                n = numel(idx.j);
                model.j_I2_vals = xblk(pos + (1:n));
                pos = pos + n;

            case 'i'
                n = numel(idx.i);
                model.i_J_vals = xblk(pos + (1:n));
                pos = pos + n;

            case 'k'
                n = numel(idx.k);
                model.k_I1_vals = xblk(pos + (1:n));
                pos = pos + n;

            case 'l'
                n = numel(idx.l);
                model.l_I2_vals = xblk(pos + (1:n));
                pos = pos + n;

            otherwise
                error('Unknown block name: %s', blockNames{k});
        end
    end

    model = rebuild_all_splines(model);
end


function model = rebuild_all_splines(model)
    if model.use_I1
        model.WI1_spline = spapi(model.k, model.I1bar_knots, model.I1bar_vals);
        model.WI1_spline = make_fnxtr_spline(model.WI1_spline, model.k);
    end

    if model.use_I2
        model.WI2_spline = spapi(model.k, model.I2bar_knots, model.I2bar_vals);
        model.WI2_spline = make_fnxtr_spline(model.WI2_spline, model.k);
    end

    if model.use_J
        model.WJ_spline = spapi(model.k, model.J_knots, model.J_vals);
        model.WJ_spline = make_fnxtr_spline(model.WJ_spline, model.k);
    end

    if model.use_I1J
        model.h_spline = spapi(model.k, model.h_I1_knots, model.h_I1_vals);
        model.h_spline = make_fnxtr_spline(model.h_spline, model.k);

        model.g_spline = spapi(model.k, model.g_J_knots, model.g_J_vals);
        model.g_spline = make_fnxtr_spline(model.g_spline, model.k);
    end

    if model.use_I2J
        model.j_spline = spapi(model.k, model.j_I2_knots, model.j_I2_vals);
        model.j_spline = make_fnxtr_spline(model.j_spline, model.k);

        model.i_spline = spapi(model.k, model.i_J_knots, model.i_J_vals);
        model.i_spline = make_fnxtr_spline(model.i_spline, model.k);
    end

    if model.use_I1I2
        model.k_spline = spapi(model.k, model.k_I1_knots, model.k_I1_vals);
        model.k_spline = make_fnxtr_spline(model.k_spline, model.k);

        model.l_spline = spapi(model.k, model.l_I2_knots, model.l_I2_vals);
        model.l_spline = make_fnxtr_spline(model.l_spline, model.k);
    end
end

function spline = make_fnxtr_spline(spline, order)
    spline = fnxtr(spline, order);
end