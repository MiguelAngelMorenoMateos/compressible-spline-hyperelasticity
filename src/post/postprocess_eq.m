% ===================== postprocess_eq.m =====================
function postprocess_eq(result_or_file)
% Postprocess equilibrium fit results: plot spline activations
%
% Usage:
%   postprocess_eq('fit_result_eq.mat')
%   postprocess_eq(result_struct)

% Uncoupled:
% Color_B = '#7556AD';
% w/ coupling term I1,J
%Color_B = '#9C6345';
% w/ coupling term I2,J
%Color_B = '#356823';
% w/ coupling term I1,I2
%Color_B = '#34689E';
% w/ coupling term I1,I2 & without other coupling terms
Color_B = '#850000';


Fontsize = 12.1;
Fontsize_axes = 11;

    if ischar(result_or_file) || isstring(result_or_file)
        tmp = load(result_or_file);
        if isfield(tmp,'result')
            result = tmp.result;
        else
            error('MAT-file does not contain variable "result".');
        end
    else
        result = result_or_file;
    end

    % load components
    x_opt    = result.x_opt;
    residual = result.residual;
    model    = result.model;
    prot     = result.prot;
    data     = result.data;


    % ----------------------------
    % 1) Plot fits to data
    % ----------------------------
    % Initialize storage for UT and UC modes
    P11_all_data   = [];  % measured P11
    P11_all_fit    = [];  % model fit P11
    P22_all = [];        % will store all P22 values
    lambda_all = [];     % corresponding lambda1 values
    modes_all = {};      % optional: keep track of mode names
    for i = 1:numel(data.modes)
    %for i=1:0
        mode = data.modes{i};
        %print mode:
        fprintf('Mode: %s\n', mode);
        D = data.(mode);
        p = prot.(mode);

        Pfit = zeros(size(D.P11));
        P22fit = zeros(size(D.P11));
        lam2 = zeros(size(D.lambda));
        detJ = zeros(size(D.lambda));
        I1bar = zeros(size(D.lambda));
        I2bar = zeros(size(D.lambda));
        for k = 1:numel(D.P11)
            F = p.Fbar(D.lambda(k)); % Initial guess is F for incompressible case, then NR will adjust lambda2 and 3 to satisfy P(2,2)=0.
            [a,b,c] = model.eval_ab(model, F);
            P = eq_P_from_F(a,b,c,F);
            if strcmpi(mode,'SS')
                Pfit(k) = P(1,2); % for SS, we want to fit P12 instead of P11
            else
                Pfit(k) = P(1,1);
            end
            P22fit(k) = P(2,2);
            lam2(k) = F(2,2); % This time, I need to save the adjusted lambda2 values to plot them and to see how far they are from the incompressible assumption. I can also compute the effective J = det(F) to see how much volume change is happening.
            detJ(k) = det(F);
            C  = F.' * F;
            I1 = trace(C);
            I2 = 0.5*(I1^2 - trace(C*C));
            I1bar(k) = I1 * (sqrt(det(C)))^(-2/3);
            I2bar(k) = I2 * (sqrt(det(C)))^(-4/3);
        end

        % print statistics for I1bar and I2bar and det(F) ranges
        fprintf('  I1bar range: [%.4f, %.4f]\n', min(I1bar), max(I1bar));
        fprintf('  I2bar range: [%.4f, %.4f]\n', min(I2bar), max(I2bar));
        fprintf('  det(F) range: [%.4f, %.4f]\n', min(detJ), max(detJ));

        % print integral of g(J) and h(I1) if they exist to check their scales
        if isfield(model, 'g_spline')
            I = fnint(model.g_spline); % integral of g(J) spline
            g_integral = fnval(I, model.g_J_knots(end)) - fnval(I, model.g_J_knots(1));
            fprintf('  Integral of g(J) over domain: %.4f\n', g_integral);
        end
        if isfield(model, 'h_spline')
            I = fnint(model.h_spline); % integral of h(I1) spline
            h_integral = fnval(I, model.h_I1_knots(end)) - fnval(I, model.h_I1_knots(1));
            fprintf('  Integral of h(I1) over domain: %.4f\n', h_integral);
        end

        % print spline slopes at J=1 for WJ and g(J) if they exist
        if isfield(model, 'WJ_spline')
            fprintf('  Spline slope at J=1 for WJ: %.4f\n', fnval(fnder(model.WJ_spline, 1), 1));
        end
        if isfield(model, 'g_spline')
            fprintf('  Spline slope at J=1 for g(J): %.4f\n', fnval(fnder(model.g_spline, 1), 1));
        end
        % print spline slopes at I1bar=3 for WI1 and h(I1bar) if they exist
        if isfield(model, 'WI1_spline')
            fprintf('  Spline slope at I1bar=3 for WI1: %.4f\n', fnval(fnder(model.WI1_spline, 1), 3));
        end
        if isfield(model, 'h_spline')
            fprintf('  Spline slope at I1bar=3 for h(I1bar): %.4f\n', fnval(fnder(model.h_spline, 1), 3));
        end
        if isfield(model, 'j_spline')
            fprintf('  Spline slope at I2bar=3 for j(I2bar): %.4f\n', fnval(fnder(model.j_spline, 1), 3));
        end

        % figure with sampled invariants (I1bar, J) across modes to check coverage of spline domains
        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I1bar, detJ, 'k.-', 'MarkerSize', 11); % data
        ax = gca; ax.FontSize = Fontsize_axes; ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_1$', 'Interpreter','latex','FontSize', Fontsize); ylabel('Effective $J = det(F)$', 'Interpreter','latex','FontSize', Fontsize); grid on; title('Effective $J$ vs $\bar{I}_1$ across all modes', 'Interpreter','latex');

        % Note, in one line: set(gcf,'Units','centimeters','Position',[2 2 8 3])

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(D.lambda, Pfit.*1000, 'LineWidth',1.8, 'Color', Color_B);
        %plot(D.lambda, D.P11, 'k.', 'MarkerSize', 11); % data
        idx = 1:10:length(D.lambda); % show marker every 10 points
        plot(D.lambda, D.P11.*1000, 'LineStyle','none','Marker','o','MarkerIndices',idx,'MarkerSize',6,'MarkerEdgeColor', 'k','LineWidth',1.0);
        ax = gca; ax.FontSize = Fontsize_axes; ax.FontName = 'Times New Roman';

        % add annotation for RMSE
        %rmse = sqrt(mean((D.P11 - Pfit).^2));
        %str_rmse = sprintf('RMSE = %.4f', rmse);
        %annotation('textbox', [0.15, 0.75, 0.1, 0.1], 'String', str_rmse, 'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor','k');
        SS_res = sum(((D.P11 - Pfit).*1000).^2);          % residual sum of squares
        SS_tot = sum(((D.P11 - mean(D.P11)).*1000).^2);   % total sum of squares
        R2 = 1 - SS_res/SS_tot;
        str_R2 = sprintf('R^2 = %.4f', R2);
        annotation('textbox', [0.40, 0.83, 0.1, 0.05], 'String', str_R2, 'FitBoxToText','on', 'BackgroundColor','w', 'EdgeColor','none', 'HorizontalAlignment','center', 'FontSize',9);

        

        if strcmpi(mode,'SS')
            xlabel('$\gamma_{12}$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$P_{12}$ [kPa]', 'Interpreter','latex','FontSize', Fontsize);
            title(['Fit: ', mode, ' (fitting $P_{12}$)'], 'Interpreter','latex');      
            xlim([0 inf]); 
        else
            xlabel('$\lambda_1$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$P_{11}$ [kPa]', 'Interpreter','latex','FontSize', Fontsize);
            title(['Fit: ', mode], 'Interpreter','latex');
            if strcmpi(mode,'UT')
                ylim([0 inf])
            end
        end
        title(['Fit: ', mode], 'Interpreter','latex'); grid on;
        legend('model prediction','experimental data','Location','southeast','Box','off','Color','none');
        print(gcf, ['eq_fit_', mode, '.svg'], '-dsvg');
        savefig(gcf, ['eq_fit_', mode, '.fig']);

        % Collect P22 only for UT or UC modes
        if strcmpi(mode,'UT') || strcmpi(mode,'UC')
            lambda_all = [lambda_all; D.lambda(:)];
            modes_all = [modes_all; repmat({mode}, numel(D.lambda), 1)];
            P22_all = [P22_all; P22fit(:)];
            P11_all_data     = [P11_all_data; D.P11(:)];
            P11_all_fit      = [P11_all_fit; Pfit(:)];
        end

        %fig = figure('Color','w'); hold on;
        %fig.Units = 'centimeters';
        %fig.Position(3:4) = [8, 6];
        %plot(D.lambda, P22fit, 'LineWidth',1.8);
        %xlabel('$\lambda_1$', 'Interpreter','latex'); ylabel('Resulting $P_{22}$', 'Interpreter','latex');
        %title(['Resulting $P_{22}$ for mode: ', mode], 'Interpreter','latex'); grid on;
        %print(gcf, ['eq_fit_P22_', mode, '.svg'], '-dsvg');

        %fig = figure('Color','w'); hold on;
        %fig.Units = 'centimeters';
        %fig.Position(3:4) = [8, 6];
        %plot(D.lambda, lam2, 'LineWidth',1.8);
        %xlabel('$\lambda_1$', 'Interpreter','latex'); ylabel('Adjusted $\lambda_2$', 'Interpreter','latex');
        %title(['Adjusted $\lambda_2$ for mode: ', mode], 'Interpreter','latex'); grid on;
        %print(gcf, ['eq_fit_lambda2_', mode, '.svg'], '-dsvg');

        %fig = figure('Color','w'); hold on;
        %fig.Units = 'centimeters';
        %fig.Position(3:4) = [8, 6];
        %plot(D.lambda, detJ, 'LineWidth',1.8);
        %xlabel('$\lambda_1$', 'Interpreter','latex'); ylabel('Effective $J = det(F)$', 'Interpreter','latex');
        %title(['Effective $J$ for mode: ', mode], 'Interpreter','latex'); grid on;
        %print(gcf, ['eq_fit_J_', mode, '.svg'], '-dsvg');

    end

    %% Plot all P22 together
    % After the loop where you collect UT/UC modes:
    if ~isempty(P22_all)
        % Sort by lambda1
        [lambda_all_sorted, sort_idx] = sort(lambda_all);

        % Reorder all UT/UC vectors accordingly
        P22_all_sorted      = P22_all(sort_idx);
        P11_all_data_sorted = P11_all_data(sort_idx);
        P11_all_fit_sorted  = P11_all_fit(sort_idx);
        modes_all_sorted    = modes_all(sort_idx); % optional, for coloring/legend
    end
    fig = figure('Color','w'); hold on;
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 3];
    plot(lambda_all_sorted, P22_all_sorted.*1000, 'LineWidth',1.0, 'Color',[0.4 0.4 0.4]); % or plot(lambda_all, P22_all,'.-')
    ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
    xlabel('$\lambda_1$', 'Interpreter','latex','FontSize', Fontsize);
    ylabel('$P_{22}$ [kPa]', 'Interpreter','latex','FontSize', Fontsize);
    grid on;
    %title('$P_{22}$ for UT and UC', 'Interpreter','latex');
    print(gcf, ['eq_fit_P22_', mode, '.svg'], '-dsvg');
    savefig(gcf, ['eq_fit_P22_', mode, '.fig']);

    fig = figure('Color','w'); hold on;
    fig.Units = 'centimeters';
    fig.Position(3:4) = [8, 4];
    plot(lambda_all_sorted, P11_all_fit_sorted.*1000, 'LineWidth',1.8, 'Color', Color_B);
    %plot(lambda_all_sorted, P11_all_data_sorted, 'k.', 'MarkerSize', 11); % data
    idx = 1:10:length(D.lambda); % show marker every 10 points
    plot(lambda_all_sorted, P11_all_data_sorted.*1000, 'LineStyle','none','Marker','o','MarkerIndices',idx,'MarkerSize',6,'MarkerEdgeColor', 'k','LineWidth',1.0);
    ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
    xlabel('$\lambda_1$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$P_{11}$ [kPa]', 'Interpreter','latex','FontSize', Fontsize); 
    grid on;
    %title('Fit: UC \& UT', 'Interpreter','latex'); 
    %legend('model prediction','experimental data','Location','southeast','Box','off','Color','none');
    print(gcf, ['eq_fit_UC_UT', '.svg'], '-dsvg');
    savefig(gcf, ['eq_fit_UC_UT','.fig']);

    

    % ----------------------------
    % 2) Plot splines
    % ----------------------------
    plot_univariate(model);

   
%{
  % ----------------------------
    % 3) Postprocess adaptivity (only for univariate)
    % ----------------------------
    if strcmp(model.type, 'univariate')
        %postprocess_adaptivity_eq(x_opt, residual, J, model, prot, data, 'amr');
    else
        fprintf('Skipping adaptivity postprocessing for model.type="%s".\n', model.type);
    end

    % ----------------------------
    % 4) Compare with phen models (only for univariate)
    % ----------------------------
    if strcmp(model.type, 'univariate')
        compare_spline_phen_eq(result);
    else
        fprintf('Skipping phenomological model comparison for model.type="%s".\n', model.type);
    end 
%}


end

% ---------------------------- helpers ----------------------------
function plot_univariate(model)
    % uncoupled
    %color_A1 = '#7556AD';
    %color_A2 = '#FF8282';
    % w/ coupling term I1,J
    %color_A1 = '#9C6345';
    %color_A2 = '#FF964D';
    % w/ coupling term I2,J
    %color_A1 = '#356823';
    %color_A2 = '#DEF5B3';
    % w/ coupling term I1,I2
    %color_A1 = '#34689E';
    %color_A2 = '#B0E3EB';
    % w/ coupling term I1,I2 & without other coupling terms
    color_A1 = '#850000';
    color_A2 = '#FFC2C2';

    Fontsize = 12.1;
    Fontsize_axes = 11;

    lw = 2;
    if model.use_I1
        I1_plot     = linspace(min(model.I1bar_knots), max(model.I1bar_knots), 100);
        W_I1_plot   = fnval(model.WI1_spline, I1_plot);
        W_I1_ctrl   = fnval(model.WI1_spline, model.I1bar_knots);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I1_plot, W_I1_plot, '-', 'LineWidth', lw, 'Color', color_A1);
        plot(model.I1bar_knots, W_I1_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_1$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$\Psi(\bar{I}_1)$ [mJ/mm$^3$]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$\Psi(\bar{I}_1)$', 'Interpreter','latex');
        xlim([3 inf]); ylim([0 inf]);
        if model.use_I2
            ymax = max([W_I1_plot(:); W_I1_ctrl(:)]);   % store max y
        end
        print(gcf, 'eq_fit_WI1.svg', '-dsvg');
        savefig(gcf, ['eq_fit_WI1', '.fig']);
    end
    if model.use_I2
        I2_plot     = linspace(min(model.I2bar_knots), max(model.I2bar_knots), 100);
        W_I2_plot   = fnval(model.WI2_spline, I2_plot);
        W_I2_ctrl   = fnval(model.WI2_spline, model.I2bar_knots);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I2_plot, W_I2_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.I2bar_knots, W_I2_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_2$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$\Psi(\bar{I}_2)$ [mJ/mm$^3$]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$\Psi(\bar{I}_2)$', 'Interpreter','latex');
        xlim([3 inf]);
        if model.use_I1
            ylim([0 ymax]);
        else
            ylim([0 inf]);
        end
        print(gcf, 'eq_fit_WI2.svg', '-dsvg');
        savefig(gcf, ['eq_fit_WI2', '.fig']);
    end

    if model.use_J
        J_plot     = linspace(min(model.J_knots), max(model.J_knots), 100);
        W_J_plot   = fnval(model.WJ_spline, J_plot);
        W_J_ctrl   = fnval(model.WJ_spline, model.J_knots);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(J_plot, W_J_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.J_knots, W_J_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$J$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$\Psi(J)$ [mJ/mm$^3$]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$\Psi(J)$', 'Interpreter','latex');
        ylim([0 inf]);
        print(gcf, 'eq_fit_WJ.svg', '-dsvg');
        savefig(gcf, ['eq_fit_WJ', '.fig']);
    end

    % plot h(I1) and g(J) for multiplicative decomposition if they exist
    if isfield(model, 'h_spline') && isfield(model, 'g_spline')
        I2_h_plot = linspace(min(model.h_I1_knots), max(model.h_I1_knots), 100);
        h_plot = fnval(model.h_spline, I2_h_plot);
        h_ctrl = fnval(model.h_spline, model.h_I1_knots);
        J_g_plot = linspace(min(model.g_J_knots), max(model.g_J_knots), 100);
        g_plot = fnval(model.g_spline, J_g_plot);
        g_ctrl = fnval(model.g_spline, model.g_J_knots);
        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I2_h_plot, h_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.h_I1_knots, h_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_1$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$h(\bar{I}_1)$ [mJ/mm$^3$]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$h(\bar{I}_1)$', 'Interpreter','latex');
        xlim([3 inf]); ylim([0 inf]);
        print(gcf, 'eq_fit_h.svg', '-dsvg');
        savefig(gcf, ['eq_fit_h', '.fig']);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(J_g_plot, g_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.g_J_knots, g_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$J$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$g(J)$ [-]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$g(J)$', 'Interpreter','latex');
        ylim([0 inf]);
        print(gcf, 'eq_fit_g.svg', '-dsvg');
        savefig(gcf, ['eq_fit_g', '.fig']);
    end

    % plot j(I2) and i(J) for multiplicative decomposition if they exist
    if isfield(model, 'j_spline') && isfield(model, 'i_spline')
        I2_j_plot = linspace(min(model.j_I2_knots), max(model.j_I2_knots), 100);
        j_plot = fnval(model.j_spline, I2_j_plot);
        j_ctrl = fnval(model.j_spline, model.j_I2_knots);
        J_i_plot = linspace(min(model.i_J_knots), max(model.i_J_knots), 100);
        i_plot = fnval(model.i_spline, J_i_plot);
        i_ctrl = fnval(model.i_spline, model.i_J_knots);
        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I2_j_plot, j_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.j_I2_knots, j_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_2$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$j(\bar{I}_2)$ [mJ/mm$^3$]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$j(\bar{I}_2)$', 'Interpreter','latex');
        xlim([3 inf]); ylim([0 inf]);
        print(gcf, 'eq_fit_j.svg', '-dsvg');
        savefig(gcf, ['eq_fit_j', '.fig']);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(J_i_plot, i_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.i_J_knots, i_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$J$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$i(J)$ [-]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$i(J)$', 'Interpreter','latex');
        ylim([0 inf]);
        print(gcf, 'eq_fit_i.svg', '-dsvg');
        savefig(gcf, ['eq_fit_i', '.fig']);
    end

    % plot k(I1) and l(I2) for multiplicative decomposition if they exist
    if isfield(model, 'k_spline') && isfield(model, 'l_spline')
        I1_k_plot = linspace(min(model.k_I1_knots), max(model.k_I1_knots), 100);
        k_plot = fnval(model.k_spline, I1_k_plot);
        k_ctrl = fnval(model.k_spline, model.k_I1_knots);
        I2_l_plot = linspace(min(model.l_I2_knots), max(model.l_I2_knots), 100);
        l_plot = fnval(model.l_spline, I2_l_plot);
        l_ctrl = fnval(model.l_spline, model.l_I2_knots);
        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I1_k_plot, k_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.k_I1_knots, k_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_1$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$k(\bar{I}_1)$ [mJ/mm$^3$]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$k(\bar{I}_1)$', 'Interpreter','latex');
        xlim([3 inf]); ylim([0 inf]);
        print(gcf, 'eq_fit_k.svg', '-dsvg');
        savefig(gcf, ['eq_fit_k', '.fig']);

        fig = figure('Color','w'); hold on;
        fig.Units = 'centimeters';
        fig.Position(3:4) = [8, 6];
        plot(I2_l_plot, l_plot, '-', 'LineWidth',lw, 'Color', color_A1);
        plot(model.l_I2_knots, l_ctrl, 'ro', 'MarkerFaceColor',color_A2, 'MarkerEdgeColor',color_A1,'MarkerSize', 5);
        ax = gca; ax.FontSize = Fontsize_axes;  ax.FontName = 'Times New Roman';
        xlabel('$\bar{I}_2$', 'Interpreter','latex','FontSize', Fontsize); ylabel('$l(\bar{I}_2)$ [-]', 'Interpreter','latex','FontSize', Fontsize); grid on; title('$l(\bar{I}_2)$', 'Interpreter','latex');
        ylim([0 inf]);
        print(gcf, 'eq_fit_l.svg', '-dsvg');
        savefig(gcf, ['eq_fit_l', '.fig']);
    end
end