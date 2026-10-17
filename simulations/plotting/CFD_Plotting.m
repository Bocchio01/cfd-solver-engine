clc
clear variables
close all

fileNames = {
    % 'result'
    % 'output'
    % '00_40_40_100_UDS_SECOND_07_07',
    % '01_40_40_400_UDS_SECOND_07_07',
    % '02_40_40_1000_UDS_SECOND_07_07',
    % '03_40_40_1000_CDS_SECOND_07_07',
    % '04_40_40_1000_QUICK_SECOND_07_07',
    % '05_80_80_1000_UDS_SECOND_07_07',
    % '06_80_80_1000_CDS_SECOND_07_07',
    % '07_80_80_1000_QUICK_SECOND_07_07',
    % '08_129_129_1000_UDS_SECOND_07_07',
    '09_129_129_1000_CDS_SECOND_07_07',
    % '10_129_129_1000_QUICK_SECOND_07_07',
    % '11_129_129_1000_UDS_FOURTH_07_07',
    % '12_129_129_1000_CDS_FOURTH_07_07',
    % '13_129_129_1000_QUICK_FOURTH_07_07',
    };

plots = true * [1 0];
% latex_img_path = 'C:/Users/Bocchio/Documents/GitHub/University_Programming_Classes/07 - ME663 Computational Fluid Dynamics/Report assignment 1/img';


%% Loading data

Ghia_Solution_U = table2struct(readtable('Ghia_Solution_U.xlsx'));
Ghia_Solution_V = table2struct(readtable('Ghia_Solution_V.xlsx'));

if(plots(1) == true)
    for file_idx = 1:length(fileNames)
        close all

        file_path = ['../output/' fileNames{file_idx} '.dat'];
        if (exist(file_path, 'file') ~= 2)
            disp(['File ' fileNames{file_idx} ' doesnt exists'])
            continue;
        end

        result = readtable(file_path, 'NumHeaderLines', 7);

        x = result.Var1;
        y = result.Var2;
        u = result.Var3;
        v = result.Var4;
        p = result.Var5;

        [x, sort_idx] = sort(x);
        y = y(sort_idx);
        u = u(sort_idx);
        v = v(sort_idx);
        p = p(sort_idx);

        % Create a grid
        [X, Y] = meshgrid(unique(x), unique(y));

        % Reshape states to match the grid
        U = reshape(u, [numel(unique(y)), numel(unique(x))]);
        V = reshape(v, [numel(unique(y)), numel(unique(x))]);
        P = reshape(p, [numel(unique(y)), numel(unique(x))]);

        axis_limit = [
            min(X, [], "all")
            max(X, [], "all")
            min(Y, [], "all")
            max(Y, [], "all")
            ];

        clear x y u v p
        clear sort_idx


        %% Plot states distribution
        % 
        % figure_state_distribution = figure( ...
        %     'Name', 'States distribution', ...
        %     'NumberTitle', 'off', ...
        %     'Position', [100, 100, 1500, 600]);
        % t = tiledlayout(1, 3);
        % title(t, strrep(fileNames{file_idx}, '_', '-'))
        % 
        % states_name = ["P", "U", "V"];
        % states_matrix = {P, U, V};
        % 
        % for state_idx = 1:length(states_matrix)
        % 
        %     nexttile
        %     hold on
        %     grid on
        % 
        %     imagesc(X(1,:), Y(:,1), states_matrix{state_idx});
        % 
        %     colorbar;
        %     title(states_name(state_idx));
        %     xlabel('X');
        %     ylabel('Y');
        %     axis equal;
        %     axis(axis_limit);
        % 
        % end
        % 
        % clear state_idx states_name states_matrix


        %% Plot velocity by streamlines

        figure_velocity = figure( ...
            'Name', 'Velocity visualization', ...
            'NumberTitle', 'off', ...
            'Position', [100, 100, 1500, 600]);
        t = tiledlayout(2, 2);
        title(t, strrep(fileNames{file_idx}, '_', '-'))

        % Streamlines
        plot_streamline = nexttile(1, [2, 1]);
        hold on
        grid on

        contourf(X, Y, sqrt(U.^2 + V.^2), 30);
        colormap(jet);
        colorbar
        shading interp;


        % Subsample for cleaner quiver plot        
        stream_obj = streamslice(X, Y, U, V, 3);
        set(stream_obj, ...
            'color', [0.0 0.0 0.0], ...
            'linewidth', 1);

        title('Velocity streamlines');
        xlabel('X');
        ylabel('Y'  );
        axis equal;
        axis(axis_limit);

        clear verts averts

        %% Plots for comparison with Ghia's solution

        % figure_ghia_comparison = figure( ...
        %     'Name', 'Ghia solution comparison', ...
        %     'NumberTitle', 'off', ...
        %     'Position', [100, 100, 1500, 600]);
        % t = tiledlayout(1, 2);
        % title(t, strrep(fileNames{file_idx}, '_', '-'))

        % Comparison of u component with Ghia's solution
        nexttile
        hold on;
        grid on;

        plot(U(:, ceil(size(X, 1)/2)), Y(:,1), 'k', 'LineWidth', 2);
        plot([Ghia_Solution_U.Re100], [Ghia_Solution_U.y], 'go');
        plot([Ghia_Solution_U.Re400], [Ghia_Solution_U.y], 'bo');
        plot([Ghia_Solution_U.Re1000], [Ghia_Solution_U.y], 'ro');

        title("u velocity @x/LX=0.5")
        xlabel("u");
        xlabel("y");
        axis tight;
        axis equal;
        legend("u velocity @x/LX=0.5", ...
            "Ghia (Re = 100)", ...
            "Ghia (Re = 400)", ...
            "Ghia (Re = 1000)", ...
            "Location", "best")

        % Comparison of v component with Ghia solution
        nexttile
        hold on;
        grid on;

        plot(X(1,:), V(ceil(size(Y, 1)/2), :), 'k', 'LineWidth', 2);
        plot([Ghia_Solution_V.x], [Ghia_Solution_V.Re100], 'go');
        plot([Ghia_Solution_V.x], [Ghia_Solution_V.Re400], 'bo');
        plot([Ghia_Solution_V.x], [Ghia_Solution_V.Re1000], 'ro');

        title("v velocity @y/LY=0.5")
        xlabel("v");
        xlabel("x");
        axis tight;
        axis equal;
        legend("v velocity @y/LY=0.5", ...
            "Ghia (Re = 100)", ...
            "Ghia (Re = 400)", ...
            "Ghia (Re = 1000)", ...
            "Location", "best")

        disp(['Current visalized: ' fileNames{file_idx}])
        if (length(fileNames) > 1)
            disp('Press a key for the next result')
            pause;
        end

        if(exist('latex_img_path', 'var') == 1)
            saveas(figure_state_distribution, [latex_img_path '/states/' strrep(fileNames{file_idx}, '_', '-') '.png']);
            saveas(figure_velocity, [latex_img_path '/velocity/' strrep(fileNames{file_idx}, '_', '-') '.png']);
            saveas(plot_streamline, [latex_img_path '/streamline/' strrep(fileNames{file_idx}, '_', '-') '.png']);
            saveas(figure_ghia_comparison, [latex_img_path '/ghia_comparison/' strrep(fileNames{file_idx}, '_', '-') '.png']);
        end
    end
end

clear t plot_streamline


%% Residual and CPU analysis

if(plots(2) == true)

    figure_residual_analysis = figure( ...
        'Name', 'Residual analysis', ...
        'NumberTitle', 'off', ...
        'Position', [100, 100, 1500, 600]);

    nexttile
    hold on
    grid on

    legends = {};
    CPU_times = [];
    for file_idx = 1:length(fileNames)

        file_path = ['../output/' fileNames{file_idx} '.dat'];
        if (exist(file_path, 'file') ~= 2)
            disp(['File ' fileNames{file_idx} ' doesnt exists'])
            continue;
        end

        fileID = fopen(file_path, 'r');

        input_json_name = fgetl(fileID);

        RE = sscanf(fgetl(fileID), 'RE=%f');
        iterations = sscanf(fgetl(fileID), 'ITERATIONS=%d');

        CPU_times = [CPU_times sscanf(fgetl(fileID), 'CPU_TIME=%f')];

        residuals_str = strsplit(fgetl(fileID), ' ');
        % Skip 'RESIDUALS:'
        residuals = cellfun(@str2double, residuals_str(2:end));

        fclose(fileID);

        it = [0:50:iterations iterations];
        plot(it(1:length(residuals)), residuals, "LineWidth", 2);

        legends = [legends fileNames{file_idx}];

    end

    yline(1e-5,'--', 'Tolerance')

    title("Residual analysis");
    xlabel("Iteration #")
    ylabel("Residual")
    legend(legends, 'Location','best', "Interpreter", "none")



    %%
    figure_CPU_analysis = figure( ...
        'Name', 'CPU analysis', ...
        'NumberTitle', 'off', ...
        'Position', [100, 100, 800, 600]);

    nexttile
    hold on
    grid on

    bar(cellfun(@(x) x(1:2), legends, 'UniformOutput', false), CPU_times);
    % legend(legends,'Location','best', "Interpreter", "none")
    xlabel("Simulation ID");
    ylabel('CPU Time [s]');


    clear residuals_str

    if(exist('latex_img_path', 'var') == 1)
        saveas(figure_residual_analysis, [latex_img_path '/residual.png']);
        saveas(figure_CPU_analysis, [latex_img_path '/CPU_time.png']);
    end

end






