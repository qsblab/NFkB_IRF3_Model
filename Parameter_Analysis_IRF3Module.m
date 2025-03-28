% =============================================================================
% Date: 25 March 2025
% Purpose: This MATLAB script performs analysis to estimate active IRF-3 
% levels inside the nucleus for varying model parameters.
% =============================================================================

% Initial concentration values for all species in the IRF-3 module
input.Conc = [124.58, 62.29, 0, 0, 124.58, 0, 62.29, 0, 0, 124.58, 0, 0, 0, 0]; 

% Poly I:C input values at different time points (in nM)
POLYICvalues = [500 176.79 62.5 22.10 7.81 2.76 0.98 0.35 0.12 0]; 

% Time points corresponding to the Poly I:C input values (in minutes)
POLYICtime = [0 30 60 90 120 150 180 210 240 2*36*60]; 

% Sets simulation time (in minutes)
input.simtime = 500;

% Interpolate Poly I:C values over time for simulation using a piecewise cubic Hermite interpolant (pchip)
input.Polyic = interp1(POLYICtime, POLYICvalues, 0:input.simtime, 'pchip');

% List of parameter names used in the model
paramNames = {'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'k7', 'k8', 'k9', 'k10', 'k11', 'k12', 'k13', 'k14', 'k15', 'k16', 'k17', 'k18', 'k19', 'k20', 'k21', 'k22', 'k23', 'k24'};

% Set initial values for parameters (from the function 'set_parameters')
pOriginal = set_parameters;

% Loop for the first 3 figures (for visualization of different parameter sets)
for figIdx = 1:3
    figure;  % Creates a new figure for each set of parameters
    for i = 1:8
        paramIdx = (figIdx-1)*8 + i; % Calculates the parameter index to plot
        if paramIdx > length(paramNames)
            break; 
        end
        
        paramName = paramNames{paramIdx}; 
        
        % Creates copies of the original parameter set for variations
        pSmaller = pOriginal;
        pLarger = pOriginal;
        pSmallest = pOriginal;
        pLargest = pOriginal;
        
        % Modifies the selected parameter by a factor of 100 and 1000 smaller and larger
        pSmaller.(paramName) = pOriginal.(paramName) / 100;
        pLarger.(paramName) = pOriginal.(paramName) * 100;
        pSmallest.(paramName) = pOriginal.(paramName) / 1000;
        pLargest.(paramName) = pOriginal.(paramName) * 1000;
        
        % Simulates the system with different parameter sets
        [t0, y0] = simulate(pOriginal, input); % Original parameter set
        [tSmaller, ySmaller] = simulate(pSmaller, input); % 100 times smaller
        [tLarger, yLarger] = simulate(pLarger, input); % 100 times larger
        [tSmallest, ySmallest] = simulate(pSmallest, input); % 1000 times smaller
        [tLargest, yLargest] = simulate(pLargest, input); % 1000 times larger
        
        % Creates a subplot for visualizing the results
        subplot(3, 3, i); 
        plot(t0, y0(:, 14), 'k', 'LineWidth', 2.5); hold on; % Plot original simulation
        plot(tSmaller, ySmaller(:, 14), 'b', 'LineWidth', 2.5); % Plot 100 times smaller simulation
        plot(tLarger, yLarger(:, 14), 'r', 'LineWidth', 2.5); % Plot 100 times larger simulation
        plot(tSmallest, ySmallest(:, 14), 'b', 'LineWidth', 2.5, 'LineStyle', '-.'); % Plot 1000 times smaller simulation
        plot(tLargest, yLargest(:, 14), 'r', 'LineWidth', 2.5, 'LineStyle', '-.'); % Plot 1000 times larger simulation
        
        % Adds titles, labels, and legends to the plot
        title(['Impact of ', paramName, ' Variations']);
        xlabel('Time (in min)', 'FontWeight', 'bold', 'FontSize', 12);
        ylabel('IRF-3 (in nM)', 'FontWeight', 'bold', 'FontSize', 12);
        ylim([0 80]); % Set the y-axis limit for IRF-3 concentration
        legend('Original', '100x Smaller', '100x Larger', '1000x Smaller', '1000x Larger', 'Location', 'best', 'FontWeight', 'bold', 'FontSize', 10);
    end
    
    % Creates an empty subplot for displaying the legend on the third figure
    if figIdx == 3
        subplot(3, 3, 9);
        axis off;
    end
end

% Function to set model parameters
function p = set_parameters
    % Assigns initial values to all the parameters for the model
    p.k1 = 2.4e-04;                           
    p.k2 = 1.5e-05;                           
    p.k3 = 2.4e-03;                           
    p.k4 = 0.015;                             
    p.k5 = 1.0;                               
    p.k6 = 0.6;                               
    p.k7 = 0.9;                               
    p.k8 = 0.001;                             
    p.k9 = 1;                                 
    p.k10 = 0.06;                             
    p.k11 = 0.006;                            
    p.k12 = 0.06;                             
    p.k13 = 0.06;                             
    p.k14 = 0.001;                            
    p.k15 = 0.0005;                           
    p.k16 = 0.0003;                           
    p.k17 = 0.000000498;                      
    p.k18 = 0.006;                            
    p.k19 = 0.0085;                           
    p.k20 = 0.00000021;                       
    p.k21 = 0.00000022;                       
    p.k22 = 0.00000012;                       
    p.k23 = 0.0000000069;                     
    p.k24 = 0.000000088;                      
end

% Function to simulate the model dynamics using ODE solver
function [tt, yy] = simulate(p, input)
    y0 = input.Conc;    % Initial concentrations
    T = 0:0.1:input.simtime; % Time vector for simulation
    [t, y] = ode15s(@(t, y) f(t, y, p, input), T, y0); % Solve ODEs using ode15s
    tt = t; % Time output
    yy = y; % Concentration output
end

% ODE function defining the model dynamics
function dydt = f(t, y, p, input)
    % Interpolate Poly I:C concentration at current time t
    if t >= input.simtime
        POLYIC = 0.00; % If time exceeds simulation time, set Poly I:C to 0
    else
        POLYIC = interp1((0:input.simtime), input.Polyic, t, 'pchip'); % Interpolate Poly I:C value at time t
    end
    
    % System of ODEs representing the dynamics of the model (rate of change for each species)
    dydt = zeros(14, 1); % Initializes the derivative vector
    dydt(1) = -p.k1 * y(1) * y(2) + p.k2 * y(3) - p.k20 * y(1);
    dydt(2) = -p.k1 * y(1) * y(2) + p.k2 * y(3) - p.k21 * y(2);
    dydt(3) = p.k1 * y(1) * y(2) - p.k2 * y(3) - p.k3 * y(3) * POLYIC + p.k4 * y(4) + p.k12 * y(6) + p.k13 * y(9);
    dydt(4) = p.k3 * y(3) * POLYIC - p.k4 * y(4) - p.k6 * y(4) * y(5) + p.k7 * y(6);
    dydt(5) = -p.k6 * y(4) * y(5) + p.k7 * y(6) + p.k12 * y(6) + p.k13 * y(9) - p.k23 * y(5);
    dydt(6) = p.k6 * y(4) * y(5) - p.k7 * y(6) - p.k8 * y(6) * y(7) - p.k9 * y(8) * y(6) + p.k10 * y(9) - p.k12 * y(6);
    dydt(7) = -p.k5 * y(7) - p.k8 * y(7) * y(6) + p.k11 * y(8) - p.k22 * y(7);
    dydt(8) = p.k5 * y(7) - p.k9 * y(8) * y(6) + p.k13 * y(9) + p.k10 * y(9) - p.k11 * y(8);
    dydt(9) = p.k8 * y(7) * y(6) + p.k9 * y(8) * y(6) - p.k10 * y(9) - p.k13 * y(9);
    dydt(10) = -p.k14 * y(10) * y(9) - p.k24 * y(10);
    dydt(11) = p.k14 * y(10) * y(9) - (p.k15 * y(11) * y(11)) * 2 + (p.k16 * y(12)) * 2;
    dydt(12) = p.k15 * y(11) * y(11) - p.k16 * y(12) - p.k19 * y(12);
    dydt(13) = -(p.k17 * y(13) * y(13)) * 2 + (p.k18 * y(14)) * 2;
    dydt(14) = p.k17 * y(13) * y(13) * 2 - p.k18 * y(14) + p.k19 * y(12);
end
