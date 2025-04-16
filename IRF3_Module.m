% Interferon Regulatory Factor-3 Module: Stimulated with Multiple PolyI:C Values %
% ==============================================================================

p = set_parameters;
input.Conc = [124.58, 62.29, 0, 0, 124.58, 0, 62.29, 0, 0, 124.58, 0, 0, 0, 0];
input.simtime = 300;

% Characterization of Poly_IC Input:
% =============================================================================
% Variation in Amplitude of Poly I:C Signal
POLYICvalue_one = [50 1.56 0.0488 0.00153 0.0000477 0.00000149 0.0000000466 0.00000000146 0.0000000000455 0];
POLYICvalue_two = [500 176.79 62.5 22.10 7.81 2.76 0.98 0.35 0.12 0];
POLYICvalue = [POLYICvalue_one;POLYICvalue_two];
POLYICtime= [0 30 60 90 120 150 180 210 240 2*36*60];

% Vector to store Active IRF-3 (Phosphorylated Nuclear Dimer) Values
% =========================================================================
y14_values = zeros(size(POLYICvalue,1),size(POLYICtime,1),input.simtime+1);

% Plotting POLY I:C :
% =========================================================================
% Desired Order of (PolyIC Index (i),Time Index (j)) Combination
combinations = [
    1, 1; 2, 1; 
];

% Creating External File
% ==============================================
%Output Filename
output_filename = 'extracted_values.txt';

% Initialize plot index
plot_index = 1;

%File Open
fileID = fopen(output_filename,"w");

%Column Headings
fprintf(fileID, 'PolyI:C_Index\tTime_Index\tIRF3_Values\tNFkB_Values\n');

figure;
% Loop over the defined combinations
for idx = 1:size(combinations, 1)
    i = combinations(idx, 1);
    k = 1;
    % Verify lengths match
    if length(POLYICtime) ~= length(POLYICvalue(i,:))
        error('Time and value vectors must have same length');
    end
    % Create a subplot grid of 1 row and 2 columns
    subplot(1, 2, plot_index); 
    current_POLYICvalues = POLYICvalue(i, :);
    current_POLYICtimes = POLYICtime(k, :);
        
    % Interpolate POLYIC values
     input.Polyic = interp1(POLYICtime, POLYICvalue(i,:), 0:input.simtime, 'pchip', 0);
    
    % Simulate
    [t0, y0] = simulate(p, input);
    
    % Active IRF-3 value: For Plotting all the simulations
     y14_values(i, k, :) = y0(:, 14);
    
    % Active IRF-3 Value: For Extracting to External File
    current_y14_values = y0(:, 14);
    
    %Storing IRF3 values to file
    fprintf(fileID,'%d\t%d\t',i,k);
    fprintf(fileID,'%f ',current_y14_values); % Appending IRF3 Values as a row
    fprintf(fileID,'\t'); % Placeholder for NFkB Values
    fprintf(fileID,'\n');
    
    % Plot the current combination
    plot(0:input.simtime, input.Polyic, 'LineWidth', 2.5);
    title(['Poly I:C Value for Input ' num2str(i) ', Time' num2str(k)]);
    xlabel('Time (in min)');
    ylabel('Concentration (in nM)');
    xlim([0 300]);
    ylim([0 500]);
    
    % Increment plot index
    plot_index = plot_index + 1;
end

%Close File
fclose(fileID);

%Display Message
disp(['IRF3 Values succesfully extracted and saved to',output_filename]);

% Plot Active IRF-3
% =========================================================================
figure;
plot_index = 1;
for idx = 1:size(combinations,1)
    i = combinations(idx,1);
    k = combinations(idx, 2);
    
    subplot(5, 4, plot_index); 
    plot(0:1:input.simtime, squeeze(y14_values(i,k,:)),'LineWidth',2.5);
    title(['Active IRF3 Value for Poly I:C Input ' num2str(i) ', Time Values ' num2str(k)]);
    xlabel('Time (in min)');
    ylabel('Concentration (in nM)');
    ylim([0 100])
    xlim([0 300]);
    plot_index = plot_index + 1;
end

function p = set_parameters
    % ========================Parameters===================================
    p.k1=2.4e-04; %RIG+MAVS=>C1
    p.k2=1.5e-05; %C1=>RIG+MAVS
    p.k3=2.4e-03; %C1+POLYIC=>C2
    p.k4=0.015; %C2=>C1+POLYIC
    p.k5=1.0; %TBK=>pTBK
    p.k6=0.6; %C2+TRAF=>C3
    p.k7=0.9; %C3=>C2+TRAF
    p.k8=0.001; %TBK+C3=>C4
    p.k9=1; %pTBK+C3=>C4
    p.k10=0.06; %C4=>C3+pTBK
    p.k11=0.006; %pTBK=>TBK
    p.k12=0.06; %C3=>C1+TRAF+POLYIC
    p.k13=0.06; %C4=>C1+TRAF+POLYIC+pTBK
    p.k14=0.001; %IRFc+C4=>C4+pIRFc
    p.k15=0.0005; %pIRFc+pIRFc=>pIRFcd
    p.k16=0.0003; %pIRFcd=> pIRFc+pIRFc
    p.k17=0.000000498; %pIRFn+pIRFn=>pIRFnd
    p.k18=0.006; %pIRFnd=> pIRFn+pIRFn
    p.k19=0.0085; %pIRFcd=>pIRFnd
    p.k20=0.00000021; %RIG => (Assumed)
    p.k21=0.00000022; %MAVS => (Assumed)
    p.k22=0.00000012; %TBK => (Assumed)
    p.k23=0.0000000069; %TRAF => (Assumed)
    p.k24=0.000000088; %IRFc => (Assumed)
    p.k25=0.5; % => RIG (Assumed)
    p.k26=0.5; % => MAVS (Assumed)
    p.k27=0.5; % => TBK (Assumed)
    p.k28=0.5; % => TRAF (Assumed)
    p.k29=0.5; % => IRFc (Assumed)
end

% Function to simulate
function [tt, yy] = simulate(p, input)
    y0 = input.Conc;
    T = 0:1:input.simtime;
    [t, y] = ode15s(@(t, y) f(t, y, p, input), T, y0);
    tt = t;
    yy = y;
end

% ODE function
function dydt = f(t, y, p, input)

% ======================= Input =========================================
if t >= input.simtime
    POLYIC = 0.00;
else
    POLYIC = interp1((0:input.simtime), input.Polyic, t, 'pchip');
end

% ====================== ODEs ===========================================
dydt = zeros(14, 1);
dydt(1) =-p.k1*y(1)*y(2) +p.k2*y(3)-p.k20*y(1)+p.k25;
dydt(2) =-p.k1*y(1)*y(2) + p.k2*y(3)-p.k21*y(2)+p.k26;
dydt(3) =p.k1*y(1)*y(2) - p.k2*y(3) - p.k3*y(3)*POLYIC +p.k4*y(4)+p.k12*y(6)+p.k13*y(9);
dydt(4) =p.k3*y(3)*POLYIC -p.k4*y(4)-p.k6*y(4)*y(5)+p.k7*y(6);
dydt(5) =-p.k6*y(4)*y(5) + p.k7*y(6)+p.k12*y(6)+p.k13*y(9) - p.k23*y(5)+p.k27;
dydt(6) =p.k6*y(4)*y(5)- p.k7*y(6)-p.k8*y(6)*y(7)-p.k9*y(8)*y(6)+p.k10*y(9)-p.k12*y(6);
dydt(7) =-p.k5*y(7)-p.k8*y(7)*y(6)+p.k11*y(8) -p.k22*y(7)+p.k28;
dydt(8) =p.k5*y(7)-p.k9*y(8)*y(6)+p.k13*y(9)+p.k10*y(9)-p.k11*y(8);
dydt(9) =p.k8*y(7)*y(6)+p.k9*y(8)*y(6)-p.k10*y(9)-p.k13*y(9);
dydt(10) =-p.k14*y(10)*y(9)-p.k24*y(10)+p.k29;
dydt(11) =p.k14*y(10)*y(9) - (p.k15*y(11)*y(11))*2 +(p.k16*y(12))*2;
dydt(12) =p.k15*y(11)*y(11) - p.k16*y(12) - p.k19*y(12);
dydt(13) =-(p.k17*y(13)*y(13))*2 + (p.k18*y(14))*2;
dydt(14) =p.k17*y(13)*y(13)*2 - p.k18*y(14) +p.k19*y(12);
end
