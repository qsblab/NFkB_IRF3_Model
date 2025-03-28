% Nuclear Factor kappa B Module: Stimulated with Multiple PolyI:C Values %
% =========================================================================

% Set simulation parameters (function call)
p = set_parameters;

% Equilibriation phase (initial setup phase for the system)
% ===================
input.Conc= zeros(1,49); % Initial concentrations (set to zeros)
input.simtime=30000; % Total simulation time for equilibrium phase (30,000 minutes)
input.Phase=1;% Phase 1: Equilibrium phase
input.NikPhase=1; % Initial phase for NIK signaling
[t0, y0] = simulate(p,input);% Runs simulation for the equilibration phase

% Simulation phase (running the actual simulation)
% ===================
input.simtime=300;% Shorter simulation time for the active phase (300 minutes)
input.Phase=2;% Switch to Phase 2: Active phase
input.Conc=y0(end,:);% Used the final concentrations from the equilibration phase
input.NikPhase=2; % Sets phase for NIK signaling during the simulation
input.Conc(27)=0;% Sets an Active NFkB concentration (index 27) to zero
NKvalues =[0.0000 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];% No canonical signaling for NK
NKtime = [0 30 60 180 300 480 960 1440 1440*2.5];% Time points for NIK signaling
input.Nik = interp1(NKtime, NKvalues, 0:input.simtime,'pchip');% Interpolates NIK values

% Characterization of Poly_IC Input:
% =========================================================================
% Defining the variations in PolyI:C signal amplitude
POLYICvalue_one = [50 1.56 0.0488 0.00153 0.0000477 0.00000149 0.0000000466 0.00000000146 0.0000000000455 0];
POLYICvalue_two = [500 176.79 62.5 22.10 7.81 2.76 0.98 0.35 0.12 0];
POLYICvalue = [POLYICvalue_one;POLYICvalue_two]; % Matrix of different PolyI:C signal variations

% Defining time points for PolyI:C signal variation
POLYICtime = [0 30 60 90 120 150 180 210 240 2*36*60];

% Vector to store Active NFkB (Nuclear Concentration) Values
% =========================================================================
Y27_values = zeros(size(POLYICvalue,1),size(POLYICtime,1),input.simtime+1);

% Plotting POLY I:C:
% =========================================================================
% Defining combinations of PolyIC signal and time to simulate
combinations = [
    1, 1;
    2, 1; 
];

% Creating External File
% ==============================================
% Output Filename
output_filename = 'extracted_values.txt';

% Reading existing data
fileID = fopen(output_filename,'r');
existingData = textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);

% Opening File for Rewriting
fileID = fopen(output_filename,'w');

% Header
fprintf(fileID, '%s\n', existingData{1}{1});

% Create a figure for PolyI:C concentration
figure;
plot_index = 1;

% Loop through combinations to simulate and plot
for idx = 1:size(combinations, 1)
    i = combinations(idx, 1);% PolyIC signal variation index
    k = 1;% Time point index
    
    % Create a subplot grid of 1 row and 2 columns
    subplot(1, 2, plot_index); 
    current_POLYICvalues = POLYICvalue(i, :);% Get current PolyIC values
    current_POLYICtimes = POLYICtime(k, :);% Get current PolyIC time points
    
    % Interpolate PolyIC values for the simulation
    input.Polyic = interp1(current_POLYICtimes, current_POLYICvalues, 0:input.simtime, 'pchip');
    
    % Runs the simulation
    [T1, Y1] = simulate(p,input);
    
    % Stores the NFkB values for the current simulation 
    Y27_values(i,k,:) = Y1(:, 27);
    
     % Active NFkB value: For Extracting to the external file
    current_Y27_values = Y1(:, 27);
    
    % Find the corresponding row in existing data
    rowToEdit = existingData{1}{idx+1}; %Skips Header Row
    
    %Split the row into components
    rowParts = strsplit(rowToEdit, '\t');
    
    %Append NFkB Value to the respective line in the file
    rowParts{4} = sprintf('%f ', current_Y27_values);
    
    %Writing Updated Row to file
    fprintf(fileID, '%s\t%s\t%s\t%s\n', rowParts{1}, rowParts{2}, rowParts{3}, rowParts{4});

    % Plots the PolyI:C concentration curve
    plot(0:input.simtime, input.Polyic, 'LineWidth', 2.5);
    title(['Value Iteration ' num2str(i) ', Time Iteration ' num2str(k)]);
    xlabel('Time (in min)');
    ylabel('Poly I:C Concentration (in nM)');
    xlim([0 300]);
    ylim([0 500]);
    
    % Increments plot index for the next subplot
    plot_index = plot_index + 1;
end

% Close the file
fclose(fileID);

% Display a message indicating the operation was succesfule
disp(['NFkB values successfully appended to',output_filename]);

% Plotting Active NFkB (Nuclear Concentration of RelAp52):
% =========================================================================
figure;
plot_index = 1;

% Loops through combinations to plot NFkB values
for idx = 1:size(combinations, 1)
    i = combinations(idx,1);% PolyIC signal variation index
    k = 1;% Time point index

    subplot(1, 2, plot_index);% Create subplot
    
    % Plots the Active NFkB values over time
    plot(0:1:input.simtime, squeeze(Y27_values(i,k,:)),'LineWidth',2.5);
    title(['Active NFkB Value Iteration ' num2str(i) ', Time Iteration ' num2str(k)]);
    xlabel('Time (in min)');
    ylabel('Concentration (in nM)');
    ylim([0 300])
    xlim([0 300]);

    % Increments plot index for the next subplot
    plot_index = plot_index + 1;
end

% Functions to set simulation parameters
function p = set_parameters
    % =================== Define model parameters =====================
    % Kinetic parameters for various molecules and processes
    p.k01=1.3e-2; % tIκBα (basal) 
    p.k02=6.6e-4; % tIκBβ (basal) 
    p.k03=1.2e-4; % tIκBε (basal) 
    p.k04=3.0e-5; % tRelA (basal)
    p.k05=6.0e-5;	%	tp50 (basal)
    p.k06=8.0e-7;	%	tRelB (basal)
    p.k07=3.0e-5;	%	tp100 (basal)
    p.k10=12;	%	mRNA => mRNA + protein
    p.k10a=12;	%	mRNA => mRNA + protein for p100
    p.k11=25;	%	=> tIκBα (RelAp50-induced)
    p.k12=150;	%	Hill Kd (RelAp50-induced)
    p.k13=125;	%	=>tIκBε (A50-induced, 37 min. delay)
    p.k14=150;	%	Hill Kd (RelAp50-induced)
    p.k15=1.2e-3;	%	p100 + p100 => IκBδ
    p.k16=1.2e-2;	%	IκBδ => p100 + p100
    p.k17=6.0e-2;	%	IκBα(c) => IκBα(n)
    p.k18=9.0e-3;	%	IκBβ(c) => IκBβ(n)
    p.k18a=4.5e-2;	%	IκBε(c) => IκBε(n), IκBδ(c) => IκBδ(n)
    p.k19=1.2e-2;	%	IκB[α/β/ε/δ](n) => IκB[α/β/ε/δ](c)
    p.k20=4.4e-2;	%	tIκBα =>
    p.k21=2.9e-3;	%	tIκBβ =>
    p.k22=3.8e-3;	%	tIκBε =>
    p.k23=0.12;	%	IκBα =>
    p.k24=0.18;	%	IκBβ => , IκBε =>
    p.k25=1.1e-2;	%	IκBδ =>
    p.k26=2.4e-4;	%	IκB[α/β/ε/δ]-NFκB => NFκB
    p.k27=1.4e-3;	%	IκBα => (NEMO-mediated), IκBα-NFκB => NFκB (NEMO-mediated)
    p.k28=4.5e-4;	%	IκBβ => (NEMO-mediated), IκBβ-NFκB => NFκB (NEMO-mediated)
    p.k29=9.0e-4;	%	IκBε => (NEMO-mediated), IκBε-NFκB => NFκB (NEMO-mediated)
    p.k30=1.2;	%	IκBδ => (NIK-mediated), IκBδ-NFκB => NFκB (NIK-mediated)
    p.k31=100;	%	IκBδ => (NIK-mediated, Km)
    p.k32=5.0e-2;	%	p100 => p52 (NIK-mediated)
    p.k33=10;	%	p100 => p52 (NIK-mediated, p100 Km)
    p.k34=100;	%	tRelB (A50-induced, 1 hr delay)
    p.k35=50;	%	Hill Kd (RelAp50-induced)
    p.k36=50;	%	tp100 (A50-induced, 1 hr delay)
    p.k37=50;	%	Hill Kd (RelAp50-induced)
    p.k38=2.9e-3;	%	tRelA =>, tp50 =>, tRelB =>
    p.k39=1.9e-3;	%	tp100 =>
    p.k40=5.7e-3;	%	RelA =>, p50 =>, RelB =>, p100 =>, p52 =>
    p.k41=9.6e-4;	%	RelA + p50 => RelAp50, RelB + p52 => RelBp52
    p.k42=9.6e-5;	%	RelB + p50 => RelBp50
    p.k43=1.3e-4;	%	RelB + p100 => RelBp100
    p.k44=1.4e-3;	%	RelAp50 => RelA + p50, RelBp52 => RelB + p52
    p.k45=2.2e-3;	%	RelBp50 => RelB + p50
    p.k46=1.1e-3;	%	RelBp100 => RelB + p100
    p.k47=5.4;	%	RelANFkB(c) => RelANFkB(n), RelBp50(c) => RelBp50(n), RelBp100(c) => RelBp100(n)
    p.k48=4.8e-3;	%	RelBp100(c) => RelBp100(n), NFκB(n) => NFκB(c)
    p.k49=2.4e-4;	%	RelAp50 =>, RelBp50 =>, RelBp100 =>, IκB[α/β/ε/δ]-NFκB => IκB[α/β/ε/δ]
    p.k50=9.6e-4;	%	RelBp52 =>
    p.k51=0.2;	%	IκBα/IκBβ/IκBε/IκBδ + RelAp50 => IκB-RelAp50
    p.k52=1.7e-2;	%	IκBα/IκBε + RelBp50 => IκBα/IκBε-RelBp50
    p.k53=1.2e-1;	%	IκBδ + RelBp50 => IκBδ-RelBp50
    p.k53a=0.5;	%	IκBδ + RelBp52 => IκBδ-RelBp50
    p.k54=8.4e-3;	%	IκBα/IκBε/IκBδ-RelAp50 => IκBα/IκBε/IκBδ + RelAp50
    p.k55=3.4e-2;	%	IκBβ-RelAp50 => IκBβ + RelAp50
    p.k56=3.2e-2;	%	IκBα/IκBε-RelBp50 => IκBα/IκBε + RelBp50
    p.k57=1.4e-2;	%	IκBδ-RelBp50 => IκBδ + RelBp50
    p.k58=0.28;	%	IκBα-NFκB(c) => IκBα-NFκB(n), IκBβ-NFκB(c) => IκBβ-NFκB(n)
    p.k58a=0.028;	%	IκBδ-NFκB(c) => IκBδ-NFκB(n)
    p.k59=0.14;	%	IκBδ-NFκB(c) => IκBδ-NFκB(n)
    p.k60=0.84;	%	IκBα-NFκB(n) => IκBα-NFκB(c)
    p.k61=0.42;	%	IκBβ-NFκB(n) => IκBβ-NFκB(c), IκBε-NFκB(n) => IκBε-NFκB(c), IκBδ-NFκB(n) => IκBδ-NFκB(c)
end

% Functions to simulate the system
function [tt,yy] = simulate(p,input)
    y0 = input.Conc; % Initial concentrations from the input
    % Time vector for simulation
    T = 0:1:input.simtime;
    % Solve the ODE system using ode15s 
    [t,y] = ode15s(@f,T,y0,[],p,input);
    tt = t; % Output time vector
    yy = y; % Output concentration data
end

% ODE function for the system of equations
function [dydt]=f(t,y,p,input)
persistent itxn_eps1 itxn_relB1 itxn_p1001; 
persistent itxn_eps itxn_relB itxn_p100;

 % Inputs for the model
if input.Phase == 1 
    POLYIC = 0.00;
    if input.NikPhase== 1
        NIK = 0.00;
    elseif input.NikPhase== 2 
        if t >= input.simtime
            NIK = input.Nik(end);
        else	
            start = floor(t)+1;stop = start + 1;
            NIK = ((input.Nik(stop) - input.Nik(start)) * (t-floor(t)) + input.Nik(start)); 
        end
    end
else
    % PolyIC and NIK signaling for Phase 2
    if t >= input.simtime
        POLYIC = 0.00;
        NIK = input.Nik(end);
    else
        start = floor(t)+1;stop = start + 1;
        POLYIC = interp1((0:input.simtime), input.Polyic, t, 'pchip');
        NIK = ((input.Nik(stop) - input.Nik(start)) * (t-floor(t)) + input.Nik(start));
    end
end
%=================== Delay and Interaction Calculations ==================================
if input.Phase == 1
    itxn_eps = p.k03*(1+(p.k13*y(27)/p.k14))/(1+(y(27)/p.k14)); 
    itxn_p100 = p.k07*(1+(p.k36*y(27)/p.k37))/(1+(y(27)/p.k37)); 
    itxn_relB = p.k06*(1+(p.k34*y(27)/p.k35))/(1+(y(27)/p.k35));
else
    if t==0
        itxn_eps1 = p.k03*(1+(p.k13*y(27)/p.k14))/(1+(y(27)/p.k14)); 
    elseif t < 37
        itxn_eps=itxn_eps1; 
    elseif t >= 37 && t < 37 + 9
        c = (t - 37); 
        r = interp1(0:1:9 , .1:.1:1, c ,'pchip'); 
        itxn_eps = p.k03*(1+(r*p.k13*y(27)/p.k14))/(1+(r*y(27)/p.k14)); 
    else
        itxn_eps = p.k03*(1+(p.k13*y(27)/p.k14))/(1+(y(27)/p.k14)); 
    end

    if t==0
        itxn_p1001 = p.k07*(1+(p.k36*y(27)/p.k37))/(1+(y(27)/p.k37)); 
        itxn_relB1 = p.k06*(1+(p.k34*y(27)/p.k35))/(1+(y(27)/p.k35)); 
    elseif t < 60
        itxn_p100=itxn_p1001; itxn_relB=itxn_relB1;
    elseif t >= 60 && t < 60 + 9
    c = (t - 60); r = interp1(0:1:9 , .1:.1:1, c ,'pchip'); 
    itxn_p100 = p.k07*(1+(r*p.k36*y(27)/p.k37))/(1+(y(27)/p.k37)); 
    itxn_relB = p.k06*(1+(r*p.k34*y(27)/p.k35))/(1+(y(27)/p.k35)); 
    else
    itxn_p100 = p.k07*(1+(p.k36*y(27)/p.k37))/(1+(y(27)/p.k37)); 
    itxn_relB = p.k06*(1+(p.k34*y(27)/p.k35))/(1+(y(27)/p.k35)); 
    end
end
%==============================Model States=================================
tIkBa = y(1);	
tIkBb = y(2);	
tIkBe = y(3);	
tRelA = y(4);
tRelB = y(5);	
tp50 = y(6);	
tp100 = y(7);	
IkBa = y(8);
IkBan = y(9);	
IkBb = y(10);	
IkBbn = y(11);	
IkBe = y(12);
IkBen = y(13);	
IkBd = y(14);	
IkBdn = y(15);	
RelA = y(16);
RelAn = y(17);	
RelB = y(18);	
RelBn = y(19);	
p50 = y(20);
p50n = y(21);	
p100 = y(22);	
p100n = y(23);	
p52 = y(24);
p52n = y(25);	
RelAp50 = y(26);	
RelAp50n = y(27);	
RelBp52 = y(28);
RelBp52n = y(29);	
RelBp100 = y(30);	
RelBp100n = y(31);	
IkBaRelAp50= y(32);
IkBaRelAp50n= y(33);	
IkBbRelAp50= y(34);	
IkBbRelAp50n= y(35);	
IkBeRelAp50= y(36);
IkBeRelAp50n= y(37);	
IkBdRelAp50= y(38);	
IkBdRelAp50n= y(39);	
RelBp50=y(40);
RelBp50n=y(41);	
IkBaRelBp50= y(42);	
IkBaRelBp50n= y(43);	
IkBeRelBp50= y(44);
IkBeRelBp50n= y(45);
IkBdRelBp52n= y(49);	
IkBdRelBp50= y(46);	
IkBdRelBp50n= y(47);	
IkBdRelBp52= y(48);
%==============================System of Ordinary Differential Equations ===================================
dydt(1,1)=	p.k01*(1+(p.k11*RelAp50n/p.k12))/(1+(RelAp50n/p.k12))-p.k20*tIkBa;	%tIkBa
dydt(2,1)=	p.k02-p.k21*tIkBb;	%tIkBb
dydt(3,1)= itxn_eps-p.k22*tIkBe;	%tIkBe 
dydt(4,1)= p.k04-p.k38*tRelA; %tRelA
dydt(5,1)= itxn_relB-p.k38*tRelB; %tRelB
dydt(6,1)= p.k05-p.k38*tp50; %tp50
dydt(7,1)= itxn_p100-p.k39*tp100; %tp100
dydt(8,1)= tIkBa*p.k10-p.k17*IkBa+p.k19*IkBan-p.k23*IkBa-p.k27*POLYIC*IkBa-p.k51*IkBa*RelAp50...
+p.k54*IkBaRelAp50-p.k52*IkBa*RelBp50+p.k56*IkBaRelBp50+p.k49*IkBaRelAp50...
+p.k49*IkBaRelBp50;%IkBa
dydt(9,1)= p.k17*IkBa-p.k19*IkBan-p.k23*IkBan-p.k51*IkBan*RelAp50n+p.k54*IkBaRelAp50n...
-p.k52*IkBan*RelBp50n+p.k56*IkBaRelBp50n+p.k49*IkBaRelAp50n+p.k49*IkBaRelBp50n;%IkBan
dydt(10,1)= tIkBb*p.k10-p.k18*IkBb+p.k19*IkBbn-p.k24*IkBb-p.k28*POLYIC*IkBb-p.k51*IkBb*RelAp50...
+p.k55*IkBbRelAp50+p.k49*IkBbRelAp50;	%IkBb
dydt(11,1)=  p.k18*IkBb-p.k19*IkBbn-p.k51*IkBbn*RelAp50n+p.k55*IkBbRelAp50n+p.k49*IkBbRelAp50n...
-p.k24*IkBbn; %IkBbn
dydt(12,1)=  tIkBe*p.k10-p.k18a*IkBe+p.k19*IkBen-p.k24*IkBe-p.k29*POLYIC*IkBe-p.k51*IkBe*RelAp50...
+p.k54*IkBeRelAp50-p.k52*IkBe*RelBp50+p.k56*IkBeRelBp50+p.k49*IkBeRelAp50...
+p.k49*IkBeRelBp50;%IkBe
dydt(13,1)=  p.k18a*IkBe-p.k19*IkBen-p.k51*IkBen*RelAp50n+p.k54*IkBeRelAp50n-p.k52*IkBen*RelBp50n...
+p.k56*IkBeRelBp50n+p.k49*IkBeRelAp50n+p.k49*IkBeRelBp50n-p.k24*IkBen;%IkBen
dydt(14,1)=  p.k15*p100*p100-p.k16*IkBd-p.k25*IkBd-p.k18a*IkBd+p.k19*IkBdn-NIK*p.k30*IkBd/(p.k31+IkBd)...
-p.k51*IkBd*RelAp50+p.k54*IkBdRelAp50-p.k53*IkBd*RelBp50-p.k53a*IkBd*RelBp52...
+p.k57*IkBdRelBp50+p.k57*IkBdRelBp52+p.k49*IkBdRelAp50+p.k49*IkBdRelBp50+p.k49*IkBdRelBp52;%IkBd 
dydt(15,1)= p.k15*p100n*p100-p.k16*IkBdn+p.k18a*IkBd-p.k19*IkBdn-p.k51*IkBdn*RelAp50n+p.k54*IkBdRelAp50n...
-p.k53*IkBdn*RelBp50n-p.k53a*IkBdn*RelBp52n+p.k57*IkBdRelBp50n+p.k57*IkBdRelBp52n...
+p.k49*IkBdRelAp50n+p.k49*IkBdRelBp50n+p.k49*IkBdRelBp52n-p.k25*IkBdn;%IkBdn 
dydt(16,1)= tRelA*p.k10-p.k40*RelA-p.k41*RelA*p50+p.k44*RelAp50;%RelA
dydt(17,1)= p.k44*RelAp50n-p.k40*RelAn-p.k41*RelAn*p50n;%RelAn
dydt(18,1)=  tRelB*p.k10-p.k40*RelB-p.k41*RelB*p52-p.k42*RelB*p50-p.k43*RelB*p100+p.k44*RelBp52...
+p.k45*RelBp50+p.k46*RelBp100;%RelB
dydt(19,1)=  p.k44*RelBp52n+p.k45*RelBp50n+p.k46*RelBp100n-p.k40*RelBn-p.k41*RelBn*p52n...
-p.k42*RelBn*p50n-p.k43*RelBn*p100n;%RelBn
dydt(20,1)= tp50*p.k10-p.k41*RelA*p50-p.k42*RelB*p50+p.k44*RelAp50+p.k45*RelBp50-p.k40*p50;%p50 
dydt(21,1)= p.k44*RelAp50n+p.k45*RelBp50n-p.k41*RelAn*p50n-p.k42*RelBn*p50n-p.k40*p50n;%p50n
dydt(22,1)=  tp100*p.k10-p.k40*p100-NIK*p.k32*p100/(p.k33+p100)-p.k15*p100*p100+p.k16*IkBd-p.k43*RelB*p100...
+p.k46*RelBp100;%p100
dydt(23,1)= p.k46*RelBp100n-p.k40*p100n-p.k15*p100n*p100n+p.k16*IkBdn-p.k43*RelBn*p100n;%p100n 
dydt(24,1)= NIK*p.k32*p100/(p.k33+p100)+p.k44*RelBp52-p.k41*RelB*p52-p.k40*p52;%p52 
dydt(25,1)= p.k44*RelBp52n-p.k41*RelBn*p52n-p.k40*p52n;%p52n
dydt(26,1)=  p.k41*RelA*p50-p.k44*RelAp50-p.k47*RelAp50+p.k48*RelAp50n-p.k49*RelAp50-p.k51*IkBa*RelAp50...
-p.k51*IkBb*RelAp50-p.k51*IkBe*RelAp50-p.k51*IkBd*RelAp50+p.k54*IkBaRelAp50+p.k55*IkBbRelAp50...
+p.k54*IkBeRelAp50+p.k54*IkBdRelAp50+NIK*p.k30*IkBdRelAp50/(p.k31+IkBdRelAp50)...
+p.k27*POLYIC*IkBaRelAp50+p.k28*POLYIC*IkBbRelAp50+p.k29*POLYIC*IkBeRelAp50;%RelAp50
dydt(27,1)=  p.k41*RelAn*p50n-p.k44*RelAp50n+p.k47*RelAp50-p.k48*RelAp50n-p.k49*RelAp50n-p.k51*IkBan*RelAp50n...
-p.k51*IkBbn*RelAp50n-p.k51*IkBen*RelAp50n-p.k51*IkBdn*RelAp50n+p.k54*IkBaRelAp50n+p.k55*IkBbRelAp50n...
+p.k54*IkBeRelAp50n+p.k54*IkBdRelAp50n;%RelAp50n
dydt(28,1)=  p.k41*RelB*p52-p.k44*RelBp52-p.k47*RelBp52+p.k48*RelBp52n-p.k50*RelBp52+p.k57*IkBdRelBp52...
+NIK*p.k30*IkBdRelBp52/(p.k31+IkBdRelBp52)-p.k53a*IkBd*RelBp52;%RelBp52
dydt(29,1)=  p.k41*RelBn*p52n-p.k44*RelBp52n+p.k47*RelBp52-p.k48*RelBp52n-p.k50*RelBp52n+p.k57*IkBdRelBp52n...
-p.k53a*IkBdn*RelBp52n;%RelBp52n
dydt(30,1)= p.k43*RelB*p100-p.k46*RelBp100-p.k49*RelBp100-p.k48*RelBp100+p.k48*RelBp100n;%RelBp100 
dydt(31,1)= p.k43*RelBn*p100n-p.k46*RelBp100n-p.k49*RelBp100n+p.k48*RelBp100-p.k48*RelBp100n;%RelBp100n 
dydt(32,1)= p.k51*IkBa*RelAp50-p.k54*IkBaRelAp50+p.k60*IkBaRelAp50n-p.k58*IkBaRelAp50-p.k49*IkBaRelAp50...
-p.k27*POLYIC*IkBaRelAp50;%IkBaRelAp50
dydt(33,1)=  p.k51*IkBan*RelAp50n-p.k54*IkBaRelAp50n-p.k60*IkBaRelAp50n+p.k58*IkBaRelAp50...
-p.k49*IkBaRelAp50n;%IkBaRelAp50n
dydt(34,1)=  p.k51*IkBb*RelAp50-p.k55*IkBbRelAp50+p.k61*IkBbRelAp50n-p.k58a*IkBbRelAp50-p.k49*IkBbRelAp50...
-p.k28*POLYIC*IkBbRelAp50;%IkBbRelAp50
dydt(35,1)=  p.k51*IkBbn*RelAp50n-p.k55*IkBbRelAp50n-p.k61*IkBbRelAp50n+p.k58a*IkBbRelAp50...
-p.k49*IkBbRelAp50n;%IkBbRelAp50n
dydt(36,1)=  p.k51*IkBe*RelAp50-p.k54*IkBeRelAp50+p.k61*IkBeRelAp50n-p.k59*IkBeRelAp50-p.k49*IkBeRelAp50...
-p.k29*POLYIC*IkBeRelAp50;%IkBeRelAp50
dydt(37,1)=  p.k51*IkBen*RelAp50n-p.k54*IkBeRelAp50n-p.k61*IkBeRelAp50n+p.k59*IkBeRelAp50...
-p.k49*IkBeRelAp50n;%IkBeRelAp50n
dydt(38,1)=  p.k51*IkBd*RelAp50-p.k54*IkBdRelAp50+p.k61*IkBdRelAp50n-p.k58a*IkBdRelAp50-p.k49*IkBdRelAp50...
-NIK*p.k30*IkBdRelAp50/(p.k31+IkBdRelAp50);%IkBdRelAp50
dydt(39,1)=  p.k51*IkBdn*RelAp50n-p.k54*IkBdRelAp50n-p.k61*IkBdRelAp50n+p.k58a*IkBdRelAp50...
-p.k49*IkBdRelAp50n;%IkBdRelAp50n
dydt(40,1)=  p.k42*RelB*p50-p.k45*RelBp50-p.k47*RelBp50+p.k48*RelBp50n-p.k49*RelBp50-p.k52*IkBa*RelBp50...
-p.k52*IkBe*RelBp50-p.k53*IkBd*RelBp50+p.k56*IkBaRelBp50+p.k56*IkBeRelBp50+p.k57*IkBdRelBp50...
+NIK*p.k30*IkBdRelBp50/(p.k31+IkBdRelBp50)+p.k27*POLYIC*IkBaRelBp50+p.k29*POLYIC*IkBeRelBp50;%RelBp50 
dydt(41,1)= p.k42*RelBn*p50n-p.k45*RelBp50n+p.k47*RelBp50-p.k48*RelBp50n-p.k49*RelBp50n-p.k52*IkBan*RelBp50n...
-p.k52*IkBen*RelBp50n-p.k53*IkBdn*RelBp50n+p.k56*IkBaRelBp50n+p.k56*IkBeRelBp50n...
+p.k57*IkBdRelBp50n;%RelBp50n
dydt(42,1)=  p.k52*IkBa*RelBp50-p.k56*IkBaRelBp50+p.k60*IkBaRelBp50n-p.k58*IkBaRelBp50-p.k49*IkBaRelBp50...
-p.k27*POLYIC*IkBaRelBp50;%IkBaRelBp50
dydt(43,1)=  p.k52*IkBan*RelBp50n-p.k56*IkBaRelBp50n-p.k60*IkBaRelBp50n+p.k58*IkBaRelBp50...
-p.k49*IkBaRelBp50n;%IkBaRelBp50n
dydt(44,1)=  p.k52*IkBe*RelBp50-p.k56*IkBeRelBp50+p.k61*IkBeRelBp50n-p.k59*IkBeRelBp50-p.k49*IkBeRelBp50...
-p.k29*POLYIC*IkBeRelBp50;%IkBeRelBp50
dydt(45,1)=  p.k52*IkBen*RelBp50n-p.k56*IkBeRelBp50n-p.k61*IkBeRelBp50n+p.k59*IkBeRelBp50...
-p.k49*IkBeRelBp50n;%IkBeRelBp50n
dydt(46,1)=  p.k53*IkBd*RelBp50-p.k57*IkBdRelBp50+p.k61*IkBdRelBp50n-p.k58a*IkBdRelBp50-p.k49*IkBdRelBp50...
-NIK*p.k30*IkBdRelBp50/(p.k31+IkBdRelBp50);%IkBdRelBp50
dydt(47,1)=  p.k53*IkBdn*RelBp50n-p.k57*IkBdRelBp50n-p.k61*IkBdRelBp50n+p.k58a*IkBdRelBp50...
-p.k49*IkBdRelBp50n;%IkBdRelBp50n
dydt(48,1)=  p.k53a*IkBd*RelBp52-p.k57*IkBdRelBp52+p.k61*IkBdRelBp52n-p.k58a*IkBdRelBp52...
-p.k49*IkBdRelBp52-NIK*p.k30*IkBdRelBp52/(p.k31+IkBdRelBp52);%IkBdRelBp52 
dydt(49,1)= p.k53a*IkBdn*RelBp52n-p.k57*IkBdRelBp52n-p.k61*IkBdRelBp52n+p.k58a*IkBdRelBp52...
-p.k49*IkBdRelBp52n;%IkBdRelBp52n
end
