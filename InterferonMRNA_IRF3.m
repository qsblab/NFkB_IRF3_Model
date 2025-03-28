% Estimation of Interferon mRNA Level Module: Promoter Binding Activity 
% when only IRF-3 is bound upstream of Type I IFN genes
% =========================================================================

% Filename
filename = 'extracted_values.txt';

% Open File
fid = fopen(filename,'r');

% Reading Header
header = fgetl(fid);

% Calculating the number of Line (numLines) in the external file
% =========================================================================
% Counter Initialization
numLines = 0;

% Counting the number of lines in External File
while ~feof(fid)
    fgetl(fid);
    numLines = numLines+1;
end

% Rewind the file to the beginning
frewind(fid);

% Skip Header Line
fgetl(fid);

% Cell Array Initialization
IRF3_values = cell(numLines,1);
NFkB_values = cell(numLines,1);
polyIC_indices = zeros(numLines,1);
time_indices = zeros(numLines,1);

% Reading External File to Store Data in Desired Format
lineCounter = 1;
while ~feof(fid)
    % Reading Data Lines
    line = fgetl(fid);
    
    % Split Lines into Component 
    parts = strsplit(line,'\t'); % For tab-seperated values

    %Extract PolyIC and Time Indices
    polyIC_indices(lineCounter) = str2double(parts{1});
    time_indices(lineCounter) = str2double(parts{2});

    % Extract IRF3 (3rd Column in 'extracted_files.txt') and convert from
    % string to numerical array
    IRF3_values{lineCounter} = str2num(parts{3});
    
    % Extract NFkB (4th Column in 'extracted_files.txt') and convert from
    % string to numerical array
    if length(parts) >= 4
        NFkB_values{lineCounter} = str2num(parts{4});
    else
        NFkB_values{lineCounter} = [];
    end

    % Increment line counter
    lineCounter = lineCounter + 1;
end

% Close the file
fclose(fid);

% Calculating IFN mRNA Level through Promoter Binding Activity using
% extracted values
results = cell(numLines,1);

for idx = 1:numLines
    
    if ~isempty (NFkB_values{idx})
        
        % To check that both arrays are of the same length
        min_length = min(length(IRF3_values{idx}),length(NFkB_values{idx}));
        irf3 = cell(numLines,1);
        for a = min_length
            irf3 = IRF3_values{idx,1}(1:a);
        end
        nfkb = cell (numLines,1);
        for b = 1: min_length
            nfkb = NFkB_values{idx,1}(1:b);
        end
        k_nfkb = 20;
        k_irf = 10;
        n = 1;

        % Calculating IFN mRNA Level and Storing in "results" array

        results{idx,1} = (((irf3.^n)./((k_irf.^n)+(irf3.^n))).*(1-((nfkb.^n)./((k_nfkb.^n)+(nfkb.^n))))).*100;
    else
        results{idx}=[];
    end

end

% Result Visualization
% =========================================================================

time = 1:301;   
figure;
plot_index = 1;
for i = 1:length(polyIC_indices) % PolyI:C indices
    for j = 1:2 %Time indices
        result_idx = find(polyIC_indices == i & time_indices == j);
        if ~isempty(result_idx)
            ifn_value=results{result_idx,1};                         % Command to extract cell array's row value: row{x,1}
            subplot(1,2,plot_index);                                 % Command to create subplots in the same figure present subplot represent current index of cell array accessed 
            plot(time,ifn_value,'LineWidth',3);
            title(['PolyI:C Index' num2str(i) ', Time Value ' num2str(j)]);
            xlabel('Time (in min)');
            ylabel('Interferon mRNA Level (in nM)');
            xlim([0 300]);
            ylim([0 100]);
            plot_index = plot_index + 1;
        end
    end
end
disp('Results succesfully plotted');



