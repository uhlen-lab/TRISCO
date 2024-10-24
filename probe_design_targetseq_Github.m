function [summary_table] = probe_design_targetseq_Github(FASTA_dir, out_dir)

% ChatGPT checked my original code and modified to compatible with MATLAB2023a.

% FASTA_dir should has folders for each gene. The folders for gene should have one text FASTA file.
% FASTA_dir = 'E:\Shigeaki\Probe_design\Probes';

% Probelist is a table of [Gene, Length, hairpin, Start1, End1, Start2, End2 .... StartN, EndN];
% [Gene_name, length_of_gene, hairpin(B1-B5), the area you want to remove showed by pairs of (start-end)]

% set max_probe_number = -1 if you don't want to reduce the number of probes

% Error log
% 1. FASTA_dir should be the folder which contain gene folder such as "Mouse_Ang2" etc. directly
% 2. Mismatch of gene list between Probelist and folders in FASTA_dir
% 3. Hairpin type is strange (B14 etc.)

%% Make gene list from folder structure

file_list = struct2table(dir(FASTA_dir));

folder_list = file_list(file_list.isdir, :);
folder_list = folder_list(3:end, :);
num_of_gene = length(folder_list{:,1});

%% preallocation

% All_probe_list= {}; 
% hairpin_type = strings(1, num_of_gene);
% Fasta_header = strings(1, num_of_gene);
% Fasta_seq = strings(1, num_of_gene);
% error_flag = 0;


targets_name_all = table2array(folder_list(:,1));
hairpin_type_all = cell(num_of_gene, 1);
target_seq_all   = cell(num_of_gene, 1);
length_all       = cell(num_of_gene, 1);

% All_probe_list= {}; 
% hairpin_type = strings(num_of_gene, 1);
% Fasta_header = strings(num_of_gene, 1);
% Fasta_seq = strings(1, num_of_gene);
% error_flag = 0;


for i=1:num_of_gene
%% Reading folder structure and files
    % Folder specifying 
    target_name = folder_list.name{i};    
    gene_dir = fullfile(FASTA_dir, char(target_name)); 

    % Exclude file reading
    gene_excludetxt = dir(fullfile(gene_dir, '*.csv'));
    fname_exclude = fullfile(gene_excludetxt.folder, gene_excludetxt.name); 
    exclude_table = readtable(fname_exclude);
    exclude_array = table2array(exclude_table);

    if isempty(exclude_array) == 1
        rem_start = [];
        rem_end   = [];
    else
        rem_start = exclude_array(:,1);
        rem_end   = exclude_array(:,2);
    end

    % FASTA file reading
    gene_fastatxt = dir(fullfile(gene_dir, '*.fasta'));
    fname_fasta = fullfile(gene_fastatxt.folder, gene_fastatxt.name); 
    Fasta_struct = fastaread(fname_fasta);
    Fasta_cell = struct2cell(Fasta_struct);
    %Fasta_length = length(Fasta_cell{2,1});
    Fasta_length = length(Fasta_struct.Sequence);

    % Hairpin type reading
    gene_hairpintxt = dir(fullfile(gene_dir, '*.txt'));
    fname_hairpin = fullfile(gene_fastatxt.folder, gene_hairpintxt.name);
    hairpin_type = fileread(fname_hairpin);

 %% Get the targetting sequence region 
    position = true(Fasta_length,1); % position is logical (now all true and change the seq to false later)

    for p=1:length(rem_start) % the number of the pair of 'rem_start' and 'rem_end'
    position(rem_start(p):rem_end(p)) = 0; % change to false if sequence is removed
    
    end

    % The script to retrieve the position of start and end
    a = find(position==1); % get sequence area which is true

    D = diff([0;diff(a)==1;0]); % get subtraction from next position and convert =1 (start position), =-1 (end position)
    first = a(D>0);  % to get start position
    last = a(D<0);   % to get end position
    target_seq = [first last];
    target_length = 1 + find(D<0) - find(D>0); % Calculate sequence length

    output_list = [target_seq target_length];
    target_seq_column = reshape(target_seq.',1,[]);  %  array of Start1, End1, Start2, End2, ..

 %% Collect the calculated information and store it in the table

    hairpin_type_all{i} = hairpin_type;
    target_seq_all{i}   = target_seq;
    length_all{i}       = Fasta_length;

end 

%% Create the Matlab table and save
summary_table = table(targets_name_all, length_all, hairpin_type_all, target_seq_all); % continue adding all columns
% Label the table headers
summary_table.Properties.VariableNames = {'Gene', 'Length', 'Hairpin', 'Target sequence region'};
% Save the Matlab table
save([out_dir '/TargetSequence_Summary.mat'], 'summary_table');

%% Create the excel table
target_seq_all_vector = cellfun(@(x) x(:), cellfun(@(x) x.', target_seq_all, 'UniformOutput', false), 'UniformOutput', false);
summary_table_excel = table(targets_name_all, length_all, hairpin_type_all, target_seq_all_vector); % continue adding all columns

% Label the table headers
summary_table_excel.Properties.VariableNames = {'Gene', 'Length', 'Hairpin', 'Target sequence region'};
 
% Save the table
writetable(summary_table_excel, [out_dir '/TargetSequence_Summary.csv']);

% % Set the base headers
% baseHeaders = {'Gene', 'Length', 'hairpin'};
% 
% % Determine the number of Start/End pairs dynamically
% % For a 1xN cell array, where each cell might have a different number of columns
% rowLengths = cellfun(@(x) size(x, 1), target_seq_all);
% max_size = max(rowLengths);
% 
% % Initialize the headers with base headers
% headers = baseHeaders;
% 
% % Dynamically add 'Start' and 'End' headers with numbering
% for i = 1:max_size
%     headers{end+1} = sprintf('Start%d', i);
%     headers{end+1} = sprintf('End%d', i);
% end
% % Now headers will be:
% % {'Gene', 'Length', 'hairpin', 'Start1', 'End1', 'Start2', 'End2', 'Start3', 'End3', 'Start4', 'End4', 'Start5', 'End5'}
% 
% summary_table.Properties.VariableNames = headers;

target_seq_all_vector = target_seq_all{1}.';
target_seq_all_vector = target_seq_all_vector(:);

target_seq_all_vector2 = cellfun(@(x) x(:), cellfun(@(x) x.', target_seq_all, 'UniformOutput', false), 'UniformOutput', false);




