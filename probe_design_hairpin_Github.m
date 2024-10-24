function probe_design_hairpin_Github(FASTA_dir, summary_table, out_dir, max_probe_number)

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

targets_name_all = summary_table{:,1};
length_all = summary_table{:,2};
hairpin_type_all = summary_table{:,3};
target_seq_all = summary_table{:,4};

num_of_gene = length(summary_table{:,1});

% file_list = struct2table(dir(FASTA_dir));
% folder_list = file_list(file_list.isdir, :);
% folder_list = folder_list(3:end, :);

%% preallocation

 All_probe_list= {}; 
 Fasta_header = strings(1, num_of_gene);
 Fasta_seq = strings(1, num_of_gene);
% error_flag = 0;
%% Probe sequence design (Tiling)

for i=1:num_of_gene

    % Get 
    target_name = targets_name_all{i};
    hairpin_type = hairpin_type_all{i};
    target_seq = target_seq_all{i}; % The number of removal sequence area

    % Get FASTA sequence infor
    gene_dir = fullfile(FASTA_dir, char(target_name)); 
    gene_fastatxt = dir(fullfile(gene_dir, '*.fasta'));
    fname = fullfile(gene_fastatxt.folder, gene_fastatxt.name); 
    Fasta_struct = fastaread(fname);
    Fasta_cell = struct2cell(Fasta_struct);

    Fasta_header(i) = Fasta_cell{1, 1};
    Fasta_seq(i) = Fasta_cell{2, 1};
    
    Probe_list = [];
            
       for s=1:size(target_seq,1)
    
           start_seq = target_seq(s, 1); 
           end_seq   = target_seq(s, 2); 
                     
           Fasta_seq_part = char(Fasta_seq(i));
           Fasta_seq_use = Fasta_seq_part(start_seq:end_seq);
            
           [Export_list] = probe_design_version3_HCR(Fasta_seq_use, hairpin_type, start_seq);

            Probe_list = cat(1, Probe_list, Export_list);
              
        end
    
     All_probe_list{i} = Probe_list;
     
end

%% save data
    
% preallocation
    probe_pair_number = zeros(length(targets_name_all),1); % The number of hairpin pair
        
for i=1:num_of_gene
    
    Odd_header = strcat('Odd probe_', hairpin_type_all(i));
    Even_hearder = strcat('Even probe_', hairpin_type_all(i));
    Export_list_Header = cat(2, Odd_header, 'Start', 'End', Even_hearder, 'Start', 'end');
    
    Fasta_header_cell = {Fasta_header(i), '', '', '', '' , ''};  % to make compatible to cat
    Fasta_header_string = string(Fasta_header_cell);
    
    probe_pair_number(i) = size(All_probe_list{i}, 1);
    
    Export_list = cat(1, Fasta_header_string, Export_list_Header, All_probe_list{i});
      
    writematrix(Export_list, fullfile(out_dir, ['Probe_list_'  targets_name_all{i} '.csv']));
    
    save(fullfile(out_dir, 'All_probe_list.mat'), 'All_probe_list');  
        
end

%% Save summary

     Header_table = ["Target_name", "Sequence_length", "Probe_pair_number", "Probe_Hapirpin"];
     ProbeSequence_length_list = probe_pair_number .* 90 ; % one pair length is 45*2 = 90
     Summary_table_pre = cat(2, string(targets_name_all), ProbeSequence_length_list, probe_pair_number, hairpin_type_all);
     
     Summary_table = cat(1, Header_table, Summary_table_pre);
     writematrix(Summary_table, fullfile(out_dir, 'Probe_design_summary.csv'));

%% Make excel sheet for opool oligo
     
    opool_header = ["Pool name", "Sequence"];
    opool_name_list = strcat(targets_name_all, '_', hairpin_type_all);  % Make order names for opool
    
    opool_excel_list = {};

    for i=1:length(opool_name_list)
        
        Gene_probe_list = All_probe_list{i};
        
        if isempty(Gene_probe_list) == 1    % If there is no target sequence available for some gene, skip the iteration. 
            continue;
        else

        probe_list_length = size(Gene_probe_list, 1);
        text_store = strings(probe_list_length .* 2, 1);
        text_store(:) = opool_name_list(i);
        
        opool_sequence_list = cat(1, Gene_probe_list(:,1), Gene_probe_list(:,4));
        opool_export_list = cat(2, text_store, opool_sequence_list);

        opool_excel_list = [opool_excel_list; opool_export_list];
        
        end

    end

     opool_excel_list_save = cat(1, opool_header, opool_excel_list);
     
     save(fullfile(out_dir, 'opool_excel_list_save.mat'), 'opool_excel_list_save');  
     
    opool_table = array2table(opool_excel_list, 'VariableNames', opool_header);
    writetable(opool_table, fullfile(out_dir, 'opool_probe_sheet.xlsx'));
 
 %% reduce probe size if necessary
    
if max_probe_number > 0
    
        All_probe_list_reduced = probe_set_reduce_ver1(All_probe_list, max_probe_number);
        mkdir(fullfile(out_dir, 'Reduced'));

    %% save data
    for i=1:num_of_gene
        probe_pair_number(i) = size(All_probe_list_reduced{i}, 1);        

        Fasta_header_cell = {Fasta_header(i), '', '', '', '' , ''};  % to make compatible to cat
        Fasta_header_string = string(Fasta_header_cell);
        
        Odd_header = strcat('Odd probe_', hairpin_type_all(i));
        Even_hearder = strcat('Even probe_', hairpin_type_all(i));
        Export_list_Header = cat(2, Odd_header, 'Start', 'End', Even_hearder, 'Start', 'end');

        Export_list = [Fasta_header_string; Export_list_Header; All_probe_list_reduced{i}];

        writematrix(Export_list, fullfile(out_dir, 'Reduced', ['Probe_list_' targets_name_all{i} '.csv']));

        save(fullfile(out_dir, 'Reduced', 'All_probe_list_reduced.mat'), 'All_probe_list_reduced');
    end

%% Save summary
    Header_table = ["Target_name", "Sequence_length", "Probe_pair_number", "Probe_Hapirpin"];
    ProbeSequence_length_list = probe_pair_number .* 90 ; % one pair length is 45*2 = 90
    Summary_table_pre = [string(targets_name_all), ProbeSequence_length_list, probe_pair_number, hairpin_type_all];

    Summary_table = [Header_table; Summary_table_pre];
    writematrix(Summary_table, fullfile(out_dir, 'Reduced', 'Probe_design_summary_reduced.csv'));

%% Make excel sheet for opool oligo
    opool_header = ["Pool name", "Sequence"];
    opool_name_list = strcat(targets_name_all, '_', hairpin_type_all);  % Make order names for opool

    opool_excel_list = {};

    for i=1:length(opool_name_list)
        Gene_probe_list = All_probe_list_reduced{i};
        
        if isempty(Gene_probe_list) == 1    % If there is no target sequence available for some gene, skip the iteration. 
            continue;
        else

        probe_list_length = size(Gene_probe_list, 1);
        text_store = strings(probe_list_length .* 2, 1);
        text_store(:) = opool_name_list(i);

        opool_sequence_list = [Gene_probe_list(:,1); Gene_probe_list(:,4)];
        opool_export_list = [text_store, opool_sequence_list];

        opool_excel_list = [opool_excel_list; opool_export_list];

        end 
    end

    opool_excel_list_save = [opool_header; opool_excel_list];

    save(fullfile(out_dir, 'Reduced', 'opool_excel_list_save_reduced.mat'), 'opool_excel_list_save');

    opool_table = array2table(opool_excel_list, 'VariableNames', opool_header);
    writetable(opool_table, fullfile(out_dir, 'Reduced', 'opool_probe_sheet_reduced.xlsx'));

end

    
end


