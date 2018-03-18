function [WT_FWD,WT_REV,NoBias_FWD,NoBias_REV,WithBias_FWD,WithBias_REV] = PromoterPredictionInEcoliGenome_MAIN(start_gene_index)

% PARAMS 
permutations_num = 1000;
genes_num_save_pause = 400;
temp_work_path = '/tmp/AY/';
results_output_file = ['/Users/avihuy/Applications/bprom_mac/Ecoli_Results/' 'all_gene_' num2str(permutations_num) '_permutaions'];

mkdir(temp_work_path);

% load e.coli genes data struct
load EColi_Codon_Bias;
load EColi_No_Bias;
load EcoliGeneBankData;

cd '/Users/avihuy/Applications/bprom_mac/';
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');

all_gene_names = {EcoliGeneBankData.CDS.gene};
num_of_genes = length(all_gene_names);

% loading previous results from unfinished run
if (start_gene_index) > 1
    
    eval(['load ' results_output_file ';']);
    
end


for curr_gene_index=start_gene_index:num_of_genes
    
    curr_gene_name = all_gene_names{curr_gene_index};
    
    % Getting current gene sequence
    curr_gene_start_pos = EcoliGeneBankData.CDS(curr_gene_index).indices(1);
    curr_gene_end_pos = EcoliGeneBankData.CDS(curr_gene_index).indices(end);

    if curr_gene_start_pos<curr_gene_end_pos
        
        curr_gene_seq = EcoliGeneBankData.Sequence(curr_gene_start_pos:curr_gene_end_pos);
        
    else
        
        curr_gene_seq = EcoliGeneBankData.Sequence(curr_gene_end_pos:curr_gene_start_pos);
        curr_gene_seq = seqrcomplement(curr_gene_seq);
        
    end
    
    % The number of predicted promoters in the WT gene - here FWD means 'sense'
    curr_gene_num_of_WT_predicted_promoters_FWD = BPROMPromoterPredictionForASingleSeq(temp_work_path,curr_gene_seq);
    WT_FWD(curr_gene_index) = curr_gene_num_of_WT_predicted_promoters_FWD;
    
    % The number of predicted promoters in the WT REVERSE COMPLEMENT ('antisense')
    curr_gene_num_of_WT_predicted_promoters_REV = BPROMPromoterPredictionForASingleSeq(temp_work_path,seqrcomplement(curr_gene_seq));
    WT_REV(curr_gene_index) = curr_gene_num_of_WT_predicted_promoters_REV;

    % Adjusting the aa seq for the permutations of BOTH no_bias and Ecoli_bias below
    % changing the stop codon simbole from * to X for easier handeling
    curr_gene_aa_seq = upper(nt2aa(curr_gene_seq));
    curr_gene_aa_seq(curr_gene_aa_seq=='*')='X';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a 2 folders with curr gene coding permutations
    % according to a the NO BIAS STRRUCT of Ecoli 
    % both FWD and REV
    
    CreateDirectoriesWithProteinPermutationsBothFWDandREV(temp_work_path,curr_gene_name,curr_gene_aa_seq,EColi_No_Bias,permutations_num);
      
    % Run BPROM prediction program on the two folders  
    curr_gene_FWD_folder = [temp_work_path curr_gene_name '_' num2str(permutations_num) '_permutations' '_FWD'];
    curr_gene_REV_folder = [temp_work_path curr_gene_name '_' num2str(permutations_num) '_permutations' '_REV'];

    unix(['/Users/avihuy/Applications/bprom_mac/LoopBPROMonFolder.sh ' curr_gene_FWD_folder]);
    unix(['/Users/avihuy/Applications/bprom_mac/LoopBPROMonFolder.sh ' curr_gene_REV_folder]);
    
    % Collect results(right now it is just the num of predicted promoters) from the two folders 
    curr_gene_FWD_predicted_promoters_array = GetBPromResultsFromFolder(curr_gene_FWD_folder,'*.txt');
    curr_gene_REV_predicted_promoters_array = GetBPromResultsFromFolder(curr_gene_REV_folder,'*.txt');

    % Put curr gene data into all genes data structure
    NoBias_FWD(curr_gene_index,1:permutations_num) = curr_gene_FWD_predicted_promoters_array;
    NoBias_REV(curr_gene_index,1:permutations_num) = curr_gene_REV_predicted_promoters_array;
    
    % Remove work folders (FWD and REV) for current gene
    eval(['rmdir(''' curr_gene_FWD_folder '''' ',' '''s''' ');']);
    eval(['rmdir(''' curr_gene_REV_folder '''' ',' '''s''' ');']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a 2 folders with curr gene coding permutations
    % according to a the EColi BIAS STRRUCT of Ecoli 
    % both FWD and REV
    
    CreateDirectoriesWithProteinPermutationsBothFWDandREV(temp_work_path,curr_gene_name,curr_gene_aa_seq,EColi_Codon_Bias,permutations_num);
      
    % Run BPROM prediction program on the two folders  
    curr_gene_FWD_folder = [temp_work_path curr_gene_name '_' num2str(permutations_num) '_permutations' '_FWD'];
    curr_gene_REV_folder = [temp_work_path curr_gene_name '_' num2str(permutations_num) '_permutations' '_REV'];

    unix(['/Users/avihuy/Applications/bprom_mac/LoopBPROMonFolder.sh ' curr_gene_FWD_folder]);
    unix(['/Users/avihuy/Applications/bprom_mac/LoopBPROMonFolder.sh ' curr_gene_REV_folder]);
    
    % Collect results(right now it is just the num of predicted promoters) from the two folders 
    curr_gene_FWD_predicted_promoters_array = GetBPromResultsFromFolder(curr_gene_FWD_folder,'*.txt');
    curr_gene_REV_predicted_promoters_array = GetBPromResultsFromFolder(curr_gene_REV_folder,'*.txt');

    % Put curr gene data into all genes data structure
    WithBias_FWD(curr_gene_index,1:permutations_num) = curr_gene_FWD_predicted_promoters_array;
    WithBias_REV(curr_gene_index,1:permutations_num) = curr_gene_REV_predicted_promoters_array;
    
    % Remove work folders (FWD and REV) for current gene
    eval(['rmdir(''' curr_gene_FWD_folder '''' ',' '''s''' ');']);
    eval(['rmdir(''' curr_gene_REV_folder '''' ',' '''s''' ');']);
    
    disp(['Completed gene ' curr_gene_name ' - ' num2str(curr_gene_index)  ' out of ' num2str(num_of_genes)]);
    
    % Partial save of results and allow user to pause run until later resume
    if mod(curr_gene_index,genes_num_save_pause)==0
        
        % Save results
        save(results_output_file,'WT_FWD','WT_REV','NoBias_FWD','NoBias_REV','WithBias_FWD','WithBias_REV');
        disp('results are saved!');
        
    end
    
end

% Save FINAL results
save(results_output_file,'WT_FWD','WT_REV','NoBias_FWD','NoBias_REV','WithBias_FWD','WithBias_REV');
disp('Final results are saved!');

