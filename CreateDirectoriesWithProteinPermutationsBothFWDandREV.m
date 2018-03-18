function CreateDirectoriesWithProteinPermutationsBothFWDandREV(path,gene_name,aa_seq,codon_bias_struct,permutations_num)

output_folder_name_fwd = [path '/' gene_name '_' num2str(permutations_num) '_permutations_FWD'];
output_folder_name_rev = [path '/' gene_name '_' num2str(permutations_num) '_permutations_REV'];
mkdir(output_folder_name_fwd);
mkdir(output_folder_name_rev);

aa_seq = upper(aa_seq);


protein_permutations = CodonBiasPermutationsOfProtein(aa_seq,codon_bias_struct,permutations_num);


for i=1:size(protein_permutations,2)

    curr_filename = [output_folder_name_fwd '/' gene_name '_' num2str(i) '_FWD' '.fa'];
    curr_mutated_seq = protein_permutations{i};
    fileID = fopen(curr_filename,'w');
    fprintf(fileID,'%s\n',['>' gene_name '_' num2str(i) '_FWD']);

    for k=1:floor(length(curr_mutated_seq)/80)
        fprintf(fileID,'%s\n',curr_mutated_seq((k-1)*80+1:k*80));
    end
    if isempty(k) % this is in case seq len is less than 80
        k=0; 
    end  
    fprintf(fileID,'%s',curr_mutated_seq(k*80+1:end));    
    fclose(fileID);
    
    % Creating a seperate file, in a seperate folder for the REV seq
    
    curr_filename = [output_folder_name_rev '/' gene_name '_' num2str(i) '_REV' '.fa'];
    curr_mutated_seq = seqrcomplement(protein_permutations{i});
    fileID = fopen(curr_filename,'w');
    fprintf(fileID,'%s\n',['>' gene_name '_' num2str(i) '_REV']);

    for k=1:floor(length(curr_mutated_seq)/80)
        fprintf(fileID,'%s\n',curr_mutated_seq((k-1)*80+1:k*80));
    end
    if isempty(k) % this is in case seq len is less than 80
        k=0; 
    end  
    fprintf(fileID,'%s',curr_mutated_seq(k*80+1:end));    
    fclose(fileID);
    
    
end

