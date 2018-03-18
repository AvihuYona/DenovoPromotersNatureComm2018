function predicted_promoters = BPROMPromoterPredictionForASingleSeq(tmp_working_folder,seq)

% Creating the input file for BPROM 
folder_name = [tmp_working_folder '/' 'BPROM_temp_dir'];
mkdir(folder_name);
tmp_filename = [folder_name '/'  'tmp_seq_file.fa'];
fileID = fopen(tmp_filename,'w');
fprintf(fileID,'%s\n',['>' 'stam test seq']);

for k=1:floor(length(seq)/80)
    fprintf(fileID,'%s\n',seq((k-1)*80+1:k*80));
end
if isempty(k) % this is in case seq len is less than 80
    k=0; 
end  
fprintf(fileID,'%s',seq(k*80+1:end));    
fclose(fileID);

% Running BPROM on the file created
cd '/Users/avihuy/Applications/bprom_mac/';
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/');
unix(['/Users/avihuy/Applications/bprom_mac/LoopBPROMonFolder.sh ' folder_name]);

% Reading BPROM results
predicted_promoters = GetBPromResultsFromFolder(folder_name,'*.txt');
predicted_promoters = predicted_promoters(1);

% Removing temp dir and files
rmdir(folder_name, 's');
