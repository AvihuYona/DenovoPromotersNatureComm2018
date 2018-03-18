function predicted_promoters_array = GetBPromResultsFromFolder(path,file_mask)

% Determine how many files we have in the folder
files=dir([path '/' file_mask ]);

for i=1:length(files)
    
    curr_filename = files(i).name;

    [curr_predicted_promoters_num,curr_total_BS_strength,curr_total_LDF_sum, curr_total_expression_score ]= GetPredictedExpressionCapacityFromBPromOutputFile(path,curr_filename);
    predicted_promoters_array(i).num = curr_predicted_promoters_num;
    predicted_promoters_array(i).BS_sum = curr_total_BS_strength;
    predicted_promoters_array(i).LDF_sum = curr_total_LDF_sum;
    predicted_promoters_array(i).total_expression_score = curr_total_expression_score;


    
end

