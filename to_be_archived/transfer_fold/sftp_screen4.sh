#!/bin/bash

# List of files
files=("scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_16.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_17.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_18.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_19.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_20.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_21.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_22.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_23.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_24.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_25.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_26.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_27.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_28.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_29.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_30.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_31.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_32.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_33.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_34.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_35.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_36.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_37.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_38.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_39.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_40.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_41.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_42.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_43.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_44.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_45.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_46.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_47.RDS"                                          
"scratch/oscDetectPaper/transfer_fold/cvxr/solns_cvxr_48.RDS")

# Loop through each file
for file in "${files[@]}"; do
    # Create a new screen session with a unique name based on the file
    screen -d -m  bash -c "sftp beluga << EOF
    get $file
    quit
    EOF
    "
done
