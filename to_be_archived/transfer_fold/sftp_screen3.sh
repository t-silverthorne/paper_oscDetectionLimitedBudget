#!/bin/bash

# List of files
files=("scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_16.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_17.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_18.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_19.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_20.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_21.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_22.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_23.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_24.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_25.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_26.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_27.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_28.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_29.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_30.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_31.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_32.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_33.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_34.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_35.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_36.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_37.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_38.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_39.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_40.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_41.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_42.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_43.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_44.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_45.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_46.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_47.RDS"                                    
"scratch/oscDetectPaper/transfer_fold/pin2lat/solns_pin2lat_48.RDS")

# Loop through each file
for file in "${files[@]}"; do
    # Create a new screen session with a unique name based on the file
    screen -d -m  bash -c "sftp beluga << EOF
    get $file
    quit
    EOF
    "
done
