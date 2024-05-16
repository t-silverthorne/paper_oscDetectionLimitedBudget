#!/bin/bash

# List of files
files=("scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_16.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_17.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_18.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_19.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_20.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_21.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_22.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_23.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_24.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_25.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_26.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_27.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_28.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_29.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_30.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_31.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_32.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_33.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_34.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_35.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_36.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_37.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_38.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_39.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_40.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_41.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_42.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_43.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_44.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_45.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_46.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_47.RDS"   
"scratch/oscDetectPaper/figs_for_paper/support_constrained_solutions/solns_spt_cvxr_drts_18_pads_0_Nmeas_48.RDS")
# Loop through each file
for file in "${files[@]}"; do
    # Create a new screen session with a unique name based on the file
    screen -d -m  bash -c "sftp beluga << EOF
    get $file
    quit
    EOF
    "
done
