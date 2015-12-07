#these are command-line functions. You might need to set the permissions
# on the functions to enable them to be executables.
# For example: in the command line, run "chmod +x main_inputReg_template.R"
#  to set the file to be an executable.
# The only thing you need to make sure is that /usr/bin/Rscript is the correct
#  location if your R interpreter.

#first create the neighbor list based on the MNI template
#nohup ./main_inputReq_template.R --template "/home/smile/klsix/fmri_script_test/20151201_test/50002/MNI152_T1_2mm_brain_symmetric.nii.gz" --output "/home/smile/klsix/fmri_script_test/20151206_felix_output/template" --pattern 7 --details "List of neighbors (size 7) for the MNI_T1_2mm_brain_symmetric file"

#extract the 2D matrix from ABIDE
nohup ./main_inputReq_convert.R --input "/home/smile/klsix/fmri_script_test/20151201_test/50002/func2mni.nii.gz" --output "/home/smile/klsix/fmri_script_test/20151206_felix_output/ABIDE_50002_matrix" --template "/home/smile/klsix/fmri_script_test/20151206_felix_output/template_2015-12-07.RData" --details "Extracted 2D matrix from ABIDE 50002 after preprocessing with C-PAC, taking preprocessing.nii.gz and using flirt to match the dimensions of the MNI_2mm template. The mask for this comes from MNI_2mm as well."
