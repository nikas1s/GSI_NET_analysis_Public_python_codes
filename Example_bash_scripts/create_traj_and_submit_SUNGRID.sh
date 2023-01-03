#!/bin/bash                                                                                                                                                                                                 
#SBATCH --job-name=HFB_reactionrate_calcs                                                                                                                                                                   
#SBATCH --account=Project_2004116                                                                                                                                                                           
#SBATCH --time=48:00:00                                                                                                                                                                                     
#SBATCH --mem-per-cpu=2G                                                                                                                                                                                    
#SBATCH --partition=small                                                                                                                                                                                   
##SBATCH --mail-type=BEGIN #uncomment to enable mail                                                                                                                                                        
module load python-data
##Define The directories and files changes here only                                                                                                                                                        
talys_dir=/users/stynikas/bin/talys
run_dir=HFB_no_exp
work_dir=/projappl/project_2004116/HFB_reaction_rates_Talys_Dec_2022/
exec_name=Access_and_run_masses_HFB.py
##Do not change from here onwards                                                                                                                                                                           

mkdir ${work_dir}/${run_dir}_$1
cp $talys_dir ${work_dir}/${run_dir}_$1/talys
cp $work_dir/$exec_name ${work_dir}/${run_dir}_$1/$exec_name

cd ${work_dir}/${run_dir}_$1
python ${work_dir}/${run_dir}_$1/$exec_name $1
