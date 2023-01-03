##Code to submit 100 jobs parallel in a supercomputer cluster 
#!/bin/bash -l                                                                                                                                                                                              
#SBATCH --job-name=AME20CS_calcs                                                                                                                                                                            
#SBATCH --account=Project_2004116                                                                                                                                                                           
#SBATCH --time=48:00:00                                                                                                                                                                                     
#SBATCH --mem-per-cpu=2G                                                                                                                                                                                    
#SBATCH --partition=large                                                                                                                                                                                   
#SBATCH --ntasks=2                                                                                                                                                                                          
#SBATCH --array=1-100                                                                                                                                                                                       
work_dir=/scratch/project_2004116/network-master/Maximes_AG_proposal/up
y=${SLURM_ARRAY_TASK_ID}
for ((y=1;y<=100;y++));do
        cd $work_dir
        mkdir up_${y}
        mkdir up_${y}/network_data
        cp network_data/* up_${y}/network_data/
        cd up_${y}
        cp ../network .
        cp ../../msu_traj_${y}.dat .
        cp ../../control_sample .
        sed -i '$ d' control_sample
        echo "msu_traj_${y}.dat" >> control_sample
        cp ../../submit_stad.sh .
        cd /scratch/project_2004116/network-master/Maximes_AG_proposal/up/up_${y}
        ./network control_sample
        rm -r network_data

        cd
done
