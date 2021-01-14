## Differential analysis: analysis region level
using CpelNano

# IO
dataDir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/toy_example/"
fasta = "$(dataDir)/reference.fa"

# Unmatched test
min_cov=5.0; max_size_subreg=1000; size_est_reg=4000; max_em_init=10; max_em_iters=20; 
config = CpelNano.CpelNanoConfig(min_cov,max_size_subreg,size_est_reg,max_em_init,max_em_iters);
config.out_dir = "/Users/jordiabante/Desktop/"; config.out_prefix = "test"; config.matched = false;
theta_path_g1 = "/Users/jordiabante/Desktop/test_theta_g1.cpelnano";
theta_path_g2 = "/Users/jordiabante/Desktop/test_theta_g2.cpelnano";
mod_files_g1 = fill(theta_path_g1,5); mod_files_g2 = fill(theta_path_g2,5);
# CpelNano.differential_analysis(mod_files_g1,mod_files_g1,fasta,config)
# CpelNano.differential_analysis(mod_files_g2,mod_files_g2,fasta,config)
CpelNano.differential_analysis(mod_files_g1,mod_files_g2,fasta,config)

# Matched test
min_cov=5.0; max_size_subreg=1000; size_est_reg=4000; max_em_init=10; max_em_iters=20; 
config = CpelNano.CpelNanoConfig(min_cov,max_size_subreg,size_est_reg,max_em_init,max_em_iters);
config.out_dir = "/Users/jordiabante/Desktop/"; config.out_prefix = "test"; config.matched = true;
theta_path_g1 = "/Users/jordiabante/Desktop/test_theta_g1.cpelnano";
theta_path_g2 = "/Users/jordiabante/Desktop/test_theta_g2.cpelnano";
mod_files_g1 = fill(theta_path_g1,8); mod_files_g2 = fill(theta_path_g2,8);
# CpelNano.differential_analysis(mod_files_g1,mod_files_g1,fasta,config)
# CpelNano.differential_analysis(mod_files_g2,mod_files_g2,fasta,config)
CpelNano.differential_analysis(mod_files_g1,mod_files_g2,fasta,config)
