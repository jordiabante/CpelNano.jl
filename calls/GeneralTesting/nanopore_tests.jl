#####################################################################################################
# Partition nanolish file
#####################################################################################################
using CpelNano

## Regular mode
pkg_dir = dirname(dirname(pathof(CpelNano)))

# IO
dataDir = "$(pkg_dir)/examples/full_example/"
nano = "$(dataDir)/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"

# Partition nanopolish file into 5 files
CpelNano.split_nanopolish_file(nano,5)

#####################################################################################################
# Get model dictionary nanopolish model file
#####################################################################################################
using CpelNano

# IO
dataDir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/Nanopolish-Pore-Model/"
nano = "$(dataDir)/r9.4_450bps.cpg.6mer.template.model"

# Get nanopolish model dictionary
nano_dic = CpelNano.get_nano_pore_mod(nano)

#####################################################################################################
# Sample averages nano file
#####################################################################################################
using CpelNano

## Regular mode
pkg_dir = dirname(dirname(pathof(CpelNano)))

# IO
dataDir = "$(pkg_dir)/examples/full_example/"
nano = "$(dataDir)/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/reference/hg38_chr17_43023997_43145780.fa"
outdir = "$(dataDir)/cpelnano/"
outdir = "/Users/jordiabante/Desktop/"

# IO
config = CpelNano.CpelNanoConfig(); 
config.out_prefix = "full_example_sample";
config.out_dir = outdir;
config.verbose = false; # true;

# Analyze nano
CpelNano.nano_smp_est(nano,fasta,config)

#####################################################################################################
# Process nano file
#####################################################################################################
using CpelNano

## Regular mode
pkg_dir = dirname(dirname(pathof(CpelNano)))

# IO
dataDir = "$(pkg_dir)/examples/full_example/"
nano = "$(dataDir)/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/reference/hg38_chr17_43023997_43145780.fa"
outdir = "$(dataDir)/cpelnano/"

# IO
config = CpelNano.CpelNanoConfig(); 
config.mod_nse = true;
config.verbose = false; # true;
config.out_dir = outdir;
config.out_prefix = "full_example_mod_nse_$(config.mod_nse)";

# Analyze nano
CpelNano.analyze_nano(nano,fasta,config)

## Targeted mode

# IO
dataDir = "$(pkg_dir)/examples/targeted_example/"
nano = "$(dataDir)/nanopolish/targeted_example_noise2.0_sorted_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/reference/hg38_chr17_42879379_43290399.fa"
bed = "$(dataDir)/regions_of_interest/brca1.bed"
outdir = "$(dataDir)/cpelnano/"
outdir = "/Users/jordiabante/Desktop/"

# IO
config = CpelNano.CpelNanoConfig(); 
config.out_prefix = "targeted_example";
config.out_dir = outdir;
config.verbose = true;
config.bed_reg = bed; 

# Analyze nano
CpelNano.analyze_nano(nano,fasta,config)

#####################################################################################################
# ESTIMATION
#####################################################################################################

## Process chromosome
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/"
nano = "$(dataDir)/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/reference/hg38_chr17_43023997_43145780.fa"
outdir = "$(dataDir)/cpelnano/"

# Configuration
min_cov = 10.0; max_size_subreg = 500; size_est_reg = 4000; max_em_init = 10; max_em_iters = 20; 
config = CpelNano.CpelNanoConfig(min_cov, max_size_subreg, size_est_reg, max_em_init, max_em_iters);
config.out_dir = outdir; config.out_prefix = "test";

# Analyze each region in chr
CpelNano.analyze_nano(nano,fasta,config)

# Output files 
mml_path = "$(dataDir)/cpelnano/test_mml.bedGraph"
nme_path = "$(dataDir)/cpelnano/test_nme.bedGraph"
theta_path = "$(dataDir)/cpelnano/test_theta.bedGraph"

# Print files
# println(String(read(mml_path)))
# println(String(read(nme_path)))
println(String(read(theta_path)))

#####################################################################################################
# DIFFERENTIAL ANALYSIS
#####################################################################################################

# Read in model file
models = CpelNano.read_mod_file_chr(theta_path, "Reference")
CpelNano.comp_gjsd(models["Reference_1_4000"],models["Reference_1_4000"])

# Get all models
mod_files = fill(theta_path, 5)
mods = CpelNano.get_all_models_chr(mod_files, "Reference")
CpelNano.get_unique_ids(mods)
uniq_ids = unique(vcat(CpelNano.get_unique_ids(mods), CpelNano.get_unique_ids(mods)))
