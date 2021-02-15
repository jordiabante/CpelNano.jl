#####################################################################################################
# CpG groups
#####################################################################################################
using FASTX
using CpelNano

# Single CpG site
CpelNano.get_grps_from_cgs([10],10)

# Two CpG sites
CpelNano.get_grps_from_cgs([10,15],10)
CpelNano.get_grps_from_cgs([10,25],10)

# Two groups
CpelNano.get_grps_from_cgs([10,30,50,70],10)
CpelNano.get_grps_from_cgs([10,18,26,50,58,66],10)

# Different min distances
CpelNano.get_grps_from_cgs([10,18,26,50,58,66],4)

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/"
fasta = "$(dataDir)/reference/hg38_chr17_43023997_43145780.fa"
reg_data = CpelNano.RegStruct(); reg_data.chrst = 1; reg_data.chrend = 3000; reg_data.chr = "hg38_dna";

# Get FASTA record for chr
fa_reader = open(FASTA.Reader, fasta, index=fasta * ".fai")
fa_rec = fa_reader["hg38_dna"]
close(fa_reader)

# Get group info
CpelNano.get_grp_info!(reg_data,fa_rec,1)
reg_data.cpg_pos
reg_data.cpg_grps
reg_data.L

#####################################################################################################
# ANALYSIS REGIONS (aka SUBREGIONS)
#####################################################################################################
## Function used in data analysis

using CpelNano
config = CpelNano.CpelNanoConfig()
config.max_size_subreg = 121
rs = CpelNano.RegStruct()
rs.chrst = 1
rs.chrend = 1007
rs.cpg_pos = [50,100,112,150,250,600,700,900]
CpelNano.get_nls_reg_info(rs,config)

## Table from fasta
using CpelNano

# IO
fasta = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/reference/chr22.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig()
config.out_dir = "/Users/jordiabante/Desktop/"
config.size_est_reg = 3000
config.max_size_subreg = 250

# Get histogram of size
CpelNano.subregion_table(fasta,config)

#####################################################################################################
# GENOME METRICS
#####################################################################################################
using CpelNano

# CpG density
dna_seq = "AACGTAACGAA"
cpg_pos = map(x -> getfield(x, :offset), eachmatch(r"[Cc][Gg]", dna_seq))
CpelNano.get_ρn(dna_seq,cpg_pos;width=3)

#####################################################################################################
# Produce genome info table
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/reference/"
fasta = "$(dataDir)/hg38_chr17_43023997_43145780.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig(); config.out_dir = "/Users/jordiabante/Desktop/"

# Get histogram of size
CpelNano.estimation_region_table(fasta,config)

#####################################################################################################
# Produce estimation region genomic info
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/reference/"
fasta = "$(dataDir)/hg38_chr17_43023997_43145780.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig(); config.out_dir = "/Users/jordiabante/Desktop/"

# Get histogram of size
CpelNano.get_estimation_region_info(fasta,"hg38_dna",1000:2000,config)

#####################################################################################################
# Produce histogram of group lengths
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/reference/"
fasta = "$(dataDir)/hg38_chr17_43023997_43145780.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig(); config.out_dir = "/Users/jordiabante/Desktop/"

# Get histogram of size
CpelNano.hist_group_length(fasta,config)

#####################################################################################################
# Produce histogram of group num of CG
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/reference/"
fasta = "$(dataDir)/hg38_chr17_43023997_43145780.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig(); config.out_dir = "/Users/jordiabante/Desktop/"

# Get histogram of size
CpelNano.hist_group_cpg_num(fasta,config)

#####################################################################################################
# Partition of chromosome (informME)
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/reference/"
fasta = "$(dataDir)/hg38_chr17_43023997_43145780.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig(); config.informme_mode = true;

# Get partition respecting CpG groups
chr_names, chr_sizes = CpelNano.get_genome_info(fasta)
chr_part = CpelNano.get_informME_chr_part(chr_sizes[1], config)

#####################################################################################################
# Partition of chromosome
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/full_example/reference/"
fasta = "$(dataDir)/hg38_chr17_43023997_43145780.fa"

# Partition evenly spaced CpG sites
config = CpelNano.CpelNanoConfig(); config.min_grp_dist = 1; config.verbose = true;

# Get partition respecting CpG groups
chr = "hg38_dna"
chr_part = CpelNano.get_chr_part(chr, config, fasta)

#####################################################################################################
# Partition in targeted case
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/targeted_example/"
fasta = "$(dataDir)/reference/hg38_chr17_42879379_43290399.fa"
bed = "$(dataDir)/regions_of_interest/brca1.bed"

# IO
config = CpelNano.CpelNanoConfig(); config.bed_reg = bed; config.min_grp_dist = 1;

# Get dictionary
targ_regs = CpelNano.read_bed_reg(config)

# Get chr partition
chr = "hg38_dna"
chr_part = CpelNano.get_chr_part(chr, config, fasta, targ_regs[chr])
