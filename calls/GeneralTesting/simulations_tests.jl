########################################################################################
# Read size info
########################################################################################
using CpelNano

dataDir = "/Users/jordiabante/Desktop/"
fastq = "$(dataDir)/ex.fastq"
CpelNano. fastq_len_store(fastq,dataDir)
CpelNano. fastq_len_dist(fastq,dataDir)

########################################################################################
# Read size testing
########################################################################################
using Plots
using CpelNano
using Distributions

# Histogram
u = CpelNano.gen_read_len(50000);
plot(u/1000,seriestype=:histogram,bins=100,xlim=(0,205),
    xticks=[18.4,38.9,59.4,79.9,100,121,141,162,182,203],
    xlabel="Read length l (kbp)",ylabel="Counts",label="")
savefig("/Users/jordiabante/Desktop/Read-Size-Histo.png")

# Total Basecalled Bases
reso = 5000
u = Int.(round.(10000.0*rand(Gamma(2,0.5),50000)))
uniq_sizes = unique(sort(u))
weights = map(y->length(findall(x->x==y,u)),uniq_sizes)
delims = collect(1:reso:130000)
bin_id = map(u->findfirst(x->u-x<reso,delims),uniq_sizes)
uniq_bin_id = unique(bin_id) 
total_basecalled = map(i->sum(uniq_sizes[bin_id.==i].*weights[bin_id.==i]),uniq_bin_id)
delims_mid = 0.5*(delims[2:end]+delims[1:(end-1)])
plot(delims/1000,total_basecalled,seriestype=:bar,xlim=(0,205),
    xticks=[18.4,38.9,59.4,79.9,100,121,141,162,182,203],
    xlabel="Read length l (kbp)",ylabel="Total Basecalled Bases",label="")

########################################################################################
# MCMC sampling (we tried to see if it's feasible to generate entire chromosome)
########################################################################################
using Plots
using CpelNano

reps = 10
ns = [5,10,25,50,100,500,1000,2500]
total_time = zeros(length(ns))
for i=1:length(ns)
    CpelNano.print_log("n=$(ns[i])")
    for j=1:reps
        total_time[i] += @elapsed CpelNano.metro_hast_inst(ns[i],10000,zeros(ns[i]),zeros(ns[i]-1))
    end
    total_time[i] /= reps
end

plot(ns,total_time,xlabel="N",ylabel="Time (s)",label="")
plot!(ns,ns/250,xlabel="N",ylabel="Time (s)",label="")

########################################################################################
# Augmenting model file
########################################################################################
using CpelNano

## Regular mode
pkg_dir = dirname(dirname(pathof(CpelNano)))

# IO
dataDir = "$(pkg_dir)/examples/full_example/"
nano = "$(dataDir)/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/reference/hg38_chr17_43023997_43145780.fa"
outdir = "/Users/jordiabante/Desktop/"

# IO
config = CpelNano.CpelNanoConfig(); 
config.out_prefix = "full_example";
config.out_dir = outdir;
config.verbose = false; # true;

# Analyze nano
CpelNano.analyze_nano(nano,fasta,config)

# Augment file
chr = "hg38_dna"
mod_file = "$(outdir)/$(config.out_prefix)_theta.txt"
CpelNano.aug_mod_file(mod_file,chr,fasta)

########################################################################################
# Sandeep testing
########################################################################################
using CpelNano
using FASTX
using DelimitedFiles

# IO
dataDir = "/Users/sandeepk/Downloads/"
fasta = "$(dataDir)/chr22.fa"
chr22record = CpelNano.get_chr_fa_rec("chr22",fasta)
genome = convert(String,FASTX.FASTA.sequence(chr22record,35723841:35723842));
genome = convert(String,FASTA.sequence(chr22record,40688790:40688791));
model_regs = aug_mod_file("/Users/sandeepk/Downloads/GM12878_CpelNano_Theta.txt", "chr22", fasta, chr22record);

config = CpelNano.CpelNanoConfig();
chr_part = CpelNano.get_chr_part("chr22", config, fasta);
model_regs = read_reg_chr("/Users/sandeepk/Downloads/GM12878_CpelNano_Theta_run4.txt", "chr22", fasta, chr22record);
# Test constructing a methyl vector

# Double check if length of test_methylvec matches number of cpg sites in read
fa_record = chr22record
chrstart = fa_record.sequence[1]
chrend = fa_record.sequence[length(fa_record.sequence)];
read_len = 100000
genom_pos, newRead = gen_read("chr22", fasta, read_len);
startidx = linearsearch_readstart(genom_pos, chr_part);
endidx = linearsearch_readend(genom_pos+read_len-1, startidx, chr_part);
meth_lattice, out_cpg_pos = get_methyl_vector(genom_pos, genom_pos+read_len-1, startidx, endidx, model_regs);
#meth_lattice, out_cpg_pos = CpelNano.get_methyl_vector(40576278, 41576277, 13526, 13859, model_regs);
cpg_pos = map(x->getfield(x,:offset),eachmatch(r"CG",newRead));
length(meth_lattice)
length(cpg_pos)
for i in 1:length(cpg_pos)
    cpg_pos[i] = cpg_pos[i] + genom_pos - 1
end

for i in 1:length(cpg_pos)
    if cpg_pos[i] != out_cpg_pos[i]
        print(i)
        break
    end
end

dataDir = "/Users/sandeepk/Downloads/"
fasta = "$(dataDir)/chr22.fa"
config = CpelNano.CpelNanoConfig();
chr_part = CpelNano.get_chr_part("chr22", config, fasta);
chr_part[13563]
chr22record = CpelNano.get_chr_fa_rec("chr22",fasta)
genome = convert(String,FASTA.sequence(chr22record,40685791:40688791))


##### Create full GM12878 Chr22 parameter file with alpha, beta
dataDir = "/Users/sandeepk/Downloads/"
fasta = "$(dataDir)/chr22.fa"
fa_record = CpelNano.get_chr_fa_rec("chr22",fasta)
config = CpelNano.CpelNanoConfig();
chr_part = CpelNano.get_chr_part("chr22", config, fasta);

model_regs = CpelNano.aug_mod_file("/Users/sandeepk/Downloads/GM12878_CpelNano_Theta.txt", "chr22", fasta, fa_record);

output = Vector{String}()

for i in 1:length(chr_part)
    reg = model_regs[i]
    #if reg.L > 0
        ϕhat_str = string(reg.ϕhat[1], ",", reg.ϕhat[2], ",", reg.ϕhat[3])
        α,β = CpelNano.get_αβ_from_ϕ(reg.ϕhat, reg)
        α_str = join(α, ",")
        β_str = join(β, ",")
        push!(output, string(reg.chr, '\t', reg.chrst, '\t', reg.chrend, '\t', ϕhat_str, '\t', α_str, " ", β_str, "\n"))
    #end
end
output_str = join(output);

open("/Users/sandeepk/Downloads/GM12878_CpelNano_Theta_run4.txt", "w") do io
    write(io, output_str)
end;
    