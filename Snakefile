shell.executable("/bin/bash")

# on our server, the best way to run this will to navigate to the directory with 'Snakefile', and run:
# activate conda asm-pipeline-env
# 'snakemake -j 14'
#
# where -j 14 specificies that many cores to use at once as a maximum for the pipeline.  Our server is 16 cores, so this is a reasonable number.  Most steps that benefit from multithreading are set to use the provided number of cores, and won't automatically grab available cores.

# "rule medaka" may need to be modified to update the basecalling model depending on what's available.  r941_min_high_g330 means R9.4.1 flow cell, min = minion, high = high accuracy model (hac), g330 = guppy 3.3.0 or newer.# Available: r103_fast_g507, r103_fast_snp_g507, r103_fast_variant_g507, r103_hac_g507, r103_hac_snp_g507, r103_hac_variant_g507, r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360, r103_prom_snp_g3210, r103_prom_variant_g3210, r103_sup_g507, r103_sup_snp_g507, r103_sup_variant_g507, r10_min_high_g303, r10_min_high_g340, r941_min_fast_g303, r941_min_fast_g507, r941_min_fast_snp_g507, r941_min_fast_variant_g507, r941_min_hac_g507, r941_min_hac_snp_g507, r941_min_hac_variant_g507, r941_min_high_g303, r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344, r941_min_high_g351, r941_min_high_g360, r941_min_sup_g507, r941_min_sup_snp_g507, r941_min_sup_variant_g507, r941_prom_fast_g303, r941_prom_fast_g507, r941_prom_fast_snp_g507, r941_prom_fast_variant_g507, r941_prom_hac_g507, r941_prom_hac_snp_g507, r941_prom_hac_variant_g507, r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344, r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303, r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_sup_g507, r941_prom_sup_snp_g507, r941_prom_sup_variant_g507, r941_prom_variant_g303, r941_prom_variant_g322, r941_prom_variant_g360

# For example the model named r941_min_fast_g303 should be used with data from MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version 3.0.3. By contrast the model r941_prom_hac_g303 should be used with PromethION data and the high accuracy basecaller (termed "hac" in Guppy configuration files). Where a version of Guppy has been used without an exactly corresponding medaka model, the medaka model with the highest version equal to or less than the guppy version should be selected.

# If you need to add rules, FYI snakemake will automatically make directories when an output file goes into a directory that doesn't exist, but will not make directories when a directory itself is an output target (without a specific mkdir call).

# view dag output (workflow) in dot, in graphviz package (apt install graphviz): snakemake --dag | dot -Tpdf > dag.pdf

#mamba create -n medaka-test -c conda-forge -c bioconda medaka

# mamba create -n asm-pipeline-env -c conda-forge -c bioconda medaka fastqc multiqc flye circlator checkm-genome prokka filtlong

# be aware that the channel order seems to matter a lot?  conda forge must be before bioconda or medaka causes the install to fail. may be resolveable by 'conda config --set channel_priority flexible' or 'strict'


#---- Malleable variables ----#
# Set these with --config
# flyemethod = [nano-raw|nano-hq] (default: nano-raw)
# guppymodel = (default: r941_min_high_g330)
# readfmt = "fastq.gz"


#------------------------------#

#---- Fixed variables ----#
# changes here may need changes elsewhere

INPUT_DIR = "input_reads"
output_dir = "assembly"
flye_dir = os.path.join(output_dir, "02-flye_assembly")
IDS, = glob_wildcards(INPUT_DIR + "/{sample}.fastq.gz")
circlator_out_dir = os.path.join(output_dir, "04-rotated")

#--------------------------#



#---- Pipeline ----#

configfile: "config/pipeline_config.yaml"

rule all:
    input:
        "01-QC_inputs/fastqc/multiqc_report.html", # multiqc
        expand("01-QC_inputs/read_info_histograms/{sample}.txt", sample = IDS),
        "05-checkM/checkM.results.tab",
        expand("06-prokka_annotation/{sample}.polished/{sample}.gbk", sample = IDS),
        expand("07-assembly_QC/read_mapping/{sample}/{sample}.coverage.tab", sample = IDS),
        expand("07-assembly_QC/read_mapping/{sample}/{sample}.idxstats.tab", sample = IDS), 
        expand("07-assembly_QC/read_mapping/{sample}/{sample}.reads.sorted.bam", sample = IDS),
        expand("07-assembly_QC/read_mapping/{sample}/{sample}.unmapped.fastq.gz", sample = IDS)


#        expand("02-flye_assembly/{sample}/assembly.fasta", sample = IDS), # flye assembly
#        expand("01-QC_inputs/fastqc/{sample}_fastqc.html", sample = IDS), # fastqc

rule QC_reads:
    input:
        fastq = "input_reads/{sample}.fastq.gz"
    output:
        "01-QC_inputs/fastqc/{sample}_fastqc.html",
#        directory("01-QC_inputs/fastqc/")
    run:
#        shell("mkdir -p 01-QC_inputs/fastqc/")
        shell("fastqc --outdir 01-QC_inputs/fastqc {input.fastq}")


rule QC_multiqc:
    input: expand(["01-QC_inputs/fastqc/{sample}_fastqc.html"], sample = IDS)
    output: "01-QC_inputs/fastqc/multiqc_report.html"
    shell: "multiqc -o 01-QC_inputs/fastqc 01-QC_inputs/fastqc"
    


rule read_info_hists:
    input:
        fastq = "input_reads/{sample}.fastq.gz"
    output:
        file = "01-QC_inputs/read_info_histograms/{sample}.txt"
    shell:
        "bin/read_info_histograms.sh {input.fastq} > {output.file}" 

# Note: I modified this version of read_info_histograms.  It requires both histograms.py and filtlong, but assumes they are in basically the same directory or in a specific structure.  Unfortunately, read_info_histograms nor histograms.py come with the conda installation of filtlong as far as I can tell.  The modification I made was to point at a conda install of filtlong, though '$(which filtlong)'.  Also, read_info_histograms will take .gzipped or not files, and apparently will simultaniously look at illumina reads if you add them as inputs 2 and 3.
   
   
rule flye_assembly:  # Note: change '--nano-raw' to '-nano-hq' if using guppy 5+
    input:
        fastq = "input_reads/{sample}.fastq.gz"
    output:
        "02-flye_assembly/{sample}/assembly.fasta"
    threads: 
        workflow.cores
    params:
        method = expand("{param}", param = config["flyemethod"])
    shell:
        "flye --{params.method} {input.fastq} --out-dir 02-flye_assembly/{wildcards.sample} --threads {threads}" 


rule medaka:
    input: 
        assembly = "02-flye_assembly/{sample}/assembly.fasta",
        fastq = "input_reads/{sample}.fastq.gz"
    output: "03-medaka_polish/{sample}/consensus.fasta"
    threads:
        workflow.cores
    params:
        model = config["guppymodel"],
        outdir = "03-medaka_polish/{sample}"
    run:
        shell("medaka_consensus -i {input.fastq} -d {input.assembly} -o {params.outdir} -t {threads} -m {params.model}") # TO-DO: modifiy the model so it can be user provided. r941_min_high_g330


rule circlator:
    input: 
#        "02-flye_assembly/{sample}/assembly.fasta"# placeholder until I get medaka running
        "03-medaka_polish/{sample}/consensus.fasta"
    output: 
        file = "04-rotated/{sample}.assembly.fasta"
    params:
        prefix = "04-rotated/{sample}.assembly"
    shell: 
        "circlator fixstart {input} {params.prefix}"


rule checkM:
    input: 
        expand(["04-rotated/{sample}.assembly.fasta"], sample = IDS) # Critically, this is how I get snakemake to key off the results of circlator, while only running checkm once.  Had to hard-code the target directory but that's ok.
    output: # had to hard-code these values as well, doing the directory as an output 
#        dir = directory("05-checkM"),
        file = "05-checkM/checkM.results.tab"
    threads: 
        workflow.cores
    params:
        genomes = "04-rotated",
        checkm_output = "05-checkM",
        name = "checkM.results.tab"
    shell: 
        "checkm lineage_wf -t {threads} --file {params.checkm_output}/{params.name} --tab_table -x .fasta {params.genomes} {params.checkm_output}"


rule prokka:
    input:
        polished = "04-rotated/{sample}.assembly.fasta" # Medaka output, placeholder until medaka works
    output:
        polished = "06-prokka_annotation/{sample}.polished/{sample}.gbk" # medaka
    params:
        polished_outdir = "06-prokka_annotation/{sample}.polished/"
    threads:
        workflow.cores
    run:
        shell("prokka --outdir {params.polished_outdir} --force --cpus {threads} --locustag {wildcards.sample} --prefix {wildcards.sample} {input.polished}")


# add rule for mapping reads against the assembly for QC:

    
rule map_to_assembly:
    input:
        assembly = "04-rotated/{sample}.assembly.fasta",
        reads = "input_reads/{sample}.fastq.gz"
    output:
        samfile = "07-assembly_QC/read_mapping/{sample}/{sample}.aln.sam"
    threads:
        workflow.cores
    run:
        shell("bwa index {input.assembly}"),
        shell("minimap2 -ax map-ont -t {threads} {input.assembly} {input.reads} > {output.samfile}"),


rule samtools_stats:
    input: 
        samfile = "07-assembly_QC/read_mapping/{sample}/{sample}.aln.sam",
    output:
        bamfile = "07-assembly_QC/read_mapping/{sample}/{sample}.reads.sorted.bam",
        idxstats = "07-assembly_QC/read_mapping/{sample}/{sample}.idxstats.tab", 
        coverage = "07-assembly_QC/read_mapping/{sample}/{sample}.coverage.tab",
        
    threads:
        workflow.cores
    run:
        shell("samtools sort -@ {threads} -o {output.bamfile} -T reads.tmp {input.samfile}"),
        shell("samtools index -@ {threads} {output.bamfile}"),
        shell("samtools idxstats {output.bamfile} > {output.idxstats}"),
        shell("samtools coverage {output.bamfile} > {output.coverage}")

# idxstats prints: contig contig_length #mapped #unmapped


rule samtools_get_unmapped:
    input: 
        bamfile = "07-assembly_QC/read_mapping/{sample}/{sample}.reads.sorted.bam"
    output:
        bamfile = "07-assembly_QC/read_mapping/{sample}/{sample}.unmapped.bam",
        sorted_bamfile = "07-assembly_QC/read_mapping/{sample}/{sample}.unmapped.sorted.bam",
        fastq = "07-assembly_QC/read_mapping/{sample}/{sample}.unmapped.fastq.gz"
    threads:
        workflow.cores
    run:
        shell("samtools view -h -f 4 {input.bamfile} > {output.bamfile}"),
        shell("samtools sort -n -o {output.sorted_bamfile} {output.bamfile}"),
        shell("samtools fastq {output.bamfile} | gzip > {output.fastq}")
























# I can't add the ideel code here, because it requires a diamond database that would need to be independently set up, which may be too much of an ask.  Also can only have one snakemake pipeline active at once.







# IDS, = glob_wildcards("input_reads/{sample}.fastq*")
