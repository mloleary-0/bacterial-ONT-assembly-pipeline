# bacterial-ONT-assembly-pipeline
This is a simple snakemake pipeline to assemble and QC bacterial genomes from oxford nanopore reads.


Currently this pipeline requires a manually installed conda environment named asm-pipeline-env, described by environment.yml.  Install this environment with either conda/mamba (mamba is significantly faster) by running:

'''
conda env create -f environment.yml
'''

Then, activate the environment:

'''
conda activate asm-pipeline-env
'''

Finally, launch the pipeline with N cores:

'''
snakemake -j N
'''




Input requirements:
This pipeline currently looks for reads in fastq.gz format in a folder 'input_reads' in the current directory.  It will create several directories:

01-QC_inputs
02-flye_assembly
03-medaka_polish
04-rotated
05-checkM
06-prokka_annotation
07-assembly_QC

Results for each *.fastq.gz file will be in an individual folder, except for the results in 04-rotated (which is an intermediate, the rotated assembly can be found in 06-prokka_annotation) and 05-checkM (which is a summary of all polished, rotated assemblies)


Finally, right now this requires slight tinkering depending on the type of reads you feed it.  The medaka call must be adjusted to the proper basecalling model (a variable is provided at the top of the snakemake pipeline to make this easy.  Second, if using guppy 5+ for basecalling, the flye call should be changed from '--nano-raw' to '--nano-hq'.
