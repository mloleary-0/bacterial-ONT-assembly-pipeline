# bacterial-ONT-assembly-pipeline
This is a simple snakemake pipeline to assemble and QC bacterial genomes from oxford nanopore reads using flye v2.9.1. \
  \
  \
This pipeline requires a manually installed conda environment named asm-pipeline-env, described by environment.yml - it will not install the environment for you.  Install this environment with either conda/mamba (mamba is significantly faster) by running:


```
conda env create -f config/environment.yml
```

Then, activate the environment:

```
conda activate asm-pipeline-env
```

Finally, move your sequencing reads (.fastq.gz) to a folder in this directory called 'input_reads', and launch the pipeline with N cores:

```
snakemake -j N
```


 \
 \
 **Modifiable Parameters**  
 
  
 
Medaka's 'medaka_consensus' program is used for polishing the flye assembly.  To obtain the best results, Medaka needs to know which model to use to correct assembly errors.  These models are based off a combination of flow cell, platform, and guppy version (e.g., "flow cell"_"sequencer"_"guppy basecall model + version"). **The default model for this pipeline is "r941_min_high_g330"**, which is appropriate for reads from an r 9.4.1 flow cell run on a MinIon using guppy_basecaller v3.3.0+ using the high accuracy model.  See the Medaka github for more details.  The 'medaka_consensus -h" command will print a list of available models - pick the newest one that is _not newer than the basecaller version and method you used___.  The default model can be overriden to use another model (in this example, r941_min_hac_g507):


```
snakemake -j N --config guppymodel=r941_min_hac_g507
```
  \
  \
By default, this pipeline will use the --nano-raw flag for flye.  However, if you are using corrected long reads, or reads generated by Guppy v5+ or newer, you can use flye's --nano-hq mode by instead adding the '--config flyemethod' parameter as shown below.  This mode will allow flye to take advantage of the higher accuracy in newer basecalling models. 
  
  
```
snakemake -j N --config flyemethod=nano-hq
```
  \
These can of course be combined to produce a better assembly from newer data:
```
snakemake -j N --config flyemethod=nano-hq guppymodel=r941_min_hac_g507
```
  \
Both of these values could also changed by editing the respective values in config/pipeline_config.yaml.  \
  \
  \
**Input requirements:**

This pipeline currently looks for reads in fastq.gz format in a folder 'input_reads' in the current directory.  It will create several directories:

01-QC_inputs  \
02-flye_assembly  \
03-medaka_polish  \
04-rotated  \
05-checkM  \
06-prokka_annotation  \
07-assembly_QC  
  

Results for each *.fastq.gz file will be in an individual folder, except for the results in 04-rotated (which is an intermediate, the rotated assembly can be found in 06-prokka_annotation) and 05-checkM (which is a summary of all polished, rotated assemblies)
