# LongReadPipeline

## Pipeline is Currently a work in Progress 

### Current Features 

1. Pipeline appears to run to completion wihout error up to binning.
2. Pipeline allows you to download all containers using the env_initialization.py script in /scripts.

### Future Improvements 
1. Better sample loading (this has some quirks with snakemake).
2. Additional annotation and analysis steps.  
3. Better quality reports.
4. Graphs and other forms of data summarization.
5. Relevant Configuration parameters.

## Initialization steps
1. Install conda if you havent already.
2. create a conda environment with the most recent version of snakemake
3. pull this directory from github.
4. run the env_initialization.py script from the scripts folder of this program. 
5. download any datasets you would like to use/test in the samples folder. 
6. update the samples.csv file with the names of your samples as such {samples}.fasta (e.g. do not include the fasta in the name). 
    -> The location information is less important for now. I had been planning on using sample name as a key for the location to pull the 
       correct files. However, snakemake doesn't seem to like me doing that so I am looking for alternatives for easy sample loading and
       reference in the pipeline. 
8. Try running the pipeline.
9. Tell Kurtis when something inevitably doesn't work so he can update the documentation. 
