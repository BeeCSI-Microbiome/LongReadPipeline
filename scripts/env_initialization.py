# script for initializeing the environment variables

import os

flye2 = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1 -P ../rules/envs/"
medaka = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/medaka:1.6.0--py38h84d2cc8_0 -P ../rules/envs/"
metabat2 = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/metabat2:2.15--h986a166_1 -P ../rules/envs/"
minimap2 = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/minimap2:2.9--1 -P ../rules/envs/ "
porechop = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39hc16433a_3 -P ../rules/envs/"
racon = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/racon:1.5.0--h7ff8a90_0 -P ../rules/envs/"
samtools = "wget --continue --no-check-certificate https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8 -P ../rules/envs/"


os.system(flye2)
os.system(medaka)
os.system(metabat2)
os.system(minimap2)
os.system(porechop)
os.system(racon)
os.system(samtools)
