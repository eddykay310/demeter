**Demeter**

[TOC]

### Table of Contents

------



### Overview

------



![image-20220712133902287](C:\Users\DELL\AppData\Roaming\Typora\typora-user-images\image-20220712133902287.png)



### Installation



### Usage

------

**Basic**

Using paired-end Illumina reads

```
bash -i path/demeter/main.sh -i path/[SRA gzipped fasta] -r PE -c 4 -m 8 -a U50 -t all
```

Using assembled whole genome

```
bash -i path/demeter/main.sh -i path/[SRA gzipped fasta] -r PE -c 4 -m 8 -a U50 -p 5 -t all
```

**Options**

```
[ -fd | --fasta ID ] 			Path to input file (sra for Illumina paired-end reads or assembled whole genome fasta file)
[ -rt | --runtype ]
[ -c | --cpus ]					Number of cpus to be used
[ -rm | --ram ]					Amount of RAM to be usd
[ -aq | --assembly_qc ]			Type of assembly QC to use
[ -ref | --reference_genome ]	 Reference genome to be used for correcting missambles and scaffolding
[ -anc | --ancestor_genome ]	Ancestor genome for calculating selection pressure (dn/ds)
[ -d | --dnds ]					dn/ds threshold to be used for determining a pseudogene
[ -ot | --organism_tag ]		Tag used for assigning unique labels during analysis (eg. gene ids during annotation)
[ -kt | --keep_transposease ]	Keep transposeases as pseudogenes 
[ -ps | --pipeline_stage ]		Stage to start pipeline
[ -s | --scaffold_genome ]		Scaffold genome (need a reference genome)
[ -pt | --pseudogene_predictor ]	Pseudogene tool to use [Prokka (prokka), DFAST (DFAST), Psedofinder (PF), All (all) ]
[ -h | --help ]			Help
```

