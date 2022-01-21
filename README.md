# EMoMiS

EMoMiS is an Epitope-based Molecular Mimicry Search pipeline. 
As a first step, antigens extracted from the Structural Antibody Database (SAbDab) are searched for sequence regions and structural similarity with the target protein. 
Then, a pre-trained deep learning model is used to evaluate if antibodies, known to recognize the SAbDab antigens, can cross-react with the target structure.

For a detailed description of the pipeline process, refer to our paper:

Stebliankin, Baral, Nunez-Castilla, Sobhan, Cickovski, Mondal, Siltberg-Liberles, Chapagain, Mathee, and Narasimhan (2022), EMoMiS: A Pipeline for Epitope-based Molecular Mimicry Search in Protein Structures with Applications to SARS-CoV-2, (Under review, ISMB 2022).

## Inatallation

There are multiple ways to set up an EMoMiS running environment.

### 1. Download singularity container

The best way to install EMoMiS is to download a pre-built singularity container from the [Chamelion Cloud](https://www.chameleoncloud.org/) server:

    wget -P ./env/ https://chi.tacc.chameleoncloud.org:7480/swift/v1/singularity_containers/emomis.sif

The container 'emomis.sif' was built with singularity version 3.8.5.

### 2. Build singularity container

Alternatively, you have an option to build a singularity container from the definition file './env/emomis.def'.

    cd env && sudo singularity build emomis.sif emomis.def

A singularity container 'emomis.sif' should appear in ./env folder. 

### 3. Manual installation
The software requirements can be manually installed. We recommend using singularity containers (options 1-2) instead of manual installation, as the user will have to manage environment variables. Below is the list of dependencies:
* [MaSIF](https://github.com/LPDI-EPFL/masif)
* [TM-align](http://zhanglab.ccmb.med.umich.edu/TM-align)
* [DSSP v.2.3.0](https://github.com/cmbi/dssp)
* [Blast v.2.12.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [pandas](https://pandas.pydata.org/)
* [tqdm](https://github.com/tqdm/tqdm)

The EMoMiS software does not require compilation and can be used directly when all requirements are met.
## Running EMoMiS
EMoMiS pipeline can be run in forward and reverse phases.
To reproduce the results from the paper on the molecular mimicry search with SARS-CoV-2 Spike applications, run the following:
* Forward phase: 


    cd ./experiments/demo_forward
    chmod +x run_emomis.sh
    ./run_emomis.sh
    python get_unique_hits.py
  
* Reverse phase:


    cd ./experiments/demo_forward
    chmod +x run_emomis.sh
    ./run_emomis.sh
    python get_unique_hits.py

## Configuring EMoMiS
The [default configuration file](https://github.com/stebliankin/emomis/src/config_default.py) used in the demo was set for molecular mimicry search for the SARS-CoV-2 Spike.
The user may create a custom configuration to adapt the pipeline for other organisms or modify the standard parameters.

### Input files
The first step is to download the summary file of the [SAbDab database](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/search/?all=true#downloads), which contains a list of known antibody-antigen complexes.
The path to the SAbDab summary file should be specified in the config file:

    config['input']['sabdab_summary'] = '/path/to/sabdab_summary_all.tsv'

Then, identify the scientific name for the query protein. Use the same name convention as specified in SAbDab:

    config['input']['target_name'] = 'severe acute respiratory syndrome coronavirus2'

You have the option to exclude homologous structures from the database. The names of the homologous organisms should be written exactly as they appear in the SAbDab summary file:

    config['input']['db_exclude'] = ('severe acute respiratory syndrome', 'sars coronavirus', 'middle east respiratory')

As an optional parameter, you may impose the minimum protein sequence length of the target structure. 
For example, the minimum length of 900 can be used for Spike protein to include only full Spike structures and exclude isolated regions such as RBD or peptides:

    config['input']['min_seq_len_target'] = 900

Finally, specify the location of the output files:

    config['dirs']['output'] = os.getcwd() + '/output/'

Intermediate output from the sequence similarity search (EMoMiS step A):
    
    config['out_files']['raw_blast_hits'] = config['dirs']['output'] + '1-raw_blast_hits.tsv'

Intermediate output from the filtering Blast filtering process (EMoMiS step B):

    config['out_files']['filtered_blast_hits'] = config['dirs']['output'] + '2-filtered_blast_hits.tsv'

Intermediate output from the structural similarity search (EMoMiS step C):

    config['out_files']['TM_align_all'] = config['dirs']['output'] + '3-TM_align.tsv'
    config['out_files']['TM_align_filtered'] = config['dirs']['output'] + '3-TM_align_filtered.tsv' # filtered by p-value

The final output from the deep learning binding prediction (EMoMiS step D):

    config['out_files']['MaSIF_scores'] = config['dirs']['output'] + '4-DL_scores.tsv'
    config['out_files']['MaSIF_scores_filtered'] = config['dirs']['output'] + '4-DL_scores_filtered.tsv' # Filtered by p-value

Other intermediate files and parameters can be directly copied from the [default_config.py](https://github.com/stebliankin/emomis/src/config_default.py) or modified by the user.

The output with the deep learning scores can be merged with the SAbDab metadata for better interpretation 
(see details in [get_unique_hits.py](https://github.com/stebliankin/emomis/demo_forward/get_unique_hits.py)).

Finally, run the EMoMiS pipeline with `--config` parameter:

    singularity exec $EMoMiS_PATH/env/emomis.sif python $EMoMiS_PATH/src/EMoMiS.py --config config.py

As a default, the pipeline takes as input only the name of the query protein and automatically extracts all the structures from the database.
To provide custom PDB structures to the input, manually update the list of target query structures `config['target_list']`, and use the option `--skip_sabdab` to avoid automatic downloading:

    singularity exec $EMoMiS_PATH/env/emomis.sif python $EMoMiS_PATH/src/EMoMiS.py --config config.py --skip_sabdab


## Output
The final output of the pipeline is the file specified in `config['out_files']['MaSIF_scores_filtered']`, 
which contains predicted molecular mimicry epitopes that passed sequence and structural similarity, as well as the deep learning antibody binding score thresholds.
The table may contain redundant structures because a single organism may have several reference PDB IDs.
Thus, it is recommended to remove duplicates by selecting the best score among similar structures (see example in [get_unique_hits.py](https://github.com/stebliankin/emomis/demo_forward/get_unique_hits.py)).

The output folder also contains intermediate files resulted from each step of the pipeline.

Each output file contains the following columns:
* `PDB_target` - The PDB ID of the target query protein
* `PDB_subject` - The PDB ID of the database protein
* `query_target` - The amino acid query of the target protein identified by the blast-search
* `match` - The amino acid match string between the target and database proteins identified by the blast-search
* `subject` - The amino acid query of the database protein identified by the blast-search
* `eval` - The e-value of the blast alignment
* `filter1_flag` - True if the length of the sequence alignment is more than 3
* `filter2_flag` - True if proteins are surface accessible in the sequence similarity regions
* `filter3_flag` - True if the database protein is in contact with antibody in the sequence similarity region
* `hit_i` - The position of the antibody contact residue compared to the similarity sequence
* `contact_subject_ag` - The known antibody contact residue of the database antigen
* `contact_subject_ab` - The contact residue of the antibody native to the database antigen
* `target_resid_ag` - Predicted antibody contact residue in the target query protein
* `target_resid_ab` - Contact residue of the antibody native to the target protein if available.
* `cross_ab_avail` - True if both antibodies from database and target protein complexes are available
* `len_exact_match` - The length of the exact sequence similarity match
* `RMSD_motifs` - Root Mean Squared Error (RMSD) value of structurally aligned motifs
* `TM_score` - The TM-score of the motif structural alignment
* `RMSD_pval` - The probability of belonging RMSD value to the random distribution
* `RMSD_ZSCORE` - The Z-score of the RMSD value
* `score_subject_subject` - The deep learning score of binding the database antigen with its native antibody
* `score_target_subject` - The deep learning score of the antibody from the database complex to bind the target antigen
* `pval_target_subject` - The probability of belonging `score_target_subject` to the distribution of non-binders
* `zscore_target_subject` - The Z-score of the `score_target_subject`

We note that if an antigen contains multiple chains, our pipeline pre-processes those chains separately to account for memory efficiency. 
As a result, amino acids in the binding interface of several antigens will be falsely recognized by EMoMiS as surface accessible.
Thus, each final molecular mimicry epitope prediction should be manually verified on surface accessibility.
This limitation will be fixed in future EMoMiS versions.



