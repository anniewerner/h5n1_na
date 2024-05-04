# h5n1_na
h5n1 NA analysis tools for avian, non-human mammalian, and human sequences
## Downloading data and initial processing

The N1 neuraminidase (NA) data that is used in these analyses are directly available on [GISAID](https://gisaid.org/). The first dataset contains H5N1 NA sequence information with the following options:

<img width="700" alt="Screenshot 2024-05-03 at 4 23 08 PM" src="https://github.com/anniewerner/h5n1_na/assets/149183888/6e665c5d-569f-4c85-ad2d-96c38b3ffaf0">

Download data in .fasta format and rename appropriately, for example, `h5n1_na_raw.fasta` It is wise to save your .fasta files in a folder that indicates the influenza type (i.e., a_h5n1) and/or host type (i.e., a_h5n1_nonhuman_mammal).

### Data cleaning and processing

#### Removing duplicated names or sequences from the initial .fasta download.

Dependencies: [seqkit](https://bioinf.shenwei.me/seqkit/download/), a command-line tool that helps clean .fasta and .fastq files as indicated. Please refer to their documentation for installing the appropriate and most up-to-date version.

The first data processing we will do is removal of duplicates by name and by sequence. Keep in mind that if you plan to run a nextstrain build on this set of sequences, **removing sequences with duplicate names is necessary for nextstrain** to function properly, otherwise it will produce an error due to duplicate names\*\*.

\*\**If you have your own nextstrain build that does not include a parse rule in the snakefile, or if your parse rule does not alter the information in the header of each sequence in your file, then this may not apply.*

Removing sequences that are identical is not necessary, but is a good initial data down-sampling if you are working with more than a few hundred sequences.

1.  copy-and-paste the following into your CLI to form and enter your working directory, \<your_directory\>. Also, you will create a directory named `nextclade_output` and/or `augur_alignment_output` within `your_directory`, depending on which method you will use for sequence alignment.

    ``` bash
    mkdir <your_directory>
    cd <your_directory>
    mkdir <nextclade_output> # create if you are aligning via nextclade CLI
    mkdir <augur_alignment_output> # create if you are aligning via augur align
    ```

2.  to remove sequences based on the content of the sequence, use the -s flag as below (**optional**):

    ``` bash
    ~/<your_directory> seqkit rmdup h5n1_na_raw.fasta -s > ./h5n1_na_noseqdups.fasta
    # removing duplicate sequenceces by sequence identity
    ```

3.  to remove sequences based on the content of the name, use the -b flag as below, *on the file from which you have already removed identical sequences,* `h5n1_na_noseqdups.fasta`.

    ``` bash
    ~/<your_directory> seqkit rmdup h5n1_na_noseqdups.fasta -n > ./h5n1_na.fasta
    # removing duplicate sequences by sequence name or sequence header
    ```

Now, we can work with a duplicate-free set of sequences that has already been initially down-sampled.

#### Multiple sequence alignment (MSA) of downloaded sequences

For multiple sequence alignment, this tutorial details two methods that you can use: [nextclade](https://clades.nextstrain.org/) ([nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/)) or [augur](https://docs.nextstrain.org/projects/augur/en/stable/) (part of the [nextstrain CLI](https://docs.nextstrain.org/en/latest/install.html)). These applications have their own advantages and disadvantages, but are commonly used for influenza virus genome sequences and therefore provide consistency between your findings and those reported by others.

Augur is part of all nextstrain runtimes and does not require an additional installation if you already have nextstrain installed. Augur uses the MAFFT algorithm ([Katoh et al., 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC135756/)) for rapid large-scale multiple sequence alignments. MAFFT can be downloaded locally ([documentation here](https://mafft.cbrc.jp/alignment/software/)) and run on its own as well, instead of through the nextstrain runtime.

Nextclade CLI should be used if you have upwards of a few hundred sequences in your file that you wish to align, as the [web-based application](https://clades.nextstrain.org/) does not handle large files as well. Nextclade aligns sequences through a *local sequence alignment* algorithm known as the Smith-Waterman algorithm.

Please refer to the (very well annotated and detailed) documentation for nextclade CLI and/or augur, both of which are linked above, to decide which you will use.

##### Using Nextclade CLI {data-link="Using Nextclade CLI"}

1.  Because standardized reference datasets are not available for H5N1 sequences as of yet, you will have to create your own or use the one available in this repository. These data are adapted from the h1n1pdm09_na nextclade dataset and from the avian-flu nextstrain build, and can be cloned by running the following code:

    ``` bash
    git clone https://github.com/h5n1_na/reference_dataset.git
    # downloading reference dataset for h5n1_na
    ```

    This should download the following files into the specified directory: `genome_annotation.gff3`, `pathogen.json,` `reference.fasta`, `reference_h5n1_na.gb`. You can also copy-and-paste your cleaned sequence file, named in the previous section `h5n1_na.fasta,` here and change the name to `sequences.fasta`. You can alternatively keep the file containing your sequences in its directory and include the path and name when running nextclade.

2.  Run nextclade by running the following code:

    ``` bash
    nextclade run h5n1_na.fasta \
    --input-ref reference_dataset/reference.fasta \
    --input-pathogen-json reference_dataset/pathogen.json \
    --input-annotation reference_dataset/genome_annotation.gff3 \
    --output-all nextclade_output 
    ```

    alternatively, if you have copied-and-pasted your `h5n1_na.fasta` file into the reference dataset folder, you can simply write:

    ``` bash
    nextclade run --input-dataset reference_dataset \
    --output-all nextclade_output
    ```

    This creates the following files in the output directory named `nextclade_output`:

    `nextclade.tsv`

    `nextclade.csv`

    `nextclade.ndjson`

    `nextclade.json`

    `nextclade.cds_translation.UFK99042.1.fasta`

    `nextclade.aligned.fasta`\*

    \*if you have provided a sequence to nextclade that has already been aligned, it will not give you this file.

##### Using Nextstrain's Augur through the nextstrain runtime

1.  Run augur by **first entering the nextstrain runtime**. For any sub-command of nextstrain, this tutorial follows the conventions for using a **docker runtime**. This assumes you have previously installed and confirmed that nextstrain can run locally on your device via CLI-given commands.

    To get started, run the following:

    ``` bash
    nextstrain shell . # creating the nextstrain shell.

    Nextstrain  ~/build $ # entering the runtime.

    Nextstrain  ~/build $ cd/<your_directory>

    Nextstrain  ~/build/<your_directory> $ augur --help
    # confirming augur is installed and functioning
    # by displaying all possible flags and subcommands of augur.
    ```

2.  Once in the appropriate directory within the nextstrain shell, run the following augur command to align your data.

    ``` bash
    Nextstrain ~/build/<your_directory> $ augur align
    --sequences <your_directory>/h5n1_na.fasta \
    --reference-sequence <your_directory>/reference_dataset/reference_h5n1_na.gb \
    --fill-gaps \
    --output <your_directory>/augur_alignment_output/alignment.fasta
    ```

    1.  Please note that depending on the size of the `.fasta` file you provide, the alignment may take more than 30 minutes (this is for particularly large numbers of sequences).

    2.  This can be changed by adding the `--nthreads` flag and specifying how many threads you want to be used.

3.  This creates the following files in the output directory `augur_alignment_output`:

    `alignment.fasta.log`

    `alignment.fasta.insertions.csv`

    `alignment.fasta`

4.  You can now use either the `alignment.fasta` file created by running augur as outlined above the `nextclade.aligned.fasta` file created by [Using Nextclade CLI] for further sequence analyses!

Note: if you are using an application like [geneious prime](https://www.geneious.com/) for your sequence analyses, you might notice that it has its own multiple sequence alignment (MSA) tools. You are welcome to use these if you wish, but note that your sequence analyses will yield different results depending on which MSA algorithm you use. You are encouraged to research all the available MSA algorithms at your disposal, and to select one that you feel is most appropriate for the analyses you plan to do.
