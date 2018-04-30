remoVecSec
==========

remoVecSec is a library of modules that allows to remove contamination
in assembled genomes prior to NCBI WGS submission. **remower.py** is a
wrapper script removes contamination in a set of sequences

Background
----------

NCBI's Foreign Contamination Screens, November 2016

The purpose of the foreign contamination screens is to identify
contaminating sequences that may be present for artificial reasons or
for biological reasons. Artificial reasons include cloning artifacts
(vector, linker/adaptor/primer, E. coli host DNA), contamination in the
lab with human sequence, mixing of samples or sequencing runs with other
organisms, and bacterial insertion sequences that have integrated into
sequenced clones. Biological reasons include the presence of
endosymbionts, infectious agents, or microbes residing on the surface of
the organism or in the gut when the DNA prep was made.

Our suite of foreign contamination screens uses BLAST to screen the
submitted sequences against:

1.  a common contaminants database that contains vector sequences,
    bacterial insertion sequences, E. coli and phage genomes

2.  a database of adaptors linkers and primers

3.  a database of mitochondrial genomes

4.  the chromosomes of unrelated organisms

5.  a database of ribosomal RNA genes

Suspect spans are re-BLASTed against:

1.  the chromosomes of unrelated organisms

2.  the chromosomes of related organisms

3.  the NCBI nt BLAST database of nucleotide sequence from all
    traditional divisions of GenBank, EMBL, and DDBJ

4.  the NCBI htgs BLAST database of sequences from the HTG division of
    GenBank, EMBL, and DDBJ

Results similar to those obtained by NCBI could be generated by running
the screens as described below.

### Common contaminant screen

#### Databases

1.  File to screen for the common contaminants in eukaryotic sequences:

    ``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_euks.fa.gz
    ```

Contains the cloning artifacts that are likely to show up as
contaminants across all eukaryotic species: vector sequences, E.coli
genome, phage genomes, bacterial Insertion Sequences and transposons.

1.  File to screen for the common contaminants in prokaryotic sequences:

    ``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/contam_in_prok.fa
    ```

Contains phiX174.

These files need to be unzipped and the resulting FASTA sequence files
formatted as BLAST databases using the makeblastdb program.

#### Programs

blastn and makeblastdb are contained in the blast+ package which can be
installed following the instruction in the BLAST help documents.

"BLAST Command Line Applications User Manual":

``` {.example}
        https://www.ncbi.nlm.nih.gov/books/NBK279671/
```

"Standalone BLAST Setup for Windows PC":

``` {.example}
        https://www.ncbi.nlm.nih.gov/books/NBK52637/
```

"Standalone BLAST Setup for Unix":

``` {.example}
        https://www.ncbi.nlm.nih.gov/books/NBK52640/
```

#### Execution

A BLAST search is run against either the contam\_in\_euks or
contam\_in\_prok database, depending on the origin of the input
sequences. The common contaminant BLAST results are filtered for hits
over various length and percent identity cut-offs.

Command line:

1.  for screening eukaryotic sequences:

    ``` {.example}
        blastn -query _input_fasta_sequences_ -db contam_in_euks -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)'
    ```

OR with an intermediate file, these 2 commands:

``` {.example}
        blastn -query _input_fasta_sequences_ -db contam_in_euks -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out _out_file_

        awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' _out_file_
```

1.  for screening prokaryotic sequences:

    ``` {.example}
        blastn -query _input_fasta_sequences_ -db contam_in_prok -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)'
    ```

OR with an intermediate file, these 2 commands:

``` {.example}
        blastn -query _input_fasta_sequences_ -db contam_in_prok -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 90.0 -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out _out_file_

        awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' _out_file_
```

### Adaptor screen

VecScreen (<https://www.ncbi.nlm.nih.gov/tools/vecscreen/>) is run
against either the adaptors\_for\_screening\_euks.fa database or
adaptors\_for\_screening\_proks.fa database, depending on the origin of
the input sequences. Hits are filtered to retain only those matches that
VecScreen classifies as "Strong" or "Moderate" (see:
<https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Categories>).

#### Databases

The adaptors\_for\_screening databases are available here:

``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_euks.fa

        ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_proks.fa
```

These FASTA sequence files need to be formatted as BLAST databases using
the makeblastdb program.

#### Programs

The VecScreen standalone program is available here:

``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/blast/demo/vecscreen
```

The script to filter the VecScreen results is here:

``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/VSlistTo1HitPerLine.awk
```

#### Execution

Command line:

1.  for screening eukaryotic sequences:

    ``` {.example}
        vecscreen -d adaptors_for_screening_euks.fa -f3 -i _input_fasta_sequences_ -o _vs_output_file_
    ```

2.  for screening prokaryotic sequences:

    ``` {.example}
        vecscreen -d adaptors_for_screening_proks.fa -f3 -i _input_fasta_sequences_ -o _vs_output_file_
    ```

Filter out the "Weak" and "Suspect Origin" hits:

``` {.example}
        VSlistTo1HitPerLine.awk suspect=0 weak=0 _vs_output_file_ > _filtered_vs_output_file_
```

### Mitochondrial genome screen

BLAST is used to screen the input sequences against a database of the
mitochondrial genome sequences in the NCBI Reference Sequences (RefSeq)
collection.

#### Database

``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/mito.nt.gz
```

This file needs to be unzipped and the resulting FASTA sequence file
formatted as a BLAST database using the makeblastdb program.

#### Programs

blastn and makeblastdb are contained in the blast+ package (see above).

#### Execution

The BLAST hits to mitochondrial genomes are filtered for hits over 98.6%
identity and at least 120 bases long.

``` {.example}
        blastn -query _input_fasta_sequences -db mito.nt -out % -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 98.6 -soft_masking true -outfmt 7 | awk '$4>=120' > _filtered_mito_output_file_
```

### Ribosomal RNA screen

Ribosomal RNA genes are the cause of many false positives because the
include some segments that align to distantly related organisms.
Segments that match rRNA genes are identified so that such segments are
not reported as being foreign.

BLAST is used to screen the input sequences against a database of the
rRNA gene sequences .

#### Database

------------------------------------------------------------------------

``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/rrna.gz
```

This file needs to be unzipped and the resulting FASTA sequence file
formatted as a BLAST database using the makeblastdb program.

#### Programs

blastn and makeblastdb are contained in the blast+ package (see above).

#### Execution

The BLAST hits to rRNA genes are filtered for hits over 95% identity and
at least 100 bases long.

``` {.example}
        blastn -query _input_fasta_sequences_ -db rrna -task megablast -template_length 18 -template_type coding -window_size 120 -word_size 12 -xdrop_gap 20 -no_greedy -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 1E-9 -gapextend 2 -gapopen 4 -penalty -4 -perc_identity 95 -reward 3 -soft_masking true -outfmt 7 | awk '$4>=100' > _filtered_mito_output_file_
```

### Foreign chromosome screen

Screens for matches to chromosome sequences from foreign organisms.
Foreign organisms are those that belong to a different taxonomic group
compared to the organism whose sequences are being screened. The
taxonomic groups are:

arthropoda, chordata, other\_metazoa,

viridiplantae, fungi, other\_eukaryota,

bacteria, archaea, viruses\_and\_viroids

#### Databases

Our databases to detect cross-contamination detection are limited to
assemblies that have been publicly released in GenBank/ENA/DDBJ and
subsequently picked up by RefSeq. Genome centers can do better by
augmenting these databases with additional genomes that they have
sequenced but which are not yet represented in the RefSeq collection.

1.  archaea

Query in Nucleotide :

archaea\[porgn\] AND srcdb\_refseq\[prop\] AND biomol\_genomic\[prop\]
AND complete\[prop\]

1.  bacteria

Query in Nucleotide :

bacteria\[porgn\] AND srcdb\_refseq\[prop\] AND biomol\_genomic\[prop\]
AND complete\[prop\]

1.  fungi

Query in Nucleotide :

fungi\[porgn\] AND srcdb\_refseq\[prop\] AND biomol\_genomic\[prop\] AND
(NC\_000000:NC\_999999\[pacc\] OR AC\_000000:AC\_999999\[pacc\] OR
(NT\_000001:NT\_999999999\[pacc\] AND ("chromosome 2L" OR "chromosome
2R" OR "chromosome 3L" OR "chromosome 3R")))

1.  arthropoda

Query in Nucleotide :

arthropoda\[porgn\] AND srcdb\_refseq\[prop\] AND
biomol\_genomic\[prop\] AND (NC\_000000:NC\_999999\[pacc\] OR
AC\_000000:AC\_999999\[pacc\] OR (NT\_000001:NT\_999999999\[pacc\] AND
("chromosome 2L" OR "chromosome 2R" OR "chromosome 3L" OR "chromosome
3R")))

1.  chordata

Query in Nucleotide :

chordata\[porgn\] AND srcdb\_refseq\[prop\] AND biomol\_genomic\[prop\]
AND (NC\_000000:NC\_999999\[pacc\] OR AC\_000000:AC\_999999\[pacc\] OR
(NT\_000001:NT\_999999999\[pacc\] AND ("chromosome 2L" OR "chromosome
2R" OR "chromosome 3L" OR "chromosome 3R")))

1.  other\_metazoa

Query in Nucleotide :

metazoa\[porgn\] NOT (arthropoda\[porgn\] OR chordata\[porgn\]) AND
srcdb\_refseq\[prop\] AND biomol\_genomic\[prop\] AND
(NC\_000000:NC\_999999\[pacc\] OR AC\_000000:AC\_999999\[pacc\] OR
(NT\_000001:NT\_999999999\[pacc\] AND ("chromosome 2L" OR "chromosome
2R" OR "chromosome 3L" OR "chromosome 3R")))

1.  viridiplantae

Query in Nucleotide :

viridiplantae\[porgn\] AND srcdb\_refseq\[prop\] AND
biomol\_genomic\[prop\] AND (NC\_000000:NC\_999999\[pacc\] OR
AC\_000000:AC\_999999\[pacc\] OR (NT\_000001:NT\_999999999\[pacc\] AND
("chromosome 2L" OR "chromosome 2R" OR "chromosome 3L" OR "chromosome
3R")))

1.  other\_eukaryota

Query in Nucleotide :

eukaryota\[porgn\] NOT (metazoa\[porgn\] OR fungi\[porgn\] OR
viridiplantae\[porgn\]) AND srcdb\_refseq\[prop\] AND
biomol\_genomic\[prop\] AND (NC\_000000:NC\_999999\[pacc\] OR
AC\_000000:AC\_999999\[pacc\] OR (NT\_000001:NT\_999999999\[pacc\] AND
("chromosome 2L" OR "chromosome 2R" OR "chromosome 3L" OR "chromosome
3R")))

1.  viruses\_and\_viroids

Query in Nucleotide :

(viruses\[porgn\] OR viroids\[porgn\]) AND srcdb\_refseq\[prop\] AND
biomol\_genomic\[prop\] AND (NC\_000000:NC\_999999\[pacc\] OR
AC\_000000:AC\_999999\[pacc\] OR (NT\_000001:NT\_999999999\[pacc\] AND
("chromosome 2L" OR "chromosome 2R" OR "chromosome 3L" OR "chromosome
3R")))

The FASTA sequence files resulting from these queries are formatted as
nine BLAST databases using the makeblastdb program.

#### Execution

Repeats in the input FASTA sequences are soft-masked to lowercase using
WindowMasker. Then BLAST hits over 98% identity are generated to the
databases for the 8 taxonomic groups to which the organism being
screened does not belong.

``` {.example}
        blastn -query _input_fasta_sequences_ -db _distant_organism_dbs_ -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -min_raw_gapped_score 100 -penalty -5 -perc_identity 98.0 -soft_masking true
```

### First pass calls

The following heuristic rules help to get rid of most false matches.

#### Process contaminant matches from 1

Contaminant matches from (1) are merged if they are from the same class
of sequence (VECTOR, E.coli, IS, PHG) and they overlap or are separated
by 50 bases or less.

If the total coverage of contaminant matches from (1) is &gt;75% of the
sequence length then flag the sequence as a contaminant to be excluded.

If the contaminant is classed as VECTOR, E.coli, IS:./, PERM:./ or
PHG:\* and the contaminant location is within 100 bases of the the start
or end of the sequence (or gap is the sequence is not contiguous), or
within 100 bases of another contaminant match that is at an end, flag
the contaminant span for trimming.

If the contaminant is one of the above, and the match is longer than 700
bases flag the contaminant span for trimming.

Other matches may be false alarms. Treat them as suspect spans and
reBLAST the hit span plus 10 Kbp of flanking sequence on each side
against nr, HTGS, related and unrelated chromosomes (as described
below).

#### Process contaminant matches from 2

Flag all adaptor spans for trimming.

#### Process mitochondrion matches from 3

If the total coverage of mitochondrial matches from (3) is &gt;75% of
the sequence length then flag the sequence as being mitochindrial
sequence to be excluded.

#### Process unrelated chromosome matches from 5

Ignore any matches to chromosomes from unrelated organisms that lie with
a region identified as being rRNA genes from (4) (the spans matched in 4
plus 100 bases on both sides). These are likely to be false matches.

Treat other matched spans as suspect and reBLAST the hit span plus 10
Kbp of flanking sequence on each side against nr, HTGS, related and
unrelated chromosomes

#### ReBLAST against nr, HTGS, related and unrelated chromosomes

Spans identified a contamination suspects in the first pass, plus 10 Kbp
of flanking sequence on each side (up to the end of the contig), are
BLASTed against nr, HTGS, related and unrelated chromosomes to generate
additional data for calling contaminants to be excluded or trimmed.

#### Databases

chromosome databases (a) to (i) from (5) above.

``` {.example}
        ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz

        ftp://ftp.ncbi.nlm.nih.gov/blast/db/htgs.*.tar.gz
```

#### Execution

The suspect spans are BLASTed against each of the 10 databases.

``` {.example}
        blastn -query _suspect_spans_plus_flanks_ -db _reblast_db_ -task megablast -dust yes -evalue 1E-9 -searchsp 1000000000 -perc_identity 98.0 -soft_masking true
```

#### Processing the reBLAST matches

Automatically exclude sequence contigs that meet all the following
criteria:

``` {.example}
        >60% of length covered with foreign hits, or less than 200 bp that are NOT covered

        Each contributing hits must be 100 bp or longer with identity >= 98%
```

The best match to chromosomes from unrelated organisms is longer than
the best match to chromosomes from the related organism group

Some of the other hits may be reviewed manually.

Usage
-----

usage: remower.py \[-h\] --genomefile GENOMEFILE \[--dbvec DBVEC\]
\[--dbmito DBMITO\] \[--dbcont DBCONT\] \[--dist DIST\]

remower is a script that allows to remove contamination in assembled
genome. It takes as input a genome file, contamination databases returns
on stdout the corrected genome and on stderr warnings regarding vector
sequences not removed.

optional arguments: -h, --help show this help message and exit
--genomefile GENOMEFILE, -g GENOMEFILE Genome file --dbvec DBVEC, -v
DBVEC The vecscreen database --dbmito DBMITO, -m DBMITO The organelle
database --dbcont DBCONT, -c DBCONT The contaminant database --dist
DIST, -d DIST Maximal distance for merging two intervals

python3 ./remower.py -c ~contaminationdb~ -v **vector~db~** -m ~mitodb~
--genomefile **genome.fa**

\#+END~SRC~
