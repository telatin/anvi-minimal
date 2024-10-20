# anvi-minimal

Minimal example of metagenome:

* ref/contigs.fa
* reads/

The reference:

```text
┌────────────┬──────┬──────────┬───────┬───────┬───────┐
│ File       │ #Seq │ Total bp │ N50   │ Min   │ Max   │
├────────────┼──────┼──────────┼───────┼───────┼───────┤
│ contigs.fa │ 3    │ 10,285   │ 2,833 │ 2,341 │ 5,111 │
└────────────┴──────┴──────────┴───────┴───────┴───────┘
```

## Generate Anvi'o databases

1. Activate Anvio-8 environment
2. You will need `bwa` and `samtools`
3. Run  `snakemake --cores 8` (`make` is available, but only produce the databases)
   1. This wil map the reads and make the BAM files
   2. Create CONGIGS.db
   3. Create and merge the coverage profiles

## Adding metadata from text files

It would be ideal to add tracks at the:

* contigs (or split) level
* gene level

In `extras/metadata.txt` there is a possible table at the split level to feed with `-d` to `anvi-interactive`.

```text
split                           Source          Phylum
c_000000000001A_split_00001     Reference       Nucleocytoviricota
c_000000000001B_split_00001     Inference       Nucleocytoviricota
c_000000000002A_split_00001     Reference       Aliumviricota
```

:question: How to add metadata at at the split level? How should the header look like?


:question: How to add metadata at at the gene level? 
