# fasta_lazy_loader.py
======================

## Description:

A module for lazy loading fasta formatted sequences from a big
sequence file.

## Classes:

  * FastaIndex

    Stores accession and byte locations of fasta sequences. Able to save
    indexed sequence locations to disk (in a pickle file) and reload saved
    files on demand.
    
  * FastaLazyLoader

    Parses fasta-formatted sequence files, whilst managing and populating
    a FastaIndex instance. Also contains a multiprocessing-friendly method
    for retrieving sequences from disk: `pipe_sequences(...)`


## Example usage:

A simple command-line application is included, which also serves to
demonstrate usage (see the `main()` function, below).


To test usage from the command line, you can just run:-

  > python fasta_lazy_loader.py -in BIG_FASTA_FILE.fasta
