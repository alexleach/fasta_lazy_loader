#!/usr/bin/env python

# Author : Alex Leach ©2013
# Contact: beamesleach /at/ gmail /dot/ com
# GNU GPL v3 license. Please get in touch if you have issues with this.


#   Copyright (C) 2013 Alex Leach
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
fasta_lazy_loader.py

Description:

A module for lazy loading fasta formatted sequences from a big
sequence file.

Classes:

    * FastaIndex

    Stores accession and byte locations of fasta sequences. Able to save
    indexed sequence locations to disk (in a pickle file) and reload saved
    files on demand.
    
    * FastaLazyLoader

    Parses fasta-formatted sequence files, whilst managing and populating
    a FastaIndex instance. Also contains a multiprocessing-friendly method
    for retrieving sequences from disk: `pipe_sequences(...)`


Example usage:

    A simple command-line application is included, which also serves to
    demonstrate usage (see the `main()` function, below).


    To test usage from the command line, you can just run:-

    > python fasta_lazy_loader.py -in BIG_FASTA_FILE.fasta

"""

## System library imports
import os, sys
try:
    import cPickle as pickle
except ImportError:
    import pickle

# BioPython imports
from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord as BioSeqRecord
from Bio.Alphabet import single_letter_alphabet

class FastaIndex:
    """
    FastaIndex

    Manages an index of fasta-file sequence co-ordinates. Indexes are
    stored in memory until the class is deleted or garbage collected,
    at which point the index is dumped to a pickle file.

    The index is maintained as a simple dictionary of the following
    key-value pairs:-

        { 'accession' : 'location' } 

    """
    
    handle      = None  # File handle to the index.
    index_loc   = None  # The index's file path.
    indexes     = {}    # The loaded in-memory indexes.
    indexed     = False # Can be set by client code, when all indexes
                        #   are loaded.
    suffix      = "pklindex"

    def __init__(self, name=None):
        """
        @param name - File name where to save the pickled index.
                      If the file already exists and is non-empty, load
                      the indexes from there into memory.
                      If no name is given, will create a temporary 
                      file, print its location to stdout and save there.
        """
        if name is None:
            from tempfile import mkstemp
            handle, name = mkstemp( suffix=self.suffix, prefix="fasta_index" )
        else:
            if isinstance( name, str ):
                if os.path.exists( name ) and os.path.getsize( name ) > 0:
                    with file(name, 'rb') as handle:
                        self.indexes = pickle.load(handle)
                        handle = None
                    self.indexed = True
                else:
                    handle = file(name, 'wb')
            else:
                raise TypeError("When given, name must be a file name " \
                                "(a string).\n"                         \
                                "Got {0}".format(type(name)) )

        self.handle = handle
        self.index_loc = name

    def __del__(self):
        handle = self.handle
        if self.handle is not None and \
            isinstance(handle, file) and not handle.closed:
                handle.close()

    def __getitem__(self, accession):
        return self.indexes[accession]

    def __iter__(self):
        # If we're reading all sequence records back from disk, may as
        # well do it in order.
        idx = self.indexes
        for loc in sorted( idx, key=idx.get ):
            yield idx[loc]

    def __len__(self):
        return len(self.indexes)

    def __setitem__(self, accession, location):
        self.indexes.update( { accession : location } )

    def save(self):
        """Save the current index, as a pickled dictionary, to self.handle."""
        handle = self.handle
        if handle is None:
            # We've been using a pre-loaded index. It hasn't broken, so no need
            # to overwrite it.
            return
        if isinstance(handle, file):
            if not handle.closed:
                print "Pickling {0} sequence indexes in {1}" \
                      .format( len(self.indexes), self.index_loc )
                pickle.dump(self.indexes , handle, -1)
            else:
                raise IOError("{0} is already closed!".format(handle.name))
        else:
            raise IOError("{0} is not a file object".format(type(handle)))


class FastaLazyLoader( object ):
    """FastaLazyLoader

    A fasta parser designed to work with files too large for your 
    system's memory.

    Keep the file object open, and move around using the ArbIO.index
    This contains byte locations of each record, but if indexed, then
    the sequence will be returned.
    
    @param in_handle  - readable file-like object containing Fasta formatted
                        sequences.
                        Alternatively, pass a file name as a string, which 
                        will be opened automatically.                   [stdin]
    @param out        - writable file-like object.                     [stdout]
    @param index      - existing index. By default, create a new one.    [None]
    
    """
    alphabet    = single_letter_alphabet
    in_handle   = None
    indexes     = None
    index_file  = None
    info        = None
    out_handle  = None
    suffix      = 'pklindex'

    def __init__(self, in_handle=None, out=None, index=None, alphabet=None):
        """Initialise FastaLazyLoader.

        If `in_handle` is given and is a file name, it is opened, and the
        handle is saved to `self.in_handle`. If not given, defaults
        to reading from stdin.

        If `out_handle` is given, a new file is opened in read-write mode,
        unless a file already exists in that location, in which case the 
        user is asked for confirmation.
        """

        ## Deal with the in_handle option.
        if in_handle is None:
            # Default input is read from stdin
            self.in_handle = sys.stdin
        else:
            # Make sure `in_handle` is a file-like object.
            if hasattr(in_handle, "read"):
                self.in_handle = in_handle
            elif isinstance(in_handle, str) and os.path.exists(in_handle):
                self.in_handle = file(in_handle, 'r')
            else:
                raise IOError(
                        "{0} not a handle or valid filename".format(in_handle))

        ## Deal with the given `out` option.
        prefix = None
        if out is None:
            # Default output is written to stdout
            self.out_handle = sys.stdout
        else:
            # Make sure `out` is a file-like object.
            # Get prefix - the file name for the index
            if isinstance(out, file):
                self.out_handle = out
                prefix = out.name.rsplit('.', 1)[0]
            elif hasattr(out, "write") and hasattr(out, "seek"):
                self.out_handle = out
                prefix = None
            elif isinstance(out, str):
                if not os.path.exists(out):
                    print "Saving output sequences to '{0}'".format(out)
                    self.out_handle = file(out, 'w+r')
                    prefix = out.rsplit('.', 1)[0]
                else:
                    # Get a non-existing file name, or confirm to overwrite
                    while os.path.exists(out):
                        reply = raw_input("{0} exists.  "               \
                            "Are you sure you want to overwrite it?\n"  \
                            " [Enter 'y' or a new filename] :\n ".format(out))
                        if reply == 'y':
                            self.out_handle = file(out, 'w+r')
                            break
                        out = reply
                    prefix = out.rsplit('.', 1)[0]
            else:
                # Don't save any output
                prefix = self.in_handle.name.rsplit('.', 1)[0]
                self.out_handle = None

        # Create a FastaIndex instance
        if index is None:
            self.indexes = FastaIndex()
        else:
            if prefix is None:
                index_name = None
            else:
                index_name = '.'.join([prefix, self.suffix])
            self.indexes = FastaIndex(index_name)

        if alphabet is not None:
            self.alphabet = alphabet

    def __del__(self):
        self.close()

    def __getitem__(self, accession):
        """The main lazy-loading magic! Given an accession or ID number,
        parsed from a fasta-formatted sequence file, look it up in the
        internal index, go to it in the original file, and return it as
        a BioPython SeqRecord."""
        location = self.indexes[accession]
        handle = self.in_handle
        handle.seek(location)
        line = handle.readline()
        if not line.startswith('>'):
            print line
            msg = "Error at indexed file position {0} : {1}\n" \
                  "Index looks corrupt. Delete it and re-run"
            raise ValueError(msg)

        id, desc = self.parse_header( line )
        if id != accession:
            msg = "Sequence index looks corrupt! :(\n"  \
                  "\tExpected: {0}\n"                   \
                  "\tGot: {1}\n"                        \
                  "\tWhole line:-\n"                    \
                  "{2}".format( accession ,id, line )
            raise IndexError(msg)
        sequence = ''
        while True:
            line = handle.readline()
            if line[:1] == '>' or len(line) == 0:
                break
            else:
                sequence += line.rstrip()
        return BioSeqRecord(
            BioSeq(sequence, alphabet=self.alphabet),
            id=id, description=desc )

    def __iter__(self):
        for seq in self.indexes:
            yield seq

    def __len__(self):
        return len(self.indexes)

    def __setitem__(self, accession, location):
        self.indexes[accession] = location

    def close(self):
        """Closes 3 open handles:-
        self.in_handle,
        self.out_handle,
        self.index_file (if it's open, save contents of self.indexes there)"""
        self.in_handle.close()
        if isinstance(self.out_handle, file) and not self.out_handle.closed:
            self.out_handle.close()
        if self.indexes is not None:
            self.indexes.save()
            del(self.indexes)

    def fetch(self, accessions):
        """Give a list of accessions and this will yield them as SeqRecord objects
        by reading them from the rewritten file."""
        for accession in accessions:
            yield self[accession]

    def index(self):
        """ Indexes self.in_handle for sequence locations.
        Returns the internal FastaIndex instance to the caller.
        """
        handle = self.in_handle
        handle.seek(0)
        while True:
            prev_loc = handle.tell()
            line = handle.readline()
            if line.startswith('>'):
                accession, desc = self.parse_header(line)
                self[accession] = prev_loc
            elif len(line) == 0:
                break
            else:
                continue
        return self.indexes

    def index_and_info(self):
        """Indexes self.in_handle for sequence accessions and locations, whilst
        simultaneously storing a separate dictionary of FASTA descriptions.

        The dictionary of FASTA descriptions is saved into self.info.
        The number of sequences indexed is returned.
        """
        handle = self.in_handle
        handle.seek(0)
        info = self.info if self.info is not None else {}
        count = 0
        while True:
            prev_loc = handle.tell()
            line = handle.readline()
            if line[:1] == '>':
                accession, desc = self.parse_header(line[1:])
                self[accession] = prev_loc
                info.update( { accession : { 'desc' : desc } } )
                count += 1
            elif len(line) == 0:
                break
            else:
                continue
        self.info = info
        return count

    def index_and_save( self, out_handle, format='fasta'):
        """Writes input sequences to out_handle in the specified format.
        Updates self.index simultaneously, to match the new file co-ordinates.
        """
        try:
            pos = out_handle.tell()
        except IOError, e:
            # not seekable. out_handle has no context of position.
            msg = "{0}\n    {1} has no knowledge of position."\
                            .format(e, out_handle)
            raise IOError(msg)

        from Bio.SeqIO import write
        idx = self.indexes
        for seq_record in self.parse(self.in_handle):
            write( seq_record, out_handle, format )
            idx[seq_record.id] = pos
            pos = out_handle.tell()
        idx.save()


    @classmethod
    def parse(cls, handle):
        """An iterator function that yield's SeqRecord objects from a readable
        file-like object full of fasta-formatted sequences.

        This serves the same purpose as BioPython's SeqIO.FastaIO.FastaIterator 
        function. It doesn't actually do any indexing, as is intended to be
        used when indexing an output file, as in `self.index_and_save(...)`
        """
        sequence = ''
        alphabet = cls.alphabet
        line = handle.readline()
        while not line.startswith('>'):
            line = handle.readline()
            if not line:
                raise IOError("Found no fasta sequence in {0}"\
                              .format(handle.name))
        id, desc = cls.parse_header(line)

        for line in handle:
            if line.startswith('>'):
                yield BioSeqRecord( 
                        BioSeq(sequence, alphabet=alphabet),
                        id=id, description=desc )

                id, desc = cls.parse_header(line)
                sequence = ''
            else:
                sequence += line.rstrip()

        yield BioSeqRecord(
                BioSeq(sequence, alphabet=alphabet),
                id=id, description=desc )

    @classmethod
    def parse_header(cls, seq_header):
        """Given a fasta sequence header, will return the tuple:-
            (id, description).
        id is everything up until the first bit of whitespace.
        description is the rest of the header string."""
        if seq_header.startswith('>'):
            header = seq_header[1:].strip().split(' \t', 1)
        else:
            header = seq_header.strip().split(' \t', 1)

        desc = ''
        if len(header) == 1:
            id = header[0]
        elif len(header) == 2:
            id = header[0]
            desc = header[1]
        else:
            id = header[0]

        return (id, desc)


    def pipe_sequences(self, accessions, out_pipe):
        """For applications using multiprocessing, specific sequences can be
        retrieved by passing an iterable sequence of accessions and a
        multiprocessing.Pipe instance. For each accession found in the index,
        it will be sent down the pipe.
        If an accession is not found in the index, an error is raised and the
        in the current thread, and the exception sent down the pipe.
        """
        try:
            for accession in accessions:
                if len( accession ) > 0:
                    seq = self[accession] #.format('fasta')
                    out_pipe.send( seq )
        except IndexError, e:
            out_pipe.send( e )
            raise(e)
        finally:
            out_pipe.close()

def main():
    # Get some options from the command line
    if len ( sys.argv ) == 1:
        sys.argv.append('-h')
    from argparse import ArgumentParser
    parser = ArgumentParser("Index a fasta file for lazy-loading")
    parser.add_argument("-in",    type=str,  default=None, dest='in_' )
    parser.add_argument("-out",   type=str,  default=None )
    parser.add_argument("-index", type=bool, default=True )
    options = parser.parse_args( sys.argv[1:] )

    # Initialise file handles
    lazy_loader = FastaLazyLoader( in_handle=options.in_,
                                   out=options.out,
                                   index=options.index )

    # 3 ways to index the input file, doing everything we want in a single
    #  pass:-
    #

    # 1. We can index the sequence file, getting in return the FastaIndex.
    print "indexing file with lazy_loader.index()"
    fasta_index = lazy_loader.index()

    # (For testing and demonstrative purposes, seek the file back to the start
    # and arbitrarily index again.)
    lazy_loader.in_handle.seek(0)

    if options.out is not None:
        # 2. Rewrite sequences to file named by `out`, and index that file 
        #    rather than the input file.
        print "indexing file with lazy_loader.index_and_save(\"{0}\")" \
              .format(options.out)
        lazy_loader.index_and_save(lazy_loader.out_handle)

    else:
        # 3. Index the input file, storing descriptive information in
        #    `lazy_loader.info` and location information in 
        #    `lazy_loader.indexes`
        print "indexing file with lazy_loader.index_and_info()"
        n_seqs = lazy_loader.index_and_info()

    print "Indexed {0} sequences.".format(len(lazy_loader))

    print "iterating through lazy-loader"
    # test iteration
    i = 0
    for seq in lazy_loader:
        i += 1
    print "Iterated through {0} sequences.".format(i)
    print "Completed everything successfully."

if __name__ == '__main__':
    main()
