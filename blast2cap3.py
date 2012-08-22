info = """
blast2cap3 joins transcriptome contigs by finding links between
contigs via homologous proteins. It then uses CAP3 with very high
sequence similarity to merge these candidate contigs by overlap.
"""

import os
import pdb
import sys
import shutil
import re
from collections import defaultdict
import subprocess
import tempfile
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def go_interactive(**kwargs):
    try:
        import readline
        import code
    except ImportError:
        pass
    else:
        import rlcompleter
        readline.parse_and_bind("tab: complete")                    
        
    code.InteractiveConsole(locals=kwargs).interact()

def parse_CAP3_out(filename, verbose=True):
    """
    Parse CAP3's out file, getting the joined contigs.

    """
    header_matcher = re.compile("\*+ Contig (\d+) \*+")
    joined = defaultdict(list)

    try:
        with open(filename) as cap3_file:
            line = next(cap3_file)

            # eat up header until Overlaps line
            while not line.strip().startswith("Overlaps"):
                line = next(cap3_file)

            # we're in the header section that indicates joined
            # contigs
            next(cap3_file) # pop off blank line
            line = next(cap3_file)

            ## now we're in the section we need to be in: headers of
            ## joined Contig IDs and the contigs that went into them.
            while True:
                header = header_matcher.search(line)
                if header is None:
                    # no output
                    return joined
                contig_num = header.group(1)

                line = next(cap3_file)
                while not line.startswith("*******************"):
                    contig_id = "Contig" + contig_num
                    joined[contig_id].append(line.strip()[:-1])

                    # contig entries end in + or -; if not, let's
                    # error our
                    if line.strip()[-1:] not in ("+", "-"):
                        raise ValueError("error: unexpected input while parsing CAP3. "
                        "Were sequences clustered first using CD-HIT?")
                    if verbose:
                        sys.stderr.write("[parse_CAP3_out] adding '%s' to '%s'\n" %
                                         (line.strip(), contig_id))
                    line = next(cap3_file)
                    if line.strip() == "" or line.startswith("DETAILED DISPLAY"):
                        return joined
    except StopIteration:
        # empty file, most likely we may want to warn
        return joined

def run_CAP3(sequences, subject_id, percent_identity=99, clipping=False,
             min_overlap=100, clean=True, verbose=True, debug=True):
    """
    Run CAP3, which creates a ton of files, which we process and
    remove afterwards.
    """
    tempdir = tempfile.mkdtemp("temporary_cap3")

    # write sequences to temporary directory
    contigs_file = os.path.join(tempdir, "contigs.fasta")
    with open(contigs_file, 'w') as f:
        SeqIO.write(sequences.values(), f, "fasta")

    # run CAP3 in temporary directory
    cmd = "cap3 %s -p %s -k %s -o %s > cap3.out"
    cmd_with_values = cmd % (os.path.basename(contigs_file),
                             percent_identity, 
                             int(clipping), min_overlap)
    full_cmd = "cd %s && " % tempdir + cmd_with_values
    if verbose:
        sys.stderr.write("[run_CAP3] executing command: %s\n" % full_cmd)
    status = subprocess.call(full_cmd, shell=True)

    if status != 0:
        raise ValueError("error: CAP3 returned non-zero exit status")

    # joined is a dictionary, with the keys being the CAP3 contig ids,
    # and the values are a list of the original contig ids it joined
    joined = parse_CAP3_out(os.path.join(tempdir, "cap3.out"))

    # joined_contigs_file is the actual fasta file of CAP3 joined
    # contigs, with the headers that correspond to joined
    joined_contigs_file = os.path.join(tempdir, "contigs.fasta.cap.contigs")

    # We want to read the contig file produced by CAP3 into memory,
    # changing the headers to something that is informative and that
    # can be parsed downstream if necessary. 
    joined_contigs_seqs = dict()
    with open(joined_contigs_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            new_id = ("joined-contig|subject-protein-link:%s|components:%s|cap3:%s" %
                      (subject_id, ','.join(joined[record.id]), record.id))
            joined_contigs_seqs[new_id] = record.seq

            # write stats to standout
            msg = "\t".join([subject_id, record.id]) + "\n"
            if debug:
                msg += "cap3 dir: %s\n" % tempdir
                
            for orig_contig_id in joined[record.id]:
                orig_contig_len = str(len(sequences[orig_contig_id]))
                msg += "\t".join(["", orig_contig_id, orig_contig_len])
                msg += "\n"
            msg += "joined length: %s\n" % len(record.seq)
            sys.stdout.write(msg)

    if clean and not debug:
        if verbose:
            sys.stderr.write("[run_CAP3] removing temp directory %s\n" % tempdir)
        shutil.rmtree(tempdir)
    return (joined_contigs_seqs, joined)


def get_contig_links(file, ignore_contigs=set()):
    """
    This builds a dictionary of a list of query (contig) IDs per each
    subject they map to, linking contigs by subject BLASTX hit
    protein.
    """
    d = defaultdict(list)
    ignore_contigs = set(ignore_contigs) # faster lookup
    for line in file:
        line = line.strip()
        try:
            query, subject = line.split("\t")[:2]
        except ValueError:
            raise ValueError("error parsing tabular BLASTX output: format not tabular")
        if query not in ignore_contigs:
            d[subject].append(query)
    file.close()
    return d

def load_exclude_file(file):
    """
    Load the exclude file, raising error if it is formatted
    incorrectly.
    """
    matcher = re.compile("^[\d\w-]+$")
    exclude = set()
    if file is None:
        return exclude
    for line in file:
        if matcher.match(line) is None:
            raise ValueError("improperly formatted exclude file - "
                             "one identifier per line only")
        exclude.add(line.strip())
    file.close()
    return exclude

def main(args):
    """
    Run blast2cap3 with the argument input.
    """
    exclude_ids = load_exclude_file(args.exclude)
    blastx_joined_contigs = get_contig_links(args.blast)
    contigs = SeqIO.to_dict(SeqIO.parse(args.contigs, "fasta"))

    cap3_joined_contigs = dict() # for writing joined fasta file
    all_joined_contigs = set() # for subsetting original contigs

    for subject_link in blastx_joined_contigs:
        if len(blastx_joined_contigs[subject_link]) == 1:
            # only one BLASTX subject hit; ignore
            continue

        # Make a dictionary of ids and sequences
        seqs = dict([(k, contigs[k]) for k in blastx_joined_contigs[subject_link]])

        # run_CAP3 on these sequences; pass in subject protein key for
        # FASTA header creation.
        join_contigs, join_info = run_CAP3(seqs, subject_link,
                                           debug=args.debug,
                                           verbose=args.verbose)
        cap3_joined_contigs[subject_link] = join_contigs

        # make a set of all contigs that were joined, as we want to
        # remove these from the original contigs
        for v in join_info.values():
            all_joined_contigs.update(v)

    if args.verbose:
        sys.stderr.write("[blast2cap3] writing unjoined contigs\n")
    unjoined_contigs = (contigs[k] for k in contigs if k not in all_joined_contigs)
    SeqIO.write(unjoined_contigs, args.unjoined, "fasta")

    if args.verbose:
        sys.stderr.write("[blast2cap3] writing joined contigs\n")
    joined_contigs = list()
    for joined in cap3_joined_contigs.values():
        # cap3_joined_contigs is a dictionary with the subject
        # protein as key
        for seq_id, seq in joined.items():
            joined_contigs.append(SeqRecord(seq, seq_id, ''))
    SeqIO.write(joined_contigs, args.joined, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-b', '--blast', help="tabular BLASTX output.",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-c', '--contigs', help="FASTA file of contig sequences.",
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-e', '--exclude', help=("plaintext newline-separated list"
                        "of contigs to ignore (i.e. those that have full length ORFs)"),
                        required=False, default=None, type=argparse.FileType('r'))
    parser.add_argument('-v', '--verbose', help="output status verbosely",
                        default=False, action="store_true")
    parser.add_argument('-d', '--debug', help="don't delete CAP3 output",
                        default=False, action="store_true")
    parser.add_argument('-j', '--joined', help="the filename to write joined contigs to",
                        default="joined.fasta", type=argparse.FileType('w'))
    parser.add_argument('-u', '--unjoined', help="the filename to write unjoined contigs to",
                        default="unjoined.fasta", type=argparse.FileType('w'))
    
    args = parser.parse_args()

    main(args)
