# blast2cap3

`blast2cap3` merges contigs together from transcriptome assemblies
using protein homology and CAP3.

## Motivation

RNA-seq has given rise to a powerful way of exploring the genes of
non-model organisms lacking a full genome. Howver, unlike genomic DNA
sequencing, RNA-seq coverage is dependent upon gene expression levels,
which can confound the assembly process. There are
transcriptome-specific assmeblers (such as
[Trinity](http://trinityrnaseq.sourceforge.net/) and
[Oases](http://www.ebi.ac.uk/~zerbino/oases/)) that address this;
however, the best transcriptome assemblies are often those that
incorporate many assemblies with varying sized k-mers. C. Titus Brown
has [a nice
introduction](http://ivory.idyll.org/blog/the-k-parameter.html) to the
importance of *k* in transcriptome assemblies.

After multiple assemblies with different k-mer sizes, the resulting
contigs should be merged. The first step is to cluster nucleotide
sequences via a program like
[CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/), so contigs that are
fully contained in another contig are merged. The second step would be
to join contigs that may have partially assembled across differently
sized k-mers. For example:


    |----------- contig k31-a ---------------|
                                       |----------- contig k61-b -----------|
    

Assemblers like CAP3 are then often used to merge contigs based on the
overlapping region. However, this strategy often either (1) overwhelms
CAP3 with too many contigs to consider merging or (2) incorrectly
merges two contigs with nucleotide similarity.

`blast2cap3` alleviates these issues by narrowing the pool of possible
contigs to merge to those that share a BLASTX hit to a common subject
protein of a relative. `blast2cap3` uses protein similarity because
most transcriptome contigs will be protein coding, and leveraging this
increases specifity. Both the narrowing of the merging pool and the
use of similarity in protein space rather than nucleotide space leads
to post-assembly merges with more specificity (more testing is needed,
but we're using it in production).

`blast2cap3` is best used after ORF prediction. Contigs with full ORFs
shouldn't be candidates for a merge (unless you're interested in
assembly the 5' and 3' UTRs; if this is the case, use a different
approach). An iterative approach of ORF prediction, `blast2cap3`, ORF
prediction, `blast2cap3`, etc may be a good strategy.

## Citation

If you use this software, please cite:

Vince Buffalo, Ksenia Krasileva, and Jorge Dubcovsky, "blast2cap3:
transcriptome contig merging via protein homology", 2012.