# rnaPUFdomains
Perl script for designing pumilio domains to bind specific RNA targets of interest

**Note**: this is a script I generated a number of years ago, in a spare moment and is provided without any guarantees the designed fusions will be functional!

For the target ACGTACGT the command to use this is:
```perl
  perl GGPUF.pl ACGTACGT
```

The script is based around the sequences from the PUF assembly kit availible from Addgene ([link](https://www.addgene.org/kits/puf-assembly-kit/). This was described in the paper:

Modular assembly of designer PUF proteins for specific post-transcriptional regulation of endogenous RNA. Abil Z, Denard CA, Zhao H. J Biol Eng . 2014 Mar 1;8(1):7. doi: 10.1186/1754-1611-8-7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24581042) .


It provides 5 outputs:
1. It returns your input sequence
2. It provides the PUF assembly kit plasmid identities required to construct the designed binding sequence by golden gate cloning using the PUF assembly kit.
3. It provides the full length DNA sequence for the designed PUF domain to permit synthesis as an e.g. gBlock.
4. It provides the DNA sequence as a GFP10 fusion to permit split GFP imaging
5. It provides the DNA sequence as a GFP11 fusion to permit split GFP imaging.

Ed Emmott, 2019, Northeastern U.
email: e.emmott@northeastern.edu
WWW: http://edemmott.co.uk
twitter: http://twitter.com/edemmott
