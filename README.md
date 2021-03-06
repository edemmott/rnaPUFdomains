# rnaPUFdomains
Perl script for designing PUF domains to bind specific RNA targets of interest. This is based around the sequences from the PUF assembly kit availible from Addgene ([link](https://www.addgene.org/kits/puf-assembly-kit/)). This was described in the paper:

Modular assembly of designer PUF proteins for specific post-transcriptional regulation of endogenous RNA. Abil Z, Denard CA, Zhao H. J Biol Eng . 2014 Mar 1;8(1):7. doi: 10.1186/1754-1611-8-7. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/24581042).

It identifies the required kit components for a desired fusion, as well as producing DNA sequences to permit direct synthesis of PUF domains of interest, either alone or as GFP10/11 fusions for split GFP reconstitution.

This script has been successfully used on Max OSX, with Perl v5.26.2.

**Note: this is a script I generated a number of years ago, in a spare moment and is provided without any guarantees the designed fusions will be functional!**



## Use:
For the target ACGTACGT the command to use this is:
```perl
  perl GGPUF.pl ACGTACGT
```

The script will convert to upper case, as well as substitute T for U.

It provides 5 outputs:
1. It returns your input sequence
2. It provides the PUF assembly kit plasmid identities required to construct the designed binding sequence by golden gate cloning using the PUF assembly kit.
3. It provides the full length DNA sequence for the designed PUF domain to permit synthesis as an e.g. gBlock.
4. It provides the DNA sequence as a GFP10 fusion to permit split GFP imaging
5. It provides the DNA sequence as a GFP11 fusion to permit split GFP imaging.

## Example output:
Output for the above target ACGTACGT from above is:
```
Your input sequence was: ACGUACGU
Base 6 is: C
The PUF assembly modules you require are for this sequence are:
8CQ 7SYR 6SE 5NQ 4CYQ 3SR 2SE 1NQ
The DNA sequence for your PUF assembly is: 
ATGCGCAGCCGCCTTTTGGAAGATTTTCGAAACAACCGGTACCCCAATTTACAACTGCGGGAGATTGCCGGACATATAATGGAATTTTCCCAAGACCAGCATGGGAACAGATTCATTCAGCTGAAACTGGA
GCGTGCCACACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAGGCTGCCTACCAACTCATGGTGGATGTGTTTGGTAGTTACGTCATTGAGAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGTCACGTCCTGTCATTGGCACTACAGATGTATGGCTCCCGTGTTATCCGCAAAG
CTCTTGAGTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGCCATGTCTTGAAGTGTGTGAAAGATCAGAATGGCTGTTACGTGGTTCAGAAA
TGCATTGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGGTATTTGCCTTATCCACACATCCTTATGGCAACCGAGTGATTCAGAGAATCCTGGAGCAC
TGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACACAGAGCAGCTTGTACAGGATCAATATGGAAGTTATGTAATCGAACATGTACTGGAGCAC
GGTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATGTACTTGTATTGAGTCAGCACAAATTTGCAAGCTATGTTGTGCGCAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
TATACACCATGATGAAGGACCAGTATGCCTGCTACGTGGTCCAGAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCCCACATCGCAACTCTTCGTAAGTACACCTATGGCAAGCACATTCTGGCCAAGCTGGAGAAGTACTACATGAAGAACGGTGTTGACTTAGGGTAA
The DNA sequence for your GFP10 PUF assembly is: 
ATGGACCTGCCCGACGACCACTACCTGAGCACCCAGACCATCCTGAGCAAGGACCTGAACGACGTGGGCAGCGGCGGCGGCAGCCACATGCGCAGCCGCCTTTTGGAAGATTTTCGAAACAACCGGTACCCCAATTTACAACTGCGGGAGATTGCCGGACATATAATGGAATTTTCCCAAGACCAGCATGGGAACAGATTCATTCAGCTGAAACTGGA
GCGTGCCACACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAGGCTGCCTACCAACTCATGGTGGATGTGTTTGGTAGTTACGTCATTGAGAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGTCACGTCCTGTCATTGGCACTACAGATGTATGGCTCCCGTGTTATCCGCAAAG
CTCTTGAGTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGCCATGTCTTGAAGTGTGTGAAAGATCAGAATGGCTGTTACGTGGTTCAGAAA
TGCATTGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGGTATTTGCCTTATCCACACATCCTTATGGCAACCGAGTGATTCAGAGAATCCTGGAGCAC
TGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACACAGAGCAGCTTGTACAGGATCAATATGGAAGTTATGTAATCGAACATGTACTGGAGCAC
GGTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATGTACTTGTATTGAGTCAGCACAAATTTGCAAGCTATGTTGTGCGCAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
TATACACCATGATGAAGGACCAGTATGCCTGCTACGTGGTCCAGAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCCCACATCGCAACTCTTCGTAAGTACACCTATGGCAAGCACATTCTGGCCAAGCTGGAGAAGTACTACATGAAGAACGGTGTTGACTTAGGGTAA
The DNA sequence for your GFP11 PUF assembly is: 
ATGCGCAGCCGCCTTTTGGAAGATTTTCGAAACAACCGGTACCCCAATTTACAACTGCGGGAGATTGCCGGACATATAATGGAATTTTCCCAAGACCAGCATGGGAACAGATTCATTCAGCTGAAACTGGA
GCGTGCCACACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAGGCTGCCTACCAACTCATGGTGGATGTGTTTGGTAGTTACGTCATTGAGAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGTCACGTCCTGTCATTGGCACTACAGATGTATGGCTCCCGTGTTATCCGCAAAG
CTCTTGAGTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGCCATGTCTTGAAGTGTGTGAAAGATCAGAATGGCTGTTACGTGGTTCAGAAA
TGCATTGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGGTATTTGCCTTATCCACACATCCTTATGGCAACCGAGTGATTCAGAGAATCCTGGAGCAC
TGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACACAGAGCAGCTTGTACAGGATCAATATGGAAGTTATGTAATCGAACATGTACTGGAGCAC
GGTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATGTACTTGTATTGAGTCAGCACAAATTTGCAAGCTATGTTGTGCGCAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
TATACACCATGATGAAGGACCAGTATGCCTGCTACGTGGTCCAGAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCCCACATCGCAACTCTTCGTAAGTACACCTATGGCAAGCACATTCTGGCCAAGCTGGAGAAGTACTACATGAAGAACGGTGTTGACTTAGGGGGCAGCGGCGGCGGCAGCGGCGGCGGCAGCGAGAAGAGAGACCACATGGTGCTGCTGGAGTACGTGACCGCCGCCGGCATCACCGACGCCAGCTAA
```


Ed Emmott, 2019, Northeastern U.

email: e.emmott@northeastern.edu
WWW: http://edemmott.co.uk
twitter: http://twitter.com/edemmott
