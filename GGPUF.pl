#!/usr/bin/perl -w
#Generating PUF donains against an 8nt sequence of interest
#
#Suitable inputs: 8 nucleotide RNA sequence. ACGU, T will be interpreted as U.
#
#
#Outputs include:
#   1) Plasmid names for construction of a designer PUF domain using the Addgene
#       PUF assembly kit.
#
#           Addgene webpage: https://www.addgene.org/kits/puf-assembly-kit/
#           Rference: Modular assembly of designer PUF proteins for specific post-transcriptional regulation of endogenous RNA. Abil Z, Denard CA, Zhao H. J.  Biol. Eng. 2014 Mar 1;8(1):7. doi: 10.1186/1754-1611-8-7
#
#   2) DNA sequence of the assembled PUF DNA for gene syntehesis/gBlock  
#       construction
#
#
# Feb 2016, Ed Emmott (University of Cambridge)
#
# Section 1: Input and sanity checking
#******************************************************************************

my $inputsequence = $ARGV[0];
$inputsequence = uc $inputsequence;
$inputsequence =~ s/T/U/g; #swaps T for U.


if (length $inputsequence != 8){
    die "Input must only contain the characters ACGU/T and be 8 bases in length.\n";
} 

if ($inputsequence =~ tr/ACGTU//c) {
       die "Input must only contain the characters ACGU/T and be 8 bases in length.\n"; 
}

#Important sanity checks include length (8nt), composition ACGU
#If fails sanity check then print error message to screen and quit.

my $base1 = substr $inputsequence,0,1;
my $base2 = substr $inputsequence,1,1;
my $base3 = substr $inputsequence,2,1;
my $base4 = substr $inputsequence,3,1;
my $base5 = substr $inputsequence,4,1;
my $base6 = substr $inputsequence,5,1;
my $base7 = substr $inputsequence,6,1;
my $base8 = substr $inputsequence,7,1;




# Section 2: Assigning variables
#******************************************************************************
#Nucleotides 1-8 PUF GG variables
my %nt1 = ('A' => '8CQ',
            'C' => '8SR',
            'G' => '8SE',
            'U' => '8NQ',
            );
my %nt2 = ('A' => '7CQ',
            'C' => '7SYR',
            'G' => '7SE',
            'U' => '7NQ',
            );
my %nt3 = ('A' => '6CQ',
            'C' => '6SR',
            'G' => '6SE',
            'U' => '6NQ',
            );
my %nt4 = ('A' => '5CQ',
            'C' => '5SR',
            'G' => '5SE',
            'U' => '5NQ',
            );
my %nt5 = ('A' => '4CQ',
            'C' => '4SR',
            'G' => '4SE',
            'U' => '4NQ',
            );
#If nucleotide 6 is a 'C', an alternative set of modules for nucleotide 5 MUST be used.
my %nt5alt = ('A' => '4CYQ',
                'C' => '4SYR',
                'G' => '4SYE',
                'U' => '4NYQ',
                );
my %nt6 = ('A' => '3CQ',
            'C' => '3SR', # See note on required use of nt5alt set.
            'G' => '3SE',
            'U' => '3NQ',
            );
my %nt7 = ('A' => '2CQ',
            'C' => '2SR',
            'G' => '2SE',
            'U' => '2NQ',
            );
my %nt8 = ('A' => '1SQ',
            'C' => '1SR',
            'G' => '1SE',
            'U' => '1NQ',
);
            
#Nucleotides 1-8 DNA sequence for gene synthesis/gBlocks. Note As BsaI digestion #generates a 4bp overlap at the 5' and 3' of each joining segment, the 5' overlap #has been removed and 3' retained to allow segment joining.

my %nt1x = ('A' => 'TATACACCATGATGAAGGACCAGTATGCCTGCTACGTGGTCCAGAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCC', #8CQ
            'C' => 'TATACACCATGATGAAGGACCAGTATGCCAGCTACGTGGTCCGCAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCC', #8SR
            'G' => 'TATACACCATGATGAAGGACCAGTATGCCAGCTACGTGGTCGAGAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCC', #8SE
            'U' => 'TATACACCATGATGAAGGACCAGTATGCCAACTACGTGGTCCAGAAGATGATTGACGTG
GCGGAGCCAGGCCAGCGGAAGATCGTCATGCATAAGATCCGACCC', #8NQ
            );
my %nt2x = ('A' => 'TACTTGTATTGAGTCAGCACAAATTTGCATGCAATGTTGTGCAGAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
', #7CQ
            'C' => 'TACTTGTATTGAGTCAGCACAAATTTGCAAGCTATGTTGTGCGCAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
', #7SYR
            'G' => 'TACTTGTATTGAGTCAGCACAAATTTGCAAGCAATGTTGTGGAGAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
', #7SE
            'U' => 'TACTTGTATTGAGTCAGCACAAATTTGCAAACAATGTTGTGCAGAAGTGTGTTACTCAC
GCCTCACGTACGGAGCGCGCTGTGCTCATCGATGAGGTGTGCACCATGAACGACGGTCCCCACAGTGCCT
', #7NQ
            );
my %nt3x = ('A' => 'CAGAGCAGCTTGTACAGGATCAATATGGATGTTATGTAATCCAACATGTACTGGAGCACGGTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATG', #6CQ
            'C' => 'CAGAGCAGCTTGTACAGGATCAATATGGAAGTTATGTAATCCGCCATGTACTGGAGCACG
GTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATG', #6SR
            'G' => 'CAGAGCAGCTTGTACAGGATCAATATGGAAGTTATGTAATCGAACATGTACTGGAGCAC
GGTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATG', #6SE
            'U' => 'CAGAGCAGCTTGTACAGGATCAATATGGAAATTATGTAATCCAACATGTACTGGAGCAC
GGTCGTCCTGAGGATAAAAGCAAAATTGTAGCAGAAATCCGAGGCAATG', #6NQ
            );
my %nt4x = ('A' => 'TATTTGCCTTATCCACACATCCTTATGGCTGCCGAGTGATTCAGAGAATCCTG
GAGCACTGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACA', #5CQ
            'C' => 'TATTTGCCTTATCCACACATCCTTATGGCTCCCGAGTGATTCGCAGAATCCTGGAGCAC
TGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACA', #5SR
            'G' => 
'TATTTGCCTTATCCACACATCCTTATGGCTCCCGAGTGATTGAGAGAATCCTGGAGCAC
TGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACA', #5SE
            'U' => 'TATTTGCCTTATCCACACATCCTTATGGCAACCGAGTGATTCAGAGAATCCTGGAGCAC
TGTCTCCCTGACCAGACACTCCCTATTTTAGAGGAGCTTCACCAGCACA', #5NQ
            );
my %nt5x = ('A' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCTGTCACGTGGTTCAGAAATGCATTGA
ATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4CQ
            'C' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCAGTCACGTGGTTCGCAAAT
GCATTGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4SR
            'G' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCAGTCACGTGGTTGAGAAATGCATTGA
ATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4SE
            'U' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCAATCACGTGGTTCAGAAATGCAT
TGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4NQ
            );
#If nucleotide 6 is a 'C', an alternative set of modules for nucleotide 5 MUST be used.
my %nt5altx = ('A' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCTGTTACGTGGTTCAGAAA
TGCATTGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4CYQ
                'C' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCAGTTACGTGGTTCGCAAATGCATTGA
ATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4SYR
                'G' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCAGTTACGTGGTTGAGAAATGCATTGA
ATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4SYE
                'U' => 'CATGTCTTGAAGTGTGTGAAAGATCAGAATGGCAATTACGTGGTTCAGAAA
TGCATTGAATGTGTACAGCCCCAGTCTTTGCAATTTATCATCGATGCGTTTAAGGGACAGG', #4NYQ
                );
my %nt6x = ('A' => 'CACGTCCTGTCATTGGCACTACAGATGTATGGCTGCCGTGTTATCCAGAAAGCTCTTGA
GTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGC', #3CQ
            'C' => 'CACGTCCTGTCATTGGCACTACAGATGTATGGCTCCCGTGTTATCCGCAAAG
CTCTTGAGTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGC', #3SR : See note on required use of nt5alt set.
            'G' => 'CACGTCCTGTCATTGGCACTACAGATGTATGGCTCCCGTGTTATCGAGAAAGCT
CTTGAGTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGC', #3SE
            'U' => 'CACGTCCTGTCATTGGCACTACAGATGTATGGCAACCGTGTTATCCAGAAAGCTCTTGA
GTTTATTCCTTCAGACCAGCAGAATGAGATGGTTCGGGAACTAGATGGC', #3NQ
            );
my %nt7x = ('A' => 'GCTGCCTACCAACTCATGGTGGATGTGTTTGGTTGTTACGTCATTCAGAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGT', #2CQ
            'C' => 'GCTGCCTACCAACTCATGGTGGATGTGTTTGGTAGTTACGTCATTCGCAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGT', #2SR
            'G' => 'GCTGCCTACCAACTCATGGTGGATGTGTTTGGTAGTTACGTCATTGAGAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGT', #2SE
            'U' => 'GCTGCCTACCAACTCATGGTGGATGTGTTTGGTAATTACGTCATTCAGAAGTTCTTTGA
ATTTGGCAGTCTTGAACAGAAGCTGGCTTTGGCAGAACGGATTCGAGGT', #2NQ
            );
my %nt8x = ('A' => 'CATATAATGGAATTTTCCCAAGACCAGCATGGGTCCAGATTCATTCAGCTGAAACTGGAGCGTGCCA
CACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAG', #1SQ
            'C' => 'CATATAATGGAATTTTCCCAAGACCAGCATGGGTCCAGATTCATTCGCCTGAAACTGGAGCGTGCCACACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAG', #1SR
            'G' => 'CATATAATGGAATTTTCCCAAGACCAGCATGGGTCCAGATTCATTGAGCTGAAACTGGA
GCGTGCCACACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAG', #1SE
            'U' => 'CATATAATGGAATTTTCCCAAGACCAGCATGGGAACAGATTCATTCAGCTGAAACTGGA
GCGTGCCACACCAGCTGAGCGCCAGCTTGTCTTCAATGAAATCCTCCAG', #1NQ
);

#5' and 3' sequences for PUF domains alone, or including GFP10/11 fusions.
# Note current GFP 10/11 sequences do NOT include extra sequence for gibson #assembly. Note that BsaI sites currently remain. 
#Join order is 5' 1xx 2xx 3xx 4xx #5xx 6xx 7xx 8xx 3'

my $prime5 = 'ATGCGCAGCCGCCTTTTGGAAGATTTTCGAAACAACCGGTACCCCAATTTACAACTGCGGGAGATTGCCGGA';
my $prime3 = 'CACATCGCAACTCTTCGTAAGTACACCTATGGCAAGCACATTCTGGCCAAGCTGGAGAAGTACTACATGAAGAACGGTGTTGACTTAGGGTAA';
my $GFP105prime = 'ATGGACCTGCCCGACGACCACTACCTGAGCACCCAGACCATCCTGAGCAAGGACCTGAACGACGTGGGCAGCGGCGGCGGCAGCCACATGCGCAGCCGCCTTTTGGAAGATTTTCGAAACAACCGGTACCCCAATTTACAACTGCGGGAGATTGCCGGA';
my $GFP113prime = 'CACATCGCAACTCTTCGTAAGTACACCTATGGCAAGCACATTCTGGCCAAGCTGGAGAAGTACTACATGAAGAACGGTGTTGACTTAGGGGGCAGCGGCGGCGGCAGCGGCGGCGGCAGCGAGAAGAGAGACCACATGGTGCTGCTGGAGTACGTGACCGCCGCCGGCATCACCGACGCCAGCTAA';




# Section 3: Matching the variables to the input to generate a storage variable and printing outputs
#******************************************************************************
if ($base6 eq 'C') {
my @PUF = ($nt1{$base1},' ',$nt2{$base2},' ',$nt3{$base3},' ',$nt4{$base4},' ',$nt5alt{$base5},' ',$nt6{$base6},' ',$nt7{$base7},' ',$nt8{$base8});
my @PUFDNA = ($prime5,$nt8x{$base8},$nt7x{$base7},$nt6x{$base6},$nt5altx{$base5},$nt4x{$base4},$nt3x{$base3},$nt2x{$base2},$nt1x{$base1},$prime3);

my @GG10PUFDNA = ($GFP105prime,$nt8x{$base8},$nt7x{$base7},$nt6x{$base6},$nt5altx{$base5},$nt4x{$base4},$nt3x{$base3},$nt2x{$base2},$nt1x{$base1},$prime3);

my @GG11PUFDNA = ($prime5,$nt8x{$base8},$nt7x{$base7},$nt6x{$base6},$nt5altx{$base5},$nt4x{$base4},$nt3x{$base3},$nt2x{$base2},$nt1x{$base1},$GFP113prime);

print "Your input sequence was: ", $inputsequence, "\n";
print "Base 6 is: ",$base6, "\n";
print "The PUF assembly modules you require are for this sequence are:\n",@PUF, "\n";
print "The DNA sequence for your PUF assembly is: \n", @PUFDNA,"\n";
print "The DNA sequence for your GFP10 PUF assembly is: \n", @GG10PUFDNA,"\n";
print "The DNA sequence for your GFP11 PUF assembly is: \n", @GG11PUFDNA,"\n";
    }
    
elsif ($base6 ne 'C') {
my @PUF = ($nt1{$base1},' ',$nt2{$base2},' ',$nt3{$base3},' ',$nt4{$base4},' ',$nt5{$base5},' ',$nt6{$base6},' ',$nt7{$base7},' ',$nt8{$base8});
my @PUFDNA = ($prime5,$nt8x{$base8},$nt7x{$base7},$nt6x{$base6},$nt5x{$base5},$nt4x{$base4},$nt3x{$base3},$nt2x{$base2},$nt1x{$base1},$prime3);

my @GG10PUFDNA = ($GFP105prime,$nt8x{$base8},$nt7x{$base7},$nt6x{$base6},$nt5x{$base5},$nt4x{$base4},$nt3x{$base3},$nt2x{$base2},$nt1x{$base1},$prime3);

my @GG11PUFDNA = ($prime5,$nt8x{$base8},$nt7x{$base7},$nt6x{$base6},$nt5x{$base5},$nt4x{$base4},$nt3x{$base3},$nt2x{$base2},$nt1x{$base1},$GFP113prime);

print "Your input sequence was: ", $inputsequence, "\n";
print "Base 6 is: ",$base6, "\n";
 print "The PUF assembly modules you require are for this sequence are:\n",@PUF, "\n";
 print "The DNA sequence for your PUF assembly is: \n", @PUFDNA,"\n";
 print "The DNA sequence for your GFP10 PUF assembly is: \n", @GG10PUFDNA,"\n";
print "The DNA sequence for your GFP11 PUF assembly is: \n", @GG11PUFDNA,"\n";
     }

exit;