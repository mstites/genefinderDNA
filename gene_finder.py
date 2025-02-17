# -*- coding: utf-8 -*-
"""
Gene finding tool

@author: Griffith Stites

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    Added this unit test because A and C are only two of the possible four nucleotides, T is another.
    >>> get_complement('T')
    'A'

    Added this unit test because A and C are only two of the possible four nucleotides, G is another.
    >>> get_complement('G')
    'C'
    """
    if(nucleotide == 'A'):
        return 'T'
    if(nucleotide == 'C'):
        return 'G'
    if(nucleotide == 'T'):
        return 'A'
    if(nucleotide == 'G'):
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    The following two unit tests are sufficient because they test different DNA sequence lengths, all four complements, and different orders of compliments.
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna_r_c = ""  # initializing the dna reverse complementary as empty string
    for n in dna[::-1]: # goes through dna in reverse order
        dna_r_c += get_complement(n) # adds next complementary nucleotide
    return dna_r_c


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open readFirst 3ing frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'

    The following unit test tests the third stop codono (TAA) as well as if there are more codons after the stop codon. It is also longer and has codons with many like characters in a row. I also threw another start and stop codon on the end.
    >>> rest_of_ORF("ATGATTTTTCAATAAATGTAA")
    'ATGATTTTTCAA'

    The following unit test tests if there is a dna sequence with no stop codon.
    >>> rest_of_ORF("ATGGCGTGG")
    'ATGGCGTGG'
    """
    dna_n = "" # initializing new dna sequence w/o stop codon as an empty string
    for i in range(0, len(dna), 3): # count i from 0 through the dna length, by 3
        codon = dna[i:i + 3] # next codon to test from the dna sequence
        if(codon == "TGA" or codon == "TAG" or codon == "TAA"): # checks for stop codon
            return dna_n # have reached the stop codon, return
        else:
            dna_n += codon # add the codon to the dna sequence
    return dna_n


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    The following unit test was added to check to reject ORFs that are not multiples of 3.
    There are also items between the stop codon and the second start codon.
    >>> find_all_ORFs_oneframe("ATGTACTAATTTGAGAATGTTTGAGTGGATTTAG")
    ['ATGTAC']

    The following unit test was added to check to reject nested ORFs.
    There are also 3 non-nested ORFs to check for more ORFs.
    >>> find_all_ORFs_oneframe("ATGTTTATGGAGATTTAAATGTTTAGGCAATGAATGGGTGGCTAATTTTTCGTT")
    ['ATGTTTATGGAGATT', 'ATGTTTAGGCAA', 'ATGGGTGGC']

    The following unit test was added to test if there was no stop codon.
    Also if the final characters do not make a codon
    >>> find_all_ORFs_oneframe("ATGTTTATGATGATGAAATG")
    ['ATGTTTATGATGATGAAATG']

    The following unit test was added to test if the sequence did not start with ATG.
    >>> find_all_ORFs_oneframe("GCATGAATGTAGATAGATGTGCCCCATGATTGAC")
    ['ATG']
    """
    ORF_list = []
    ORF_num = 0 # number of ORFs in the list
    codons_checked = 0

    while codons_checked*3 < len(dna): # checks if all values have been tested
        codon = dna[codons_checked*3:codons_checked*3 + 3]
        if (codon == "ATG"): # checks if codon is a start codon
            ORF = rest_of_ORF(dna[codons_checked*3:]) # sends dna with current codon to end
            ORF_list.insert(ORF_num, ORF) # adds the ORF to the list
            ORF_num += 1
            codons_checked += (int(len(ORF)/3)) # adds length of ORF to checked
        codons_checked += 1 # also for stop codon if ORF found

    return ORF_list

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']

    The following unit test was added to test ignoring nested ORFs.
    >>> find_all_ORFs("ATGTTTATGAACTGGTAGCATGTAG")
    ['ATGTTTATGAACTGG', 'ATG']

    The following unit test was added to test multiple ORFs in one frame and not ending in a stop codon.
    >>> find_all_ORFs("ATGCATGAATGTAGATAGATGTGCCCCATGATTGAC")
    ['ATGCATGAATGTAGA', 'ATGTGCCCCATGATTGAC', 'ATGAATGTAGATAGATGTGCCCCA', 'ATG']
    """
    ORF_list = []
    for i in range(3): # checking the offsets 0-2
        ORF_list += find_all_ORFs_oneframe(dna[i:]) # add all the ORFs in the current frame
    return ORF_list

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']

    Test not starting with ATG and with the reverse not having an ORF, also offset.
    >>> find_all_ORFs_both_strands("CCAAAAGATGTACTGCTGG")
    ['ATGTACTGCTGG']

    Test with multiple items on both strands. Test for nesting?
    >>> find_all_ORFs_both_strands("CATCATATGTAGTTAATCATGCGGAGGCATGGG")
    ['ATG', 'ATGCGGAGGCATGGG', 'ATGGG', 'ATGCCTCCGCATGAT', 'ATGATG', 'ATGATTAACTACATA']


    With multiple items on both sides
    """
    ORF_list = find_all_ORFs(dna) # add ORFs in inital strand
    ORF_list += find_all_ORFs(get_reverse_complement(dna)) # add ORFs in reverse_complement strand
    return ORF_list

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'

    Testing if there are two ORFs with the same length. Also testing more than two ORFs.
    >>> longest_ORF("CATCATATGTAGTTAATCATGCGGAGGCATGGG")
    'ATGCGGAGGCATGGG'

    Testing if there is only one ORF.
    >>> longest_ORF("ATGCGA")
    'ATGCGA'

    Testing if there is no ORF.
    >>> longest_ORF("ACGTAAAAAAAA")
    ''
    """
    ORF_list = find_all_ORFs_both_strands(dna) # get the ORFs
    if (len(ORF_list) == 0): # if there are no ORFs
        return "" # returns an empty string
    else:
        return max(ORF_list, key=len) # finds the max length string in the list


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

    dna: a DNA sequence
    num_trials: the number of random shuffles
    returns: the maximum length longest ORF

    We cannot write doctests with expected results for this function because
    the result is randomized, so we cannot predict the outcome. Running the
    function 100 times might result in a different solution that running it
    100 times a seperate time would. A way totest the function is to run the
    function only a few times, printing each result, then checking to make
    sure the function returned the longer one.

    # >>> longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 7)
    # >>> longest_ORF_noncoding("ATGCGA", 3)"""
    longest = 0
    for x in range(num_trials): # runs the number of times of num_trials
        dna_s = shuffle_string(dna) # find a shuffled dna sequence
        # print("dna_s = ", dna_s) # print statement to make sure the string is shuffled
        longest_s = longest_ORF(dna_s) # finds the longest ORF on that
        # print("longest_ORF = ", longest_s) # print statement to check the longest (so we can see if the returned longest is actually the longest)
        if(len(longest_s) > longest): # if the longest ORF in the strand is longer than the recorded longest ORF
            longest = len(longest_s) # set longest to the longest ORF strand, int
    return longest

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'

        Checks when there is only a start codon and a stop codon.
        >>> coding_strand_to_AA("ATGTAA")
        'M'

        Checks when there are multiple sequences (need to get longest)
        >>> coding_strand_to_AA("CATCATATGTAGTTAATCATGCGGAGGCATGGG")
        'MRRHG'
    """
    sequence = longest_ORF(dna) # gets the sequence

    sequence_AA = "" # string where the AA will be added, starting with M for ATG
    for i in range(0, len(sequence), 3): # goes through by steps of 3
        if(i+3 > len(sequence)): # if the next item is not a full codon
            break
        sequence_AA += aa_table[sequence[i:i+3]] # add the amino acid to the sequence

    return sequence_AA


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
        #>>> gene_finder("CATCATATGTAGTTAATCATGCGGAGGCATGGG")
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    AA_list = [] # list of amino acid strings
    for i in range(len(ORFs)): # go through all ORFs
        if (len(ORFs[i]) > threshold): # check if ORF is longer than threshold
            AA_list.append(coding_strand_to_AA(ORFs[i])) # add the current ORF to the AA_list

    return AA_list

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    ##########################################################
    ### Edit this portion to run  a different dna sample #####
    # the final name and location relative to gene_finder.py #
    dna = load_seq("./data/X73525.fa")
    ##########################################################
    results = gene_finder(dna)
    print results
