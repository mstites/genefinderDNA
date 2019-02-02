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
    dna = dna[::-1]     # reverses the dna sequence
    dna_r = ""  # initializing the dna reverse complementary as empty string
    while (len(dna_r) < len(dna)):
        dna_r += get_complement(dna[len(dna_r)])     # adds next complementary nucleotide
    return dna_r


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
    while (len(dna_n) != len(dna)):
        codon = dna[len(dna_n):(len(dna_n) + 3)] # next codon to test from the dna sequence
        if(codon == "TGA" or codon == "TAG" or codon == "TAA"): # checks for stop codon
            return dna_n
        else:
            dna_n += codon
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
    >>> find_all_ORFs_oneframe("ATGTTTATGATGATGAAATGT")
    ['ATGTTTATGATGATGAAATGT']

    Check if all valid and in a row
    DO I NEED TO TEST ONES THAT DO NOT START WITH ATG????
    """
    ORF_list = [] # initialzing the new dna to add to the list
    ORF_num = 0 # number of ORFs in the list
    codons_checked = 0

    while codons_checked*3 < len(dna): # checks if all values have been tested
        codon = dna[codons_checked*3:codons_checked*3 + 3]
        if (codon == "ATG"): # checks if the codon is a start codon
            ORF = rest_of_ORF(dna[codons_checked*3:])
            ORF_list.insert(ORF_num, ORF)
            ORF_num += 1
            # adds the number of codons in the ORF to checked, stop codon is handled by codons_checked += 1
            ORF_codons = int(len(ORF)/3)
            codons_checked += (ORF_codons)
        else:
            dna = dna.replace(codon, "") # remove the codon
        # accounts for the checked codon or the stop codon if ORF was found
        codons_checked += 1
        #break
        #Exit if not complete codon (IE only two or one character)
        #Exit if got through whole length
    # check for nested
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
    """
    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=False)
    #doctest.run_docstring_examples(rest_of_ORF, globals(), verbose = True)
