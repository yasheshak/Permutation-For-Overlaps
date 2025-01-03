import sys
import numpy as np

def load_bed(file_path):
    """Load BED file into a list of tuples (chromosome, start, end)"""
    bed_data = []
    with open(file_path) as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3] #strip/split by white space of first 3 fields
            bed_data.append((chrom, int(start), int(end))) #store in tuple form and add to list bed_data
    return bed_data

def load_fai(file_path):
    """Load FASTA index (.fai) file to get chromosome name and lengths, and stores in dict"""
    genome_data = {}
    with open(file_path) as f:
        for line in f:
            chrom, length = line.strip().split()[:2] #stip/split by white space of first 2 fields
            genome_data[chrom] = int(length) #convert length into integer and store as name:value in dictionary 
    return genome_data

def merge_intervals(intervals):
    """Take a list of genomic intervals (tuples with chromosome, start, and end positions)
        Merge any intervals that overlap or are adjacent on the same chromosome. Result is new list of tuples with merged"""

    if not intervals:
        return []
    
    # Sort intervals by chromosome and then start position               
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
   
   #Initialize the merged list and set the first interval as the starting point for merging
    merged = []
    current_chrom, current_start, current_end = intervals[0]
    
    #iterate through rest of the intervals
    for i in range(1, len(intervals)):
        chrom, start, end = intervals[i]
        
        #If current interval is on same chromo and overlaps/touches the previous
        if chrom == current_chrom and start <= current_end:
            current_end = max(current_end, end) #merge the intervalsby exentding the end to the max of the 2
        else:
            # If they don't overlap, add the current interval to the list
            merged.append((current_chrom, current_start, current_end))
            current_chrom, current_start, current_end = chrom, start, end #move on to next interval
    
    # Add the last interval after the loop
    merged.append((current_chrom, current_start, current_end))   
    
    return merged

def calculate_overlap(set_a, set_b):
    """
    Calculate the number of overlapping bases between two sets of genomic intervals (set_a and set_b).
    
    - Takes two lists of intervals (each interval is represented as a tuple: chromosome, start, end).
    - Merges any overlapping intervals within each set to simplify the comparison.
    - Identifies where intervals from both sets overlap on the same chromosome and calculates the total number of overlapping bases.
    """

    # Initialize the variable to store the total number of overlapping bases
    overlap = 0
    
    # Merge overlapping intervals within each set to avoid redundant comparisons
    set_a_merged = merge_intervals(set_a)
    set_b_merged = merge_intervals(set_b)

    # Initialize pointers (i for set A, j for set B) to track the current interval being compared in each set
    i, j = 0, 0

    # Loop through both sets of merged intervals until all intervals in one set have been processed
    while i < len(set_a_merged) and j < len(set_b_merged):
        # Unpack the chromosome, start, and end positions for the current intervals in each set
        chrom_a, start_a, end_a = set_a_merged[i]
        chrom_b, start_b, end_b = set_b_merged[j]

        # Compare chromosomes: move forward in the set with the "earlier" chromosome
        if chrom_a < chrom_b:
            i += 1  # Move to the next interval in set A
        elif chrom_a > chrom_b:
            j += 1  # Move to the next interval in set B
        
        # If intervals are on the same chromosome, check for potential overlap
        else:
            # If the current interval from set A ends before the current interval from set B starts, no overlap, move to next interval in set A
            if end_a <= start_b:
                i += 1
            # If the current interval from set B ends before the current interval from set A starts, no overlap, move to next interval in set B
            elif end_b <= start_a:
                j += 1
            else:
                # Calculate the overlap between the two intervals by taking the intersection of the ranges
                overlap += max(0, min(end_a, end_b) - max(start_a, start_b))
                
                # Move forward in the set where the interval ends first to continue checking other overlaps
                if end_a < end_b:
                    i += 1  # Move to the next interval in set A
                else:
                    j += 1  # Move to the next interval in set B
    
    # Return the total number of overlapping bases
    return overlap


def permute_ranges(set_a, genome_data):
    """
    Randomly shifts each genomic interval in `set_a` to a new position on the same chromosome while preserving its length.
    
    This process is part of the permutation test to assess whether the observed overlap between `set_a` and `set_b` is statistically
    significant. By randomly repositioning intervals within the bounds of each chromosome, we create randomized versions of `set_a`
    that are used to estimate the expected overlap by chance.
    """
    
    # Initialize an empty list to store the permuted intervals
    permuted_set = []

    # Iterate through each interval in `set_a`
    for chrom, start, end in set_a:
        # Calculate the length of the current interval
        length = end - start

        # Retrieve the total length of the current chromosome from the genome data dictionary
        chrom_length = genome_data[chrom]

        # Randomly select a new start position for the interval, ensuring it remains within the chromosome bounds.
        # The new start position must be such that the interval doesn't extend beyond the chromosome length.
        new_start = np.random.randint(0, chrom_length - length + 1)

        # Create the new permuted interval by keeping the original length and chromosome, 
        # but shifting the start and end positions to the newly selected location
        permuted_set.append((chrom, new_start, new_start + length))

    # Return the list of permuted intervals
    return permuted_set


def permutation_test(set_a, set_b, genome_data, num_permutations=10000):
    """
    Perform a permutation test to calculate a p-value, assessing whether the observed overlap between
    two sets of genomic intervals (set_a and set_b) is significantly different from what would be expected by chance.

    - The test compares the actual overlap between set_a and set_b to overlaps obtained by randomly permuting 
      the intervals in set_a, while keeping set_b fixed.
    - The p-value is the proportion of random overlaps (from permuted set_a) that are greater than or equal to the observed overlap.
    
    Parameters:
        set_a (list): List of genomic intervals (tuples of chromosome, start, end) for set A.
        set_b (list): List of genomic intervals for set B.
        genome_data (dict): Information about chromosome lengths to ensure permuted intervals stay within valid bounds.
        num_permutations (int): The number of permutations to perform (default is 10,000).
    """
    
    # Calculate the actual overlap between the original set_a and set_b
    observed_overlap = calculate_overlap(set_a, set_b)
    
    # List to store the overlap values from permuted set_a datasets
    permuted_overlaps = []
    
    # Perform the permutation test by shuffling set_a and calculating overlaps num_permutations times
    for _ in range(num_permutations):
        # Randomly permute the intervals in set_a within the bounds of each chromosome
        permuted_a = permute_ranges(set_a, genome_data)
        
        # Calculate the overlap between the permuted set_a and the fixed set_b
        permuted_overlap = calculate_overlap(permuted_a, set_b)
        
        # Store the result of the overlap for this permutation
        permuted_overlaps.append(permuted_overlap)
    
    # Convert the list of permuted overlaps to a NumPy array for easier statistical calculations
    permuted_overlaps = np.array(permuted_overlaps)
    
    # Calculate the p-value as the proportion of permuted overlaps greater than or equal to the observed overlap
    # This tells us how often a random overlap is at least as large as the actual overlap
    p_value = np.mean(permuted_overlaps >= observed_overlap)
    
    # Return the observed overlap and the computed p-value
    return observed_overlap, p_value


def main():
    if len(sys.argv) < 4:
        print("Usage: python submission.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]")
        return #check
    
    #extract files from command line arguments
    path_set_a = sys.argv[1]
    path_set_b = sys.argv[2]
    path_genome = sys.argv[3]
    
    # Optional number of permutations argument, default to 10,000 if not provided
    num_permutations = 10000
    if len(sys.argv) == 5:
        num_permutations = int(sys.argv[4])
    
    # Load the BED files and the genome index (FAI) file
    set_a = load_bed(path_set_a)
    set_b = load_bed(path_set_b)
    genome_data = load_fai(path_genome)
    
    # Perform the permutation test
    observed_overlap, p_value = permutation_test(set_a, set_b, genome_data, num_permutations)
    
    # Print the result in the required format
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()
