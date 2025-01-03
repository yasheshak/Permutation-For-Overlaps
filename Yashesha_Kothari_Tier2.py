import sys
import numpy as np

def load_bed(file_path):
    """Load BED file into a dictionary of lists of tuples (chromosome, start, end)."""
    bed_data = {}
    with open(file_path) as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3] #strip/split by whitespace of first 3 fields
            if chrom not in bed_data:
                bed_data[chrom] = [] #initialize an empty list if chromosome not seen before
            bed_data[chrom].append((int(start), int(end))) #store in tuple form and add to list for that chromosome
    return bed_data

def load_fai(file_path):
    """Load FASTA index (.fai) file to get chromosome lengths, and stores in a dictionary."""
    genome_data = {}
    with open(file_path) as f:
        for line in f:
            chrom, length = line.strip().split()[:2] #strip/split by whitespace of first 2 fields
            genome_data[chrom] = int(length) #convert length to integer and store in genome_data dictionary
    return genome_data

def merge_intervals(intervals):
    """Merge overlapping intervals within a set of intervals on the same chromosome."""
    if not intervals:
        return []
    
    # Sort intervals by start position only (chromosome is implied by structure)
    intervals = sorted(intervals, key=lambda x: x[0]) 
    
    merged = []
    current_start, current_end = intervals[0] #initialize start and end for merging

    #iterate through the rest of the intervals
    for i in range(1, len(intervals)):
        start, end = intervals[i]
        
        #If intervals overlap/touch, merge them
        if start <= current_end:
            current_end = max(current_end, end) #merge by extending the end to the max of the two
        else:
            # If they don't overlap, add the current interval to the merged list
            merged.append((current_start, current_end))
            current_start, current_end = start, end #move on to the next interval
    
    # Add the last interval after the loop
    merged.append((current_start, current_end))
    
    return merged

def calculate_overlap(set_a, set_b):
    """Calculate the number of overlapping bases between two sets of merged ranges.
    Takes two sets of intervals set_a and set_b, for each chromosome (chromo, start, end), 
    merges any overlapping intervals within each set and calculates the total number of bases where they overlap.
    """
    overlap = 0
    
    # Iterate over chromosomes in set_a
    for chrom in set_a:
        if chrom not in set_b:
            continue #if chrom not in set_b, skip it
        
        # Merge intervals within each set for the current chromosome
        set_a_merged = merge_intervals(set_a[chrom])
        set_b_merged = merge_intervals(set_b[chrom])

        i, j = 0, 0 #initialize pointers for both sets of merged intervals

        #Loop through both sets of merged intervals for the chromosome
        while i < len(set_a_merged) and j < len(set_b_merged):
            start_a, end_a = set_a_merged[i] #extract start, end for current intervals in both
            start_b, end_b = set_b_merged[j]

            if end_a <= start_b:
                i += 1 #if interval A ends before B starts, move to next in A
            elif end_b <= start_a:
                j += 1 #if interval B ends before A starts, move to next in B
            else:
                # There is an overlap and calculate number of bases overlapped
                overlap += max(0, min(end_a, end_b) - max(start_a, start_b))
                
                # Move the interval with the smallest end
                if end_a < end_b:
                    i += 1 #interval A finishes first, move to the next in set A
                else:
                    j += 1 #interval B finishes first, move to next in set B

    return overlap

def permute_ranges(set_a, genome_data):
    """Randomly shifts each genomic interval from set_a to a new position on the same chromosome, maintaining length.
    The new start position is selected randomly within the bounds of the chromosome, ensuring the interval does not exceed chromosome length.
    """
    permuted_set = {} #initialize dictionary for permuted intervals
    
    #iterate through each chromosome in set_a
    for chrom in set_a:
        if chrom not in genome_data:
            continue #skip if chromosome not in genome_data
        
        permuted_set[chrom] = [] #initialize empty list for this chromosome
        chrom_length = genome_data[chrom] #get chromosome length from genome_data

        #iterate through each interval in set_a for the chromosome
        for start, end in set_a[chrom]:
            length = end - start #calculate length of current interval

            #randomly select new start position within bounds of chromosome
            new_start = np.random.randint(0, chrom_length - length + 1)
            permuted_set[chrom].append((new_start, new_start + length)) #append new interval to permuted set
    
    return permuted_set

def permutation_test(set_a, set_b, genome_data, num_permutations=10000):
    """Perform permutation test to calculate p-value, accounting for chromosome structure.
    Permutation test compares the observed overlap between set_a and set_b to overlaps calculated from randomly permuted versions of set_a.
    """
    #actual overlap between original set_a and set_b
    observed_overlap = calculate_overlap(set_a, set_b)
    permuted_overlaps = []

    #perform the permutations for num_permutations times
    for _ in range(num_permutations):
        permuted_a = permute_ranges(set_a, genome_data) #permute intervals in set_a by shifting positions randomly within chromosome bounds
        permuted_overlap = calculate_overlap(permuted_a, set_b) #calculate overlap between permuted set_a and set_b
        permuted_overlaps.append(permuted_overlap) #store result

    permuted_overlaps = np.array(permuted_overlaps) #convert to numpy array for easy calculation
    p_value = np.mean(permuted_overlaps >= observed_overlap) #calculate p-value based on permuted overlaps

    return observed_overlap, p_value

def main():
    if len(sys.argv) < 4:
        print("Usage: python Firstname_Lastname_Tier2.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]")
        return #check if correct number of arguments provided
    
    #extract file paths from command line arguments
    path_set_a = sys.argv[1]
    path_set_b = sys.argv[2]
    path_genome = sys.argv[3]
    
    # Optional number of permutations argument, default to 10,000 if not provided
    num_permutations = 10000
    if len(sys.argv) == 5:
        num_permutations = int(sys.argv[4]) #use provided number of permutations
    
    # Load the BED files and the genome index (FAI) file
    set_a = load_bed(path_set_a) #load genomic intervals for set A
    set_b = load_bed(path_set_b) #load genomic intervals for set B
    genome_data = load_fai(path_genome) #load chromosome lengths from genome index

    # Perform the permutation test
    observed_overlap, p_value = permutation_test(set_a, set_b, genome_data, num_permutations)

    # Print the result in the required format
    print(f"Number of overlapping bases observed: {observed_overlap}, p value: {p_value:.4f}")

if __name__ == "__main__":
    main()
