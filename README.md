"""
# **Genomic Overlap and Permutation Test Tool**

## **Overview**
This tool provides functionality to compute genomic interval overlaps, merge overlapping intervals, and perform a permutation test to assess the significance of observed overlaps. It can handle chromosome-specific structures and ensures intervals stay within chromosome boundaries.

---

## **Features**

### **Genomic Interval Operations**
1. **Load Genomic Intervals**:
   - Reads BED files to extract genomic intervals (chromosome, start, end).
   - Supports chromosome-specific structures.

2. **Merge Overlapping Intervals**:
   - Consolidates overlapping or adjacent intervals within the same chromosome.

3. **Calculate Overlap**:
   - Computes the total number of overlapping bases between two sets of genomic intervals.

### **Permutation Test**
1. **Random Permutations**:
   - Randomly shifts intervals within chromosome boundaries while maintaining their lengths.

2. **Statistical Analysis**:
   - Performs a permutation test to compare observed overlaps to randomized overlaps.
   - Outputs a p-value to determine the significance of the observed overlap.

---

## **Dependencies**
- **Python 3.x**
- **NumPy** (`pip install numpy`)

---

## **Input Files**
1. **SetA.bed**:
   - Contains the first set of genomic intervals.
2. **SetB.bed**:
   - Contains the second set of genomic intervals.
3. **Genome.fa.fai**:
   - Contains chromosome names and lengths in FASTA index format.

---

## **Usage**

### **Input Arguments**:
1. `path/to/SetA.bed`: Path to the first BED file.
2. `path/to/SetB.bed`: Path to the second BED file.
3. `path/to/genome.fa.fai`: Path to the genome index file.
4. `[num_permutations]`: Optional number of permutations (default is 10,000).

### **Execution**:
```bash
python genomic_overlap.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]
