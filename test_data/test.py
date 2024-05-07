import sys

# Function to check if two intervals overlap
def overlap(interval1, interval2):
    return not (interval1[1] < interval2[0] or interval1[0] > interval2[1])

# Load CDS intervals
cds_intervals = {}
with open(sys.argv[1], "r") as cds_file:
    for line in cds_file:
        transcript, start, end = line.strip().split("\t")
        cds_intervals[transcript] = (int(start), int(end))

# Process the 12-column TSV file
with open(sys.argv[2], "r") as input_file:
    for line in input_file:
        fields = line.strip().split("\t")
        transcript = fields[9]
        interval = (int(fields[10]), int(fields[11]))
        
        # Check if there's an overlap with CDS interval
        if transcript in cds_intervals:
            cds_interval = cds_intervals[transcript]
            overlap_start = max(interval[0], cds_interval[0])
            overlap_end = min(interval[1], cds_interval[1])
            if overlap_start <= overlap_end:
                print("\t".join(fields) + f"\t{overlap_start}\t{overlap_end}")
            else:
                print("\t".join(fields) + "\tNA\tNA")
        else:
            print("\t".join(fields) + "\tNA\tNA")
