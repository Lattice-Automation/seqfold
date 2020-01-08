import profile
import pstats

from seqfold import calc_dg

# ~7 seconds with python alone
# ~4 seconds with cython extension's output
STATS_FILE = "fold_profile.stats"

SEQ = "'GAAATAGACGCCAAGTTCAATCCGTACTCCGACGTACGATGGAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGGGTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCACCAGGGAGTTTAAGCCGAGTCAATGGAGCTCGCAATACAGAGTT'"
SCRIPT = f"calc_dg({SEQ}, 37.0)"
print(SCRIPT)
profile.run(SCRIPT, STATS_FILE)

STATS = pstats.Stats(STATS_FILE)
STATS.strip_dirs()
STATS.sort_stats("cumulative")
STATS.print_stats()

