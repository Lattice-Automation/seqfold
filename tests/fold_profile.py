import profile
import pstats

from seqfold.fold import calc_dg, _bulge, _pair, _hairpin, _internal_loop


STATS_FILE = "fold_profile.stats"

SEQ = "'GAAATAGACGCCAAGTTCAATCCGTACTCCGACGTACGATGGAACAGTGTGGATGTGACGAGCTTCATTTATACCCTTCGCGCGCCGGACCGGGGTCCGCAAGGCGCGGCGGTGCACAAGCAATTGACAACTAACCACCGTGTATTCGTTATGGCACCAGGGAGTTTAAGCCGAGTCAATGGAGCTCGCAATACAGAGTT'"
SCRIPT = f"calc_dg({SEQ}, 37.0)"
print(SCRIPT)
profile.run(SCRIPT, STATS_FILE)

STATS = pstats.Stats(STATS_FILE)
STATS.strip_dirs()
STATS.sort_stats("cumulative")
# STATS.print_stats("fold")
STATS.print_stats()

