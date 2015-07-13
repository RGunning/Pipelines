from bedtools import Interval, IntervalFile
import argparse



def main(args):
	"""
	Examples of printing each interval in an interval file.
	- Works with BED, GTF and VCF files.
	- Can be uncompressed or GZIP compressed.
	"""

	##########################################################
	# ex1. Report the coordinates of overlap b/w exons and rmsk
	#
	# Equivalent to: intersectBed -a exons -b rmsk
	# Uses:           IntervalFile.all_hits()
	##########################################################
	genes = IntervalFile(args.genefile)
	peaks  = IntervalFile(args.peakfile)

	for gene in genes:
		for peak_hit in peaks.all_hits(gene):
			print "\t".join(str(f) for f in [gene.chrom, peak_hit.o_start, peak_hit.o_end])

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Intersect a peak file with a bed file and report the length of peaks residing in each gene')
	parser.add_argument("-p", "--peak", dest="peakfile",
                    help="read peak FILE", metavar="FILE")
	parser.add_argument("-g", "--gene", dest="genefile",
                    help="read peak FILE", metavar="FILE")
	args = parser.parse_args()
	main(args)

