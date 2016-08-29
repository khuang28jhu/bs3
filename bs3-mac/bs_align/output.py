try :
    import pysam
except ImportError :
    print "[Error] It seems that you haven't install \"pysam\" package.. Please do it before you run this script."
    exit(-1)

import sys
from bs_align_utils import *
import pdb
BAM = 'bam'
SAM = 'sam'
BS_SEEKER1 = 'bs_seeker1'

formats = [BAM, SAM, BS_SEEKER1]

class outfile:
    def __init__(self, filename, format, chrom_len, cmd_line, suppress_SAM_header):
        self.filename = filename
        self.format = format
        self.chrom_ids = dict((k, i) for i, k in enumerate(sorted(chrom_len)))
        self.to_string = ''
        if format == BS_SEEKER1:
            self.f = open(filename, 'w')
        elif format in [SAM, BAM]:
            header = { 'HD' : { 'VN': '1.0'},
                       'SQ' : [ {'LN' : chrom_len[c], 'SN' : c} for c in sorted(chrom_len) ],
                       'PG' : [ { 'ID' : 1, 'PN' : 'BS Seeker 2', 'CL' : cmd_line} ]
                     }
            f = open(filename, 'w')
	    f.write('@HD\tVN:1.0\tSO:unsorted'+ '\n')
	    for c in sorted(chrom_len):
	    	f.write('@SQ\tSN:' + str(c) +'\tLN:' + str(chrom_len[c])+'\n')
 	    f.write('@PG\tID:1\tPN:BS Seeker 2\tCL:'+cmd_line+'\n')
	    f.close() 
            #self.f = pysam.Samfile(filename, 'w' + ('b' if format == BAM else ('' if suppress_SAM_header else 'h')), header = header)
	#pdb.set_trace()

    def critical(self):
	return self.chrom_ids
    def close(self):
        self.f.close()
    def restart(self):
	self.to_string = ''
    def store(self, qname, N_mismatch, FR, refname, strand, pos, cigar, original_BS, methy, STEVE, rnext = -1, pnext = -1, qual = None, output_genome = None,
              rrbs = False, my_region_serial = -1, my_region_start = 0, my_region_end = 0):

        if self.format == BS_SEEKER1:

            # remove the soft clipped bases from the read
            # this is done for backwards compatibility with the old format
            r_start, r_end, _ = get_read_start_end_and_genome_length(cigar)
            original_BS = original_BS[r_start : r_end]

            if not rrbs:
                #pdb.set_trace()
                self.f.write('%s\t%2d\t%s\t%s%s%s\t%s\t%s\t%s\t%d\n' % (qname, N_mismatch, FR, refname, strand, str(pos+1).zfill(10), output_genome, original_BS, methy, STEVE))
            else:
                #pdb.set_trace()
                self.f.write('%s\t%2d\t%s\t%s%s%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n' % (qname, N_mismatch, FR, refname, strand, str(pos+1).zfill(10), output_genome, original_BS, methy, STEVE, my_region_serial, my_region_start, my_region_end))


        elif self.format == BAM or self.format == SAM:
#	    flag =  0x10 if strand == '-' else 0
#	    cigarr = cigar if strand == '+' else list(reversed(cigar))
 #           rnextt = rnext if rnext == -1 else self.chrom_ids[rnext]
  #          seqq = original_BS if strand == '+' else reverse_compl_seq(original_BS)
	    #pdb.set_trace()
#	    self.to_string = self.to_string +  '\n' + qname + '\t' + str(flag) + '\t' + str(self.chrom_ids[refname]) + '\t' + str(pos) + '\t' + '255' + '\t' +  str(cigarr) + '\t' + str(rnextt) + '\t' + str(pnext) + '\t' + str(seqq) + '\t' + str(qual) + '\t' +  'None' +'\t' +   '[(\'XO\', ' + str(FR) +  '), (\'XS\', '+ str(STEVE) + '), (\'NM\', ' + str(N_mismatch) + '), (\'XM\', ' + str(methy) + '), (\'XG\', ' + str(output_genome) + ')]'
 #           return
	    pdb.set_trace()
            a = pysam.AlignedRead()
            a.qname = qname
            a.seq = original_BS if strand == '+' else reverse_compl_seq(original_BS)
            a.flag =  0x10 if strand == '-' else 0
            a.tid = self.chrom_ids[refname]
            a.pos = pos
            a.mapq = 255
            a.cigar = cigar if strand == '+' else list(reversed(cigar))
            a.rnext = rnext if rnext == -1 else self.chrom_ids[rnext]
            a.pnext = pnext
            a.qual= qual
            if rrbs:
                a.tags = (('XO', FR),
                          ('XS', STEVE),
                          ('NM', N_mismatch),
                          ('XM', methy),
                          ('XG', output_genome),
                          ('YR', my_region_serial),
                          ('YS', my_region_start),
                          ('YE', my_region_end)
                          )

            else:
                a.tags = (('XO', FR),
                          ('XS', STEVE),
                          ('NM', N_mismatch),
                          ('XM', methy),
                          ('XG', output_genome))
            #pdb.set_trace()
            self.f.write(a)

    def print_string(self):
	with open(self.filename , 'w') as fh:
		fh.write(self.to_string)

    def store2(self, qname, flag, N_mismatch, FR, refname, strand, pos, cigar, original_BS, methy, STEVE, rnext = -1, pnext = -1, qual = None, output_genome = None,
              rrbs = False, my_region_serial = -1, my_region_start = 0, my_region_end = 0):

        if self.format == BS_SEEKER1:

            # remove the soft clipped bases from the read
            # this is done for backwards compatibility with the old format
            r_start, r_end, _ = get_read_start_end_and_genome_length(cigar)
            original_BS = original_BS[r_start : r_end]

            if rrbs:
                self.f.write('%s\t%2d\t%s\t%s%s%s\t%s\t%s\t%s\t%d\n' % (qname, N_mismatch, FR, refname, strand, str(pos+1).zfill(10), output_genome, original_BS, methy, STEVE))
            else:
                self.f.write('%s\t%2d\t%s\t%s%s%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n' % (qname, N_mismatch, FR, refname, strand, str(pos+1).zfill(10), output_genome, original_BS, methy, my_region_serial, my_region_start, my_region_end, STEVE))


        elif self.format == BAM or self.format == SAM:
#	    flag =  0x10 if strand == '-' else 0
#	    cigarr = cigar if strand == '+' else list(reversed(cigar))
#	    rnextt = rnext if rnext == -1 else self.chrom_ids[rnext]
#	    seqq = original_BS if strand == '+' else reverse_compl_seq(original_BS)
#	    #pdb.set_trace()
#	    self.to_string = self.to_string + '\n' + qname + '\t' + str(flag) + '\t' + str(self.chrom_ids[refname]) + '\t' + str(pos) + '\t' + '255' + '\t' +  str(cigarr) + '\t' + str(rnextt) + '\t' + str(pnext) + '\t' + str(seqq) + '\t' + str(qual) + '\t' +  'None' +'\t' +   '[(\'XO\', ' + str(FR) +  '), (\'XS\', '+ str(STEVE) + '), (\'NM\', ' + str(N_mismatch) + '), (\'XM\', ' + str(methy) + '), (\'XG\', ' + str(output_genome) + ')]'
 #           return 
	    a = pysam.AlignedRead()
            a.qname = qname
            a.seq = original_BS if strand == '+' else reverse_compl_seq(original_BS)
            a.flag = flag
            a.tid = self.chrom_ids[refname]
            a.pos = pos
            a.mapq = 255
            a.cigar = cigar if strand == '+' else list(reversed(cigar))
            a.rnext = rnext if rnext == -1 else self.chrom_ids[rnext]
            a.pnext = pnext
            a.qual= qual
            if rrbs:
                a.tags = (('XO', FR),
                          ('XS', STEVE),
                          ('NM', N_mismatch),
                          ('XM', methy),
                          ('XG', output_genome),
                          ('YR', my_region_serial),
                          ('YS', my_region_start),
                          ('YE', my_region_end)
                          )

            else:
                a.tags = (('XO', FR),
                          ('XS', STEVE),
                          ('NM', N_mismatch),
                          ('XM', methy),
                          ('XG', output_genome))

            self.f.write(a)
