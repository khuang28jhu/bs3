from bs_utils.utils import *
import subprocess
import marshal

def wg_build(fasta_file, build_command, ref_path, aligner):

    # ref_path is a string that contains the directory where the reference genomes are stored with
    # the input Fasta filename appended
    ref_path = os.path.join(ref_path,
                            os.path.split(fasta_file)[1] + '_'+aligner)

    clear_dir(ref_path)
    #---------------------------------------------------------------
    # 1. First get the complementary genome (also do the reverse)
    # 2. Then do CT and GA conversions
    #---------------------------------------------------------------

    open_log(os.path.join(ref_path, 'log'))
    refd = {}
    w_c2t = open(os.path.join(ref_path, 'W_C2T.fa'),'w')
    #c_c2t = open(os.path.join(ref_path, 'ref.fa'),'w')
    chrom_conv = open(os.path.join(ref_path, 'chrom_conv_table'),'w')
    w_g2a = open(os.path.join(ref_path, 'W_G2A.fa'),'w')
    #c_g2a = open(os.path.join(ref_path, 'ref.fa'),'w')
    chrom_num = 0
    
    for chrom_id, chrom_seq in read_fasta(fasta_file):
        serialize(chrom_seq, os.path.join(ref_path, str(chrom_num)))
        marshal.dump(chrom_seq, open(os.path.join(ref_path, str(chrom_num + 1)+'.data'), 'wb'))
        print os.path.join(ref_path, str(int(chrom_num + 1))+'.data')
        #print len(chrom_seq)
        #continue
        refd[str(chrom_num)] = len(chrom_seq)
        chrom_conv.write(chrom_id + ': ' + str(chrom_num) + '\n')
        w_c2t.write('>%s_w_c\n%s\n' % (str(chrom_num), chrom_seq.replace("C","T")))
        w_g2a.write('>%s_c_c\n%s\n' % (str(chrom_num), chrom_seq.replace("G","A")))
        chrom_seq = reverse_compl_seq(chrom_seq)
        w_c2t.write('>%s_w_g\n%s\n' % (str(chrom_num), chrom_seq.replace("C","T")))

        #chrom_seq = reverse_compl_seq(chrom_seq)

        #w_g2a.write('>%s_c_c\n%s\n' % (chrom_id, chrom_seq.replace("C","T")))
        
        w_g2a.write('>%s_c_g\n%s\n' % (str(chrom_num), chrom_seq.replace("G","A")))

        elapsed('Preprocessing '+chrom_id)
	chrom_num += 1

    w_c2t.close()
    w_g2a.close()
    #for outf in [w_c2t, c_c2t, w_g2a, c_g2a]:
    #    outf.close()
    
    serialize_m(refd, os.path.join(ref_path,"refname"))
    elapsed('Genome preprocessing')
    # append ref_path to all elements of to_bowtie
    to_snap = map(lambda f: os.path.join(ref_path, f), ['W_C2T', 'W_G2A'])

    # start bowtie-build for all converted genomes and wait for the processes to finish

    #run_in_parallel([(build_command % { 'fname' : fname }) for fname in to_bowtie])
    subprocess.call(build_command % { 'fname' : os.path.join(ref_path, 'W_C2T') }, shell = True)
    subprocess.call(build_command % { 'fname' : os.path.join(ref_path, 'W_G2A') }, shell = True)
    # delete fasta files of converted genomes
    if aligner != "rmap" :
        #delete_files(f+'.fa' for f in to_bowtie)
        delete_files(os.path.join(ref_path, 'W_C2T.fa'))
	delete_files(os.path.join(ref_path, 'W_G2A.fa'))
    elapsed('Done')

