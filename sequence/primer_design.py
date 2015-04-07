#!/opt/local/bin/python2.7


def revcomp(s):
    import string
    complement = string.maketrans('ATCGNatcgnKRSBDMYWVH',
                                  'TAGCNtagcnMYWVHKRSBD')
    return s.translate(complement)[::-1]


def slice_sequence(gb_file, five_prime_res_site="BglII",
                   three_prime_res_site="BamHI"):
    # """return a sequence slice between two restriction sites"""
    from Bio import SeqIO

    gb_file = open(gb_file, 'rU')

    for record in SeqIO.parse(gb_file, 'gb'):
        seq = str(record.seq)
        for feature in record.features:
            for key, value in feature.qualifiers.iteritems():
                global start, end
                if five_prime_res_site in value:
                    start = feature.location.end.position
                if three_prime_res_site in value:
                    end = feature.location.start.position
    return seq[start: end]
    # print start, end
    # print len(seq[start: end])

def get_sequence(gb_file, gene_name):
    """The get_sequence method retrieves the nucleotide sequence of a specfieid
    gene from a genbank file.

    Parameters
    ----------
    gb_file: string
        String equivalent to the path pointing to relevant genbank file.
    gene_name: string
        Shorthand gene name such as "adhE" or "ldhA" that can be found in
        genbank file under the gene category for particular CDS region

    Returns
    -------
    string
        The nucleotide sequence of indicated gene
    """
    from Bio import SeqIO
    gb_file = open(gb_file, 'rU')



def Tm_calc(seq, oligo_conc, salt_conc, RNA=0):
    import Bio
    import Bio.SeqUtils.MeltingTemp
    Tmcalc = Bio.SeqUtils.MeltingTemp.Tm_staluc
    freeenergy = float(Tmcalc(seq, oligo_conc, salt_conc, RNA))
    return round(freeenergy, 2)


class Create_seq_list:
    def __init__(self):
        self.seq_list = []
        self.order = 0

    # def __str__(self):
    #    s = str(self.seq_list)
    #    print s

    def build_seq_list(self, seq):
        self.seq_list.append(seq)

    def return_seq_list(self):
        return self.seq_list

    def paired_seq_list(self):
        all_paired_seqs = []
        while len(self.seq_list) > 0:
            paired_seq_list = []
            for seq in self.seq_list:
                if self.seq_list.index(seq) < 2:
                    paired_seq_list.append(seq)
            if paired_seq_list == []:
                pass
            else:
                for seq in paired_seq_list:
                    if paired_seq_list.index(seq) == 0:
                        self.seq_list.pop(self.seq_list.index(seq))
                    else:
                        all_paired_seqs.append(paired_seq_list)
                        # pass
                        # print all_paired_seqs
                        # print paired_seq_list
            # all_paired_seqs.append(paired_seq_list)
        return all_paired_seqs

    def get_terminal_seqs(self):
        terminal_seqs = []
        terminal_seqs.append(self.seq_list[0])
        terminal_seqs.append(self.seq_list[-1])
        return terminal_seqs


class Get_primer_names:
    def __init__(self):
        self.new_list = []

    def get_terminal_seq_names(self, terminal_seq_names):
        self.terminal_seq_names = terminal_seq_names
        self.new_list.append(str(terminal_seq_names[1]+"_"+terminal_seq_names[0]))
        return self.new_list

    def get_paired_seq_names(self, submitted_seq_names):
        self.submitted_seq_names = submitted_seq_names
        for paired_names in self.submitted_seq_names:
            if len(paired_names) > 2:
                print "error: paired seq name list contains more than 2 elements"
            else:
                self.new_list.append(str(paired_names[0]+"_"+paired_names[1]))
        return self.new_list

    def get_fwd_paired_seq_names(self):
        fwd_primer_names = []
        for name in self.new_list:
            fwd_primer_names.append("fwd_" + name)
        return fwd_primer_names

    def get_rev_paired_seq_names(self):
        rev_primer_names = []
        for name in self.new_list:
            rev_primer_names.append("rev_" + name)
        return rev_primer_names


class Get_primers:
    def __init__(self, terminal_seqs, internal_seqs, overlap, melting_temp):
        self.terminal_seqs = terminal_seqs
        self.internal_seqs = internal_seqs
        self.overlap = overlap
        self.melting_temp = melting_temp
        self.count = 1
        self.oligo_list = []
        self.oligo_overlaps = []
        self.terminal_oligos = []

    def __str__(self):
        primer = "Primer 1: " + fwd_primer_seq1

    # def get_melting_temp(self):
    #    return self.melting_temp

    def find_overlaps(self):
        seq_list = self.internal_seqs
        for paired_seqs in seq_list:
            if len(paired_seqs) > 2:
                print "error: paired seq list contains more than 2 elements"
            else:
                for seq in paired_seqs:
                    if paired_seqs.index(seq) == 0:
                        while len(seq[:self.count]) < len(seq[:self.overlap]):
                            self.count += 1
                            if abs(Tm_calc(seq[-self.count:], 1, 110, 0) - self.melting_temp) < 2:
                                seq1_primer = seq[-self.count:]
                                break
                            else:
                                seq1_primer = seq[-self.count:]
                        print abs(Tm_calc(seq[-self.count:], 1, 110, 0) - self.melting_temp)
                        self.oligo_overlaps.append((seq1_primer, Tm_calc(seq[-self.count:], 1, 110, 0), len(seq1_primer)))
                        self.count = 0
                    if paired_seqs.index(seq) == 1:
                        while len(seq[:self.count]) < len(seq[:self.overlap]):
                            self.count += 1
                            if abs((Tm_calc(seq[-self.count:], 1, 110, 0)) - self.melting_temp) < 2:
                                seq2_primer = seq[:self.count]
                                break
                            else:
                                seq2_primer = seq[:self.count]
                        self.oligo_overlaps.append((seq2_primer, Tm_calc(seq[-self.count:], 1, 110, 0), len(seq2_primer)))
                        self.count = 0
                full_oligo = seq1_primer + seq2_primer
            self.oligo_overlaps.append(full_oligo)
        # print self.oligo_overlaps
        return self.oligo_overlaps

    def find_terminal_overlaps(self):
        for seq in self.terminal_seqs:
            if self.terminal_seqs.index(seq) == 0:
                while len(seq[:self.count]) < len(seq[:self.overlap]):
                    self.count += 1
                    if abs((Tm_calc(seq[-self.count:], 1, 110, 0)) - self.melting_temp) < 2:
                        seq1_primer = seq[:self.count]
                        break
                    else:
                        seq1_primer = seq[:self.count]
                print abs((Tm_calc(seq[-self.count:], 1, 110, 0)) - self.melting_temp)
                self.terminal_oligos.append((seq1_primer, Tm_calc(seq[-self.count:], 1, 110, 0), len(seq1_primer)))
                self.count = 1
            if self.terminal_seqs.index(seq) == 1:
                while len(seq[-self.count:]) < len(seq[-self.overlap:]):
                    self.count += 1
                    if abs((Tm_calc(seq[-self.count:], 1, 110, 0)) - self.melting_temp) < 2:
                        seq2_primer = seq[-self.count:]
                        break
                    else:
                        seq2_primer = seq[-self.count:]
                self.terminal_oligos.append((seq2_primer, Tm_calc(seq[-self.count:], 1, 110, 0), len(seq2_primer)))
                self.count = 0
        full_oligo = seq2_primer + seq1_primer
        self.terminal_oligos.append(full_oligo)
        return self.terminal_oligos


    def fwd_primers(self):
        fwd_primers_to_order = []
        for element in self.terminal_oligos:
            if isinstance(element, str) is True:
                fwd_primers_to_order.append(element)
        for element in self.oligo_overlaps:
            if isinstance(element, str) is True:
                fwd_primers_to_order.append(element)
        return fwd_primers_to_order

    def rev_primers(self):
        rev_primers_to_order = []
        for element in self.terminal_oligos:
            if isinstance(element, str) is True:
                rev_primer = revcomp(element)
                rev_primers_to_order.append(rev_primer)
        for element in self.oligo_overlaps:
            if isinstance(element, str) is True:
                rev_primer = revcomp(element)
                rev_primers_to_order.append(rev_primer)
        return rev_primers_to_order

    def get_oligo_length(self):
        oligo_length = []
        for seq, melting_temp, length in self.oligo_list:
            oligo_length.append(length)
        return oligo_length
