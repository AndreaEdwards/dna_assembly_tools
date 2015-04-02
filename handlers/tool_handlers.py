import os
from os import getcwd
from os import remove
from os.path import join
from os.path import splitext
from random import choice

import tornado.web

from sequence.primer_design import (revcomp, slice_sequence, Tm_calc,
                                    Create_seq_list, Get_primers,
                                    Get_primer_names)


class CPECHandler(tornado.web.RequestHandler):
    def get(self):
        self.render("cpec_tools.html")


class UploadHandler(tornado.web.RequestHandler):
    def post(self):
        seq_list = Create_seq_list()
        seq_names = Create_seq_list()
        ordered_fileinfo_list = []

        for num in range(1, int(self.get_argument("sequences", ""))+1):
            ordered_fileinfo_list.append(num)
            if self.get_argument(str(num) + "-Sequence") == "":
                ordered_fileinfo_list[num-1] = [self.get_argument(str(num) +
                                                "-Name"),
                                                self.get_argument(str(num) +
                                                "-Sequence"),
                                                self.get_argument(str(num) +
                                                "-Upstream Feature"),
                                                self.get_argument(str(num) +
                                                "-Downstream Feature"),
                                                self.request.files[str(num)]]
            else:
                ordered_fileinfo_list[num-1] = [self.get_argument(str(num) +
                                                "-Name"),
                                                self.get_argument(str(num) +
                                                "-Sequence")]


        for list_item in ordered_fileinfo_list:
            if list_item[1] == "":
                upstreamFeature = list_item[2]
                downstreamFeature = list_item[3]
                for fileinfo in list_item[4]:
                    fname = fileinfo['filename']
                    seq_name, extn = os.path.splitext(fname)
                    # print seq_name
                    fbody = fileinfo['body']
                    path = "results/" + fname
                    seq_names.build_seq_list(seq_name)
                    with open(path, "w") as temp_file:
                        temp_file.write(fbody)
                    seq = slice_sequence(path, five_prime_res_site=upstreamFeature,
                                         three_prime_res_site=downstreamFeature)
                    seq_list.build_seq_list(seq)
                    os.remove(path)
            else:
                seq = str(list_item[1].upper())
                print seq
                seq_name = list_item[0]
                seq_list.build_seq_list(seq)
                seq_list.return_seq_list()
                seq_names.build_seq_list(seq_name)


        terminal_seqs = seq_list.get_terminal_seqs()
        terminal_seq_names = seq_names.get_terminal_seqs()
        internal_seqs = seq_list.paired_seq_list()
        internal_seq_names = seq_names.paired_seq_list()
        oligo_labels = Get_primer_names()
        oligo_labels.get_terminal_seq_names(terminal_seq_names)
        oligo_labels.get_paired_seq_names(internal_seq_names)
        fwd_primer_names = oligo_labels.get_fwd_paired_seq_names()
        rev_primer_names = oligo_labels.get_rev_paired_seq_names()
        get_primers = Get_primers(terminal_seqs, internal_seqs, 25, 68)
        get_internal_overlaps = get_primers.find_overlaps()
        print get_internal_overlaps
        get_terminal_overlaps = get_primers.find_terminal_overlaps()
        print get_terminal_overlaps
        fwd_primers = get_primers.fwd_primers()
        rev_primers = get_primers.rev_primers()
        fwd_colorseqs = '\n'.join(colorize(fwd_primer_names, fwd_primers))
        rev_colorseqs = '\n'.join(colorize(rev_primer_names, rev_primers))
        melting_temps = '\n'.join(showTm(get_internal_overlaps, get_terminal_overlaps))
        with open("results/oligo_seqs.txt", "w") as oligo_seqs:
            fwd_primers = str(get_primers.fwd_primers())
            rev_primers = str(get_primers.rev_primers())
            oligo_seqs.write(fwd_primers + rev_primers)
        self.render("cpec_results.html", user=self.get_current_user(),
                    fwd_seqs=fwd_colorseqs, rev_seqs=rev_colorseqs,
                    melting_temps=melting_temps)


def colorize(seqnames, seqlist):
    dictionary = dict(zip(seqnames, seqlist))
    for seqname in dictionary:
        colorseq = ["<div class='row' style='width:100%; margin-left:0px; padding:0px'>",
                    "<div class='col-lg-12' style='margin:0px; padding:0px'>"]
        colorseq.append("<span style='color:%s'>%s</style></br>" % ("#2E9AFE", seqname))
        colorseq.append("<span style='color:%s'>%s</style>" % ("white", dictionary[seqname]))
        colorseq.extend(("</div></div>", "<div class='row' style='width:auto'>",
                             "<div class='col-lg-12' style='margin:0px; padding:0px'>"))
        colorseq.append("</div></div>")
        yield "".join(colorseq)
        print "".join(colorseq)


def showTm(internal_overlaps, terminal_overlaps):
    colorseq = ["<div class='row' style='width:100%; margin-left:0px; margin-top:30px; padding:0px'>", "<div class='col-lg-12' style='margin:0px; padding:0px'>"]
    colorseq.append("<span style='color:%s'>%s</style></br>" % ("white", internal_overlaps))
    colorseq.append("<span style='color:%s'>%s</style>" % ("white", terminal_overlaps))
    colorseq.extend(("</div></div>", "<div class='row' style='width:auto'>",
                         "<div class='col-lg-12' style='margin:0px; padding:0px'>"))
    colorseq.append("</div></div>")
    yield "".join(colorseq)
    print "".join(colorseq)
