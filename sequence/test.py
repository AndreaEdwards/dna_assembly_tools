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

    for record in SeqIO.parse(gb_file, 'gb'):
        seq = str(record.seq)
        #print seq
        for feature in record.features:
            feature_dict = feature.qualifiers
            if feature_dict.has_key('gene'):
                for gene in feature_dict['gene']:
                    if gene == gene_name:
                        sequence = feature.extract(record.seq)
                        return sequence #stops at first instance of gene
                    else:
                        error = "sequence not found"
            else:
                pass

    gb_file.close()




def main():
    gb_file = "/Users/Andrea/repositories/dna_assembly_tools/static/txt/mg1655.gb"
    sequence = get_sequence(gb_file, "adhE")
    print sequence

main()
