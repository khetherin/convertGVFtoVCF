
class GvfFeatureline:
    """This is to store each field of the GVF line"""
    def __init__(self, seqid, source, so_type, start, end, score, strand, phase, attributes):
        """ Initialise GvfFeature line
        :param seqid: sequence ID i.e. chromosome number
        :param source: source i.e. DGVa
        :param so_type: sequence ontology structural variant type
        :param start: start position
        :param end: end position
        :param score: score
        :param strand: strandedness
        :param phase: phase
        :param attributes: attribute key values
        """
        self.seqid = seqid
        self.source = source
        self.feature_type = so_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    def __str__(self):
        """
        Helper to print variables of the GVF feature line
        :return:line_to_print
        """
        line_to_print = (self.seqid + "\t" + self.source + "\t" + self.feature_type + "\t" + self.start + "\t" +
                         self.end + "\t" + self.score +"\t" + self.strand + "\t" + self.phase + "\t" + self.attributes)
        return line_to_print
