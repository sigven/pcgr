
## https://gist.github.com/arq5x/5408712
def make_grantham_dict(grantham_mat_file):

    """
	Citation:   http://www.ncbi.nlm.nih.gov/pubmed/4843792
	Provenance: http://www.genome.jp/dbget-bin/www_bget?aaindex:GRAR740104
	"""
    f = open(grantham_mat_file)
    header = f.next().strip().split('\t')
    idx_to_aa = dict(zip(range(0,len(header)), header))

    grantham_dict = {}
    for line in f:
        fields = line.strip().split('\t')
        from_aa = fields[0]
        
        for idx, score in enumerate(fields):
            if idx == 0:
                continue
            to_aa = idx_to_aa[idx]
            grantham_dict[(from_aa, to_aa)] = score
    
    return grantham_dict