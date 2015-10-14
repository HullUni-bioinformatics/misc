### misc ###
def get_two_senses_strs(record):
    from Bio import SeqIO
    """
    record: a SeqRecord object
    returns a list with two strings: the sequence and the revcomp sequence
    """
    return [str(record.seq), str(record.seq.reverse_complement())]


def execute_cline(cline):
    from subprocess import Popen,PIPE
    p = Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    err, out = p.communicate()
    return err, out

def sed_subset_fastq_gz(full_filename, subset, subset_filename):
    
    values = [full_filename, subset*4, subset_filename]
    
    cline = "zcat {0[0]} | sed -n 1,{0[1]}p | gzip > {0[2]}".format(values)
    
    print "A subset of %i sequences"%(values[1]/4)
    print "from the file %s"%values[0]
    print "will be written to %s"%values[2]
    print "using the command"
    print cline
    
    err, out = execute_cline(cline)
    
    print "stdout:\n%s"%out
    print "stderr:\n%s"%err
    
def parse_nhmmer_tblout(filname):
    lines = [l for l in open(filname,'r').readlines() if not l.startswith('#')]
    matches = []
    for l in lines:
        target_name,accession,query_name,accession,hmmfrom,\
        hmm_to,alifrom,ali_to,envfrom,env_to,sq_len,strand,\
        E_value,score,bias,description_of_target=l.rstrip().split()
        matches.append({
                'target_name':target_name,
                'accession':accession,
                'query_name':query_name,
                'accession':accession,
                'hmmfrom':hmmfrom,
                'hmm_to':hmm_to,
                'alifrom':alifrom,
                'ali_to':ali_to,
                'envfrom':envfrom,
                'env_to':env_to,
                'sq_len':sq_len,
                'strand':strand,
                'E_value':E_value,
                'score':score,
                'bias':bias,
                'description_of_target':description_of_target
                
            })
    return sorted(matches, key=lambda m: float(score), reverse=True)        

def is_overlapping(coords1, coords2, max_overlapp=0):
    
    """
    Check if two ranges are overlapping. The ranges are each
    a list or tupple including of two values, the begining and end.
    The begining and end are included in the range.
    
    # overlapping ranges
    
    >>> coords1 = (10, 100)
    >>> coords2 = (50, 150)
    >>> is_overlapping(coords1, coords2)
    True
    
    >>> is_overlapping(coords2, coords1)
    True
    
    # Second range is revcomp and input is
    # list instead of tuple
    
    >>> coords1 = [10, 100]
    >>> coords2 = [150, 50]
    >>> is_overlapping(coords1, coords2)
    True
    
    # One range is completely nested in the other
    
    >>> coords1 = (10, 100)
    >>> coords2 = (20, 90)
    >>> is_overlapping(coords1, coords2)
    True
    
    >>> is_overlapping(coords2, coords1)
    True
    
    # The ranges overlap by a single position
    
    >>> coords1 = (10, 100)
    >>> coords2 = (100, 200)
    >>> is_overlapping(coords1, coords2)
    True
    
    >>> is_overlapping(coords2, coords1)
    True
    
    # The ranges do not overlap
    
    >>> coords1 = (10, 100)
    >>> coords2 = (200, 300)
    >>> is_overlapping(coords1, coords2)
    False
    """
    
    # inputs need to be lists or tuples
    assert all([isinstance(coords1,(list,tuple)),
                isinstance(coords2,(list,tuple))])
    
    for lst in [coords1, coords2]:
        
        # inputs need to have two values each
        assert len(lst) == 2
        
        # inputs values have to be all int
        for i in lst:
            assert isinstance(i, int)
    
    if any([coords2[0] <= coords1[0] <= coords2[1]-max_overlapp,
            coords2[0]+max_overlapp <= coords1[1] <= coords2[1],
            coords1[0] <= coords2[0] <= coords1[1]-max_overlapp,
            coords1[0]+max_overlapp <= coords2[1] <= coords1[1]]):
        return True
    else:
        return False
    
class Indel:
    
    def __init__(self, pos, length, indel_type, aln):
        
        """
        Characterizes indels in a top sequence alignment
        compared to a bottom sequence alignment
        pos: the position in the top sequence where the indel starts
        length: the number of '-' characters following pos, either in
        the top or in the bottom sequence
        indel_type: insertion: the '-' characters are in the bottom 
        sequence. deletion: the '-' characters are in the top sequence.
        aln: the alignment object.
        """
        
        self.position = pos
        self.length = length
        self.type = indel_type
        self.ref_alignment = aln
        
    def __str__(self):
        return("start: %i,  length: %i, type: %s, sequence: %s, reference: %s"%(
                self.position,
                self.length,
                self.type,
                self.ref_alignment[0].id,
                self.ref_alignment[1].id
               ))

def is_indel(aln, i, indel_type):
    k = None
    if indel_type == 'insertion':
        k = 1
    elif indel_type == 'deletion':
        k = 0
    if aln[k,i] == '-' and (i == 0 or not aln[k,i-1] == '-'):
        # find pos in top sequence
        sub_top_seq_length=len(aln[0,:i])-str(aln[0,:i]).count('-')
        # find insertion length (number of '-' in bottom sequence)
        # or
         # find deletion length (number of '-' in top sequence)
        length = 0
        j = i
        p = aln[k,j]
        while p == '-':
            length += 1
            j += 1
            p = aln[k,j]
        return Indel(sub_top_seq_length-1, length, indel_type, aln)
    else:
        return False
        
def is_insertion(aln, i):
    return is_indel(aln, i, 'insertion')

def is_deletion(aln, i):
    return is_indel(aln, i, 'deletion')

    
def find_indels(aln):
    indels = []
    for i in range(aln.get_alignment_length()):
        indel = is_insertion(aln, i)
        if indel:
            indels.append(indel)
        else:
            indel = is_deletion(aln, i)
            if indel:
                indels.append(indel)
            
    return indels