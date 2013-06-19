############################
#            Imports             #
############################

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import sys,os

############################
#          Functions            #
############################

def nice_write(s,ind=60):
    o=''
    lines=len(s)/ind
    if lines==0:
        return s + '\n'
    for l in range(1,lines+1):
        st,end=ind*(l-1),(ind*l)-1
        o= o + s[st:end+1] + '\n'
    o= o + s[end+1:] + '\n'
    return o

def rev_complement(seq):
    # get reverse complement of a sequence
    d= {'A': 'T', 'C': 'G',
        'G': 'C', 'T': 'A',
        'a': 't', 'c': 'g', 
        'g': 'c', 't': 'a'}
    rev_compl=''.join(d.get(n, 'N') for n in seq[::-1])
    return rev_compl

def get_replicons(strain):
    # curl is required... if not present, please
    # sudo apt-get install curl
     os.system('bash get_replicons.sh %s > %s' %(strain,'tmp'))
     replicons=[line.strip()[:-4] for line in open('tmp')]
     os.remove('tmp')
     return replicons

def get_from_genbank(accession,mail='emanuele.bosi@unifi.it'):
    Entrez.email=mail
    handle = Entrez.efetch(db="nucleotide",
                        id=accession,
                        rettype="gbwithparts")
    record = SeqIO.read(handle,"genbank")
    handle.close()
    return record

def get_genes(record):
    #return only genes as list of SeqFeatures
    genes = [Gene(r) for r in record.features if r.type == 'CDS']
    return genes
    
def get_gene_name(feat):
    #return gene name(s)
    name=','.join(feat.qualifiers['locus_tag'])
    return name

def get_whole_genome_seq(record):
    #return whole genome seq
    whole_genome= record.seq.tostring()
    return whole_genome

def get_subseq(tupla,genome_seq,flag=False,b_up=400,b_down=0):
    #given some indices, return the corresponding genomic subseq
     strand=tupla[-1]
     if tupla[0] <= tupla[1]:
          if strand == '+':
                left=max(tupla[0],tupla[1]-b_up)
                right=tupla[1]+b_down
                return whole_seq[left:right]
          else:
                left=tupla[0]-b_down
                right=min(tupla[1],tupla[0]+b_up)
                return rev_complement(whole_seq[left:right])
     if flag == 'wrapped around':
          if strand=='+':
                seq= ( whole_seq[tupla[0]:]
                          +
                          whole_seq[:tupla[1]+b_down]
                          )
                if len(seq) > b_up+b_down:
                     seq=seq[-b_up:]
                return seq
          else:
                seq= (  rev_complement(whole_seq[:tupla[1]])
                     +
                     rev_complement(whole_seq[tupla[0]-b_down:])
                          )
                if len(seq) > b_up+b_down:
                     seq=seq[-b_up:]
                return seq
     if b_down > 0:
        if strand == '+':
            left=tupla[1]
            right=tupla[1]+b_down
            return whole_seq[left:right]
        else:
            left=tupla[0]-b_down
            right=tupla[0]
            return rev_complement(whole_seq[left:right]) 
     return ''

def get_all_upstreams(genes,whole_seq,b_up=400,b_down=0,is_circ=True):
    tot=len(genes)
    if tot == 0:
        return
    elif tot == 1:
        if genes[0].strand == '+':
            upstr = [(genes[0].start-b_up,genes[0].start+b_down,genes[0].name,genes[0].strand)]
        else:
            upstr = [(genes[0].end-b_down,genes[0].end+b_up,genes[0].name,genes[0].strand)]
    else:
        upstr=[  (genes[i-1].end,genes[i].start,genes[i].name,genes[i].strand) if genes[i].strand == '+'
                    else 
                 (genes[i].end,genes[i+1].start,genes[i].name,genes[i].strand)
                    for i in range(-1,tot-1)]
        upstr= upstr[1:] + [upstr[0]]
    if is_circ:
        out=[Upstream(upstr[0][2],upstr[0][1],upstr[0][0],get_subseq(upstr[0],whole_seq,'wrapped around',b_up,b_down),upstr[0][3],b_up,b_down)]
        out+=[Upstream(up[2],up[1],up[0],get_subseq(up,whole_seq,b_up=b_up,b_down=b_down),up[3],b_up,b_down) for up in upstr[1:tot - 1]]
        out+=[Upstream(upstr[-1][2],upstr[-1][1],upstr[-1][0],get_subseq(upstr[-1],whole_seq,'wrapped around',b_up,b_down),upstr[-1][3],b_up,b_down)]
    else:
        out=[Upstream(up[2],up[1],up[0],get_subseq(up,whole_seq,b_up=b_up,b_down=b_down),up[3],b_up,b_down) for up in upstr]
    return out

def export_upstreams(upstreams,nome_out,infos):
    out=open(nome_out,'w')
    organism,molecule=infos
    for up in upstreams:
        name,from_,size,start,end,strand,b_up,b_down=up.infos
        # upstream bases
        upstream = -(size-b_down)
        if upstream > 0:
            upstream=0
        to_write='>%s; upstream from %s to %s; size: %s; feature type: cds; location:%s:%s:%s:%s:%s\n' %(name,upstream,b_down,size,organism,molecule,start,end,strand)
        seq=nice_write(up.seq)
        out.write(to_write)
        out.write(seq)
    

############################
#            Classes            #
############################

class Gene(object):
    #
    def __init__(self,SeqFeature):
        self.name=get_gene_name(SeqFeature)
        self.start=int(SeqFeature.location.start)
        self.end=int(SeqFeature.location.end)
        self.strand='+' if SeqFeature.location.strand == 1 else '-'
    #
#
    
class Upstream(object):
    #
    def __init__(self,name,end,start,seq,strand,b_up,b_down):
        self.name=name
        self.start=start
        self.end  =end
        self.strand=strand
        self.seq=seq
        self.length=len(self.seq)
        self.b_up = int(b_up)
        self.b_down = int(b_down)
        self.infos=[ self.name,
                    -self.length,
                     self.length,
                     self.start,
                     self.end,
                     self.strand,
                     -self.b_up,
                     self.b_down]
                     
    #
#


############################
#             Inputs             #
############################

usage='''python upstream_getter.py ORGANISM bases_up bases_down IS_CIRCULAR
            if organism name is unknown, try using "list_genbank_id.sh"
'''

if __name__=='__main__':
    # do stuff
    args=sys.argv
    try:
        organism, b_up, b_down, is_circ = args[1:5]
        b_up = int(b_up)
        b_down = int(b_down)
        is_circ=bool(is_circ)
    except:
        print usage
        sys.exit()

############################
#              Main              #
############################

if __name__ == '__main__':
	# do stuff
	replicons = get_replicons(organism)
	print 'working on these molecules...'
	for replicon in replicons:
		print replicon
		record = get_from_genbank(replicon)
		genes=get_genes(record)
		if len(genes)==0:
			print 'no genes in', replicon
			continue
		whole_seq=get_whole_genome_seq(record)
		upstreams=get_all_upstreams(genes,whole_seq, b_up, b_down, is_circ)
		if not upstreams:continue
		file_name='%s_%s_upstreams.fasta' %(organism,replicon)
		infos=(organism,replicon)
		export_upstreams(upstreams,file_name,infos)
