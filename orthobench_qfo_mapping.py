import os
import re
import sys
import glob
import urllib
import tarfile

"""
Version 1
=========
- Get the QfO ID from "y" in ">x|y|z ..." in QfO fasta files
- Map to Ensembl_PRO ID using their mapping file
- Compare this is Orthobench

Examples or more problematic matches:
QfO - P39061

mapping:
P39061	Ensembl	ENSMUSG00000001435
P39061-1	Ensembl_PRO	ENSMUSP00000101049   (not in orthobench)
P39061-2	Ensembl_PRO	ENSMUSP00000080358   (not in orthobench)

My gene mapping file:
ENSMUSG00000001435.15	ENSMUSP00000072538

Results:
240497 QfO sequences
251378 Orthobench sequences
42904 QfO not found
197593 QfO found
197553 QfO found unique


Version 2
=========
The above suggests:
- Map to gene
- Look in my gene mapping file for that, after stripping the .X (version number?)
  part from the end

Results:
240497 QfO sequences
350451 Orthobench sequences
10472 QfO not found  = 4%
230025 QfO found

A lot of the remainder don't seem to have Ensembl IDs. There are a few where we
get e.g. STRING	9606.ENSP00000328232, and this matches.
What about the FASTA sequence for somthing that has nothing for Ensembl?
I don't find that sequence, I tried one, it was very short

Options: 
- If these are short & not too significant I could ignore these genes completely
- More important question to ask now is what proportion of the Orthobench actual 
  sequences are found (it could be 4-8% since there are 4% more sequences)

Version 3
=========
Map from the orthobench sequences (which I need) to the sequences in QfO. I will 
need to take account of the fact that the mapping is not one-to-one, at least as 
far as potential IDs are concerned.

"""

species_info_qfo="""UP000000589	10090	MOUSE	eukaryota	21989	41750	65484	Mus musculus (Mouse)
UP000008144	7719	CIOIN	eukaryota	16678	631	17317	Ciona intestinalis (Transparent sea squirt) (Ascidia intestinalis)
UP000005640	9606	HUMAN	eukaryota	20600	76467	100605	Homo sapiens (Human)
UP000002277	9598	PANTR	eukaryota	23053	25740	49243	Pan troglodytes (Chimpanzee)
UP000008143	8364	XENTR	eukaryota	17272	37604	57261	Xenopus tropicalis (Western clawed frog) (Silurana tropicalis)
UP000000539	9031	CHICK	eukaryota	18116	9678	27940	Gallus gallus (Chicken)
UP000002280	13616	MONDO	eukaryota	21225	14996	36403	Monodelphis domestica (Gray short-tailed opossum)
UP000002494	10116	RAT	eukaryota	21587	9981	32671	Rattus norvegicus (Rat)
UP000001940	6239	CAEEL	eukaryota	19819	8458	28483	Caenorhabditis elegans
UP000000803	7227	DROME	eukaryota	13811	9639	23628	Drosophila melanogaster (Fruit fly)
UP000002254	9615	CANLF	eukaryota	20649	24748	45715	Canis lupus familiaris (Dog) (Canis familiaris)
UP000000437	7955	DANRE	eukaryota	25698	21387	47544	Danio rerio (Zebrafish) (Brachydanio rerio)""".split("\n")

sp_translate = dict()
for line in species_info_qfo:
    t = line.split()
    sp_translate[t[0]] = t[2]

def get_qfo_filenames(species_info_qfo):
    """
    Require the .fasta and the .idmapping files
    """
    fns = []
    for sp_line in species_info_qfo:
        t = sp_line.split("\t")
        domain = t[3]
        domain = domain[0].upper() + domain[1:]
        fn_base = domain + "/" + t[0] + "_" + t[1]
        fns.append(fn_base + ".fasta")
        fns.append(fn_base + ".idmapping")
    return fns

def print_untar_files(fns):
    print(" ".join(fns))

def get_orthobench_seqs(d, species_exclude):
    accs = []
    for fn in glob.glob(d + "*"):
        if any([sp in fn for sp in species_exclude]):
            continue
        with open(fn, 'r') as infile:
            accs.extend([l[1:].split()[0].rstrip() for l in infile if l.startswith(">")])
    return accs

def get_orthobench_prot_2_gene(fn):
    with open(fn, 'r') as infile:
        r = [l.rstrip().split("\t") for l in infile]
        d = {p:g.split(".")[0] for g,p in r}
    return d

def get_orthobench_genes_from_map_file(fn):
    with open(fn, 'r') as infile:
        genes = [l.split("\t")[0].split(".")[0] for l in infile]
        # genes = [l.split("\t")[0] for l in infile]
    return genes

def get_qfo_seqs(d):
    accs = []
    for fn in glob.glob(d + "*fasta"):
        with open(fn, 'r') as infile:
            try:
                for l in infile:
                    if not l.startswith(">"):
                        continue
                    accs.append(l.rstrip().split()[0].split("|")[1])
            except:
                print(l)
                raise 
    # print("%d seqs" % len(accs))
    return accs

def get_mapping_2_ensembl(direct, id_type, q_rev=False):
    """
    Get the mapping from QfO seqs to orthobench IDs
    Args:
        d - directory containing the *.idmapping files
        id_type - "protein", "gene"
    """
    d = dict()
    for fn in glob.glob(direct + "*.idmapping"):
        with open(fn, 'r') as infile:
            for l in infile:
                t = l.rstrip().split()
                if id_type == "protein":
                    if t[1] == "Ensembl_PRO" or t[1] == "EnsemblGenome_PRO":
                        if q_rev:
                            d[t[2]] = t[0]
                        else:
                            d[t[0]] = t[2]
                elif id_type == "gene":
                    if t[1] == "Ensembl" or t[1] == "EnsemblGenome":
                        if q_rev:
                            d[t[2]] = t[0]
                        else:
                            d[t[0]] = t[2]
    return d

def get_mapping_2_qfo(d, id_type):
    return get_mapping_2_ensembl(d, id_type, q_rev=True)

def sample_not_found(not_found, n = 36):
    print("Sample of sequences not found:")
    N = len(not_found)
    r = int(N/n)
    print(" ".join([not_found[i*r] for i in range(n)]))

def compare(seq_ob, seqs_qfo, d_qfo_2_ob):
    print("%d QfO sequences" % len(seqs_qfo))
    print("%d Orthobench sequences" % len(seq_ob))
    qfo_trans = [d_qfo_2_ob[g] if g in d_qfo_2_ob else None for g in seqs_qfo]
    not_found = [g for g, status in zip(seqs_qfo, qfo_trans) if status is None]
    print("%d QfO not found" % qfo_trans.count(None))
    qfo_trans_unique = set(qfo_trans)
    print("%d QfO found" % (len(qfo_trans) - qfo_trans.count(None)))
    # print(qfo_trans[:10])
    print("%d QfO found unique" % (len(qfo_trans_unique) - 1))
    sample_not_found(not_found)

    print("\nReverse - What orthobench seqs found")
    d_rev = {v:k for k,v in d_qfo_2_ob.items()}
    print("%d in reverse dict compared with %d" % (len(d_rev),len(d_qfo_2_ob))) 
    ob_trans = [d_rev[g] if g in d_rev else None for g in seq_ob]
    not_found = [g for g, status in zip(seq_ob, ob_trans) if status is None]
    print("%d Orthobench not found" % ob_trans.count(None))
    ob_trans_unique = set(ob_trans)
    print("%d Orthobench found" % (len(ob_trans) - ob_trans.count(None)))
    # print(ob_trans[:10])
    print("%d Orthobench found unique" % (len(ob_trans_unique) - 1))
    sample_not_found(not_found)

def find_matches(seq_ob, seqs_qfo, d_ob_2_qfo):
    print("Find Orthobench genes")
    print("%d genes" % len(seq_ob))
    print(seq_ob[:5])
    print(list(d_ob_2_qfo.keys())[:5])
    translated = [d_ob_2_qfo[s] for s in seq_ob if s in d_ob_2_qfo]
    print("%d translated" % len(translated))
    tr_set = set(translated)
    print("%d uniquely translated" % len(tr_set))
    found = tr_set.intersection(seqs_qfo)
    print("%d found" % len(found))
    print("%d not found" % (len(tr_set) - len(found)))

if __name__ == "__main__":
    # Download data 
    url_qfo = "https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/QfO_release_2021_03.tar.gz"
    fn_qfo = url_qfo.split("/")[-1]
    url_orthobench = "https://github.com/davidemms/Open_Orthobench/releases/download/v1.1/BENCHMARKS.tar.gz"
    fn_orthobench = "BENCHMARKS.tar.gz"
    if not os.path.exists(fn_orthobench):
        print("Downloading Orthobench data")
        urllib.urlretrieve(url_orthobench, fn_orthobench)
        tar = tarfile.open(fn_orthobench)
        tar.extractall()
        tar.close()
    else:
        print("Using existing Orthobench data")
    d_fasta_ob = "BENCHMARKS/Input/"
    if not os.path.exists(fn_qfo):
        print("Downloading QfO data")
        #urllib.urlretrieve(url_qfo, fn_qfo)
        tar = tarfile.open(fn_qfo)
        mem = tar.getmembers()
        required = [f for f in mem if re.match("Eukaryota.*\d+\.fasta", f.name)]
        required += [f for f in mem if re.match("Eukaryota.*.idmapping", f.name)]
        tar.extractall(members=required)
        tar.close()
    else:
        print("Using existing QfO data")
    d_qfo = "Eukaryota/"
    url_ob_gene_ids = "https://github.com/davidemms/Open_Orthobench/blob/master/Supporting_Data/Additional_Files/gene_to_transcript_map.txt?raw=true"
    fn_ob_gene_ids = "gene_to_transcript_map.txt"
    urllib.urlretrieve(url_ob_gene_ids, fn_ob_gene_ids)

    # Params
    species_exclude = ["Tetraodon_nigroviridis", ]

    # 0. Print files to untar
    # fns = get_qfo_filenames(species_info_qfo)
    # print_untar_command(fns)

    # 1. Get all orthobench sequences
    seq_ob = get_orthobench_seqs(d_fasta_ob, species_exclude)
    d_prot_2_gene_ob = get_orthobench_prot_2_gene(fn_ob_gene_ids)
    seq_ob = [d_prot_2_gene_ob[p] for p in seq_ob]
    # seq_ob = get_orthobench_genes_from_map_file(fn_ob_gene_ids)

    # 2. Get all qfo sequences
    seqs_qfo = get_qfo_seqs(d_qfo)

    # 3. Get mapping
    # d_qfo_2_ob = get_mapping_2_ensembl(d_qfo, "gene")
    d_ob_2_qfo = get_mapping_2_qfo(d_qfo, "gene")

    # 5. Compare sequence sets
    # if q_per_species:
    #     for seqs_qfo, sp in zip(seqs_qfo_per_species, species):
    #         print("\n" + sp)
    #         compare(seq_ob, seqs_qfo, d_qfo_2_ob)
    # else:
    #     compare(seq_ob, seqs_qfo, d_qfo_2_ob)

    # 6. Find matches for Orthobench sequences
    find_matches(seq_ob, seqs_qfo, d_ob_2_qfo)

