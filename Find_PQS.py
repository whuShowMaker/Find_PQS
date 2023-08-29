import re
import pandas as pd

loop3 = ".{1,3}?"
loop5 = ".{4,5}?"
loop7 = ".{6,7}?"
loop_normal = ".{1,7}?"
loopL_any = ".{8,12}?"
loopL_middle = ".{8,21}?"

G3 = "G{3,}"
G2 = "G{2}"


PQS_loop3 = f"{G3}({loop3}{G3})({loop_normal}{G3}){{2}}|{G3}({loop_normal}{G3})({loop3}{G3})({loop_normal}{G3})|{G3}({loop_normal}{G3}){{2}}({loop3}{G3})"
PQS_loop5 = f"{G3}({loop5}{G3})({loop_normal}{G3}){{2}}|{G3}({loop_normal}{G3})({loop5}{G3})({loop_normal}{G3})|{G3}({loop_normal}{G3}){{2}}({loop5}{G3})"
PQS_loop7 = f"{G3}({loop7}{G3})({loop_normal}{G3}){{2}}|{G3}({loop_normal}{G3})({loop7}{G3})({loop_normal}{G3})|{G3}({loop_normal}{G3}){{2}}({loop7}{G3})"
PQS_loop_long = f"{G3}({loopL_any}{G3})({loop_normal}{G3}){{2}}|{G3}({loop_normal}{G3})({loopL_middle}{G3})({loop_normal}{G3})|{G3}({loop_normal}{G3}){{2}}({loopL_any}{G3})"


bulge1G = "(GG[^G]GG?|GG?[^G]GG)"
bulge5G = "(GG[^G]{{1,5}}GG?|GG?[^G]{{1,5}}GG)"
bulge7G = "(GG[^G]{{1,7}}GG?|GG?[^G]{{1,7}}GG)"


PQS_bulge_single_b7 = f"{bulge7G}({loop_normal}{G3}){{3}}|{G3}{bulge7G}({loop_normal}{G3}){{2}}|{G3}({loop_normal}{G3}){bulge7G}({loop_normal}{G3})|{G3}({loop_normal}{G3}){{2}}{bulge7G}"

PQS_bulge_mutiple_b1_n2 = f"{bulge1G}({loop_normal}{bulge1G})({loop_normal}{G3}){{2}}|{bulge1G}({loop_normal}{G3})({loop_normal}{bulge1G})({loop_normal}{G3})|{bulge1G}({loop_normal}{G3}){{2}}({loop_normal}{bulge1G})|{G3}({loop_normal}{bulge1G}){{2}}({loop_normal}{G3})|{G3}({loop_normal}{bulge1G})({loop_normal}{G3})({loop_normal}{bulge1G})|{G3}({loop_normal}{G3})({loop_normal}{bulge1G}){{2}}"
PQS_bulge_mutiple_b1_n3 = f"{bulge1G}({loop_normal}{bulge1G}){{2}}({loop_normal}{G3})|{bulge1G}({loop_normal}{bulge1G})({loop_normal}{G3})({loop_normal}{bulge1G})|{bulge1G}({loop_normal}{G3})({loop_normal}{bulge1G}){{2}}|{G3}({loop_normal}{bulge1G}){{3}}"
PQS_bulge_mutiple_b1_n4 = f"{bulge1G}({loop_normal}{bulge1G}){{3}}"
PQS_bulge_mutiple_b1 = f"{PQS_bulge_mutiple_b1_n2}|{PQS_bulge_mutiple_b1_n3}|{PQS_bulge_mutiple_b1_n4}"


PQS_bulge_mutiple_b5_n2 = f"{bulge5G}({loop_normal}{bulge5G})({loop_normal}{G3}){{2}}|{bulge5G}({loop_normal}{G3})({loop_normal}{bulge5G})({loop_normal}{G3})|{bulge5G}({loop_normal}{G3}){{2}}({loop_normal}{bulge5G})|{G3}({loop_normal}{bulge5G}){{2}}({loop_normal}{G3})|{G3}({loop_normal}{bulge5G})({loop_normal}{G3})({loop_normal}{bulge5G})|{G3}({loop_normal}{G3})({loop_normal}{bulge5G}){{2}}"
PQS_bulge_mutiple_b5_n3 = f"{bulge5G}({loop_normal}{bulge5G}){{2}}({loop_normal}{G3})|{bulge5G}({loop_normal}{bulge5G})({loop_normal}{G3})({loop_normal}{bulge5G})|{bulge5G}({loop_normal}{G3})({loop_normal}{bulge5G}){{2}}|{G3}({loop_normal}{bulge5G}){{3}}"
PQS_bulge_mutiple_b5_n4 = f"{bulge5G}({loop_normal}{bulge5G}){{3}}"
PQS_bulge_mutiple_b5 = f"{PQS_bulge_mutiple_b5_n2}|{PQS_bulge_mutiple_b5_n3}|{PQS_bulge_mutiple_b5_n4}"


PQS_tetrad_G2 = f"{G2}({loop_normal}{G2}){{3}}"

with open("/media/whushowmaker/简单卷/rG4/PQS.fa", "r") as file:
    lines = file.readlines()

pos = ""
strand = ""
results = []
transcript_id = ""
for line in lines:

    if line.startswith(">"):
        pos = line.strip()[1:]
        transcript_id, chrom_pos, strand_info = pos.split("\t")
        chrom, positions = chrom_pos.split(":")
        start, end = map(int, positions.split("-"))
        strand = strand_info.strip()
        continue


    def process_matches(matches, pqs_type, results):
        if matches:
            for match in matches:
                outbase = match.group()
                outpos = match.start()
                #out_start_end = f"{chrom}:{start + match.start()}-{start + match.end()}"
                results.append(
                    [transcript_id, strand, pqs_type, outbase, outpos])


    process_matches(re.finditer(PQS_loop3, line), "loop3", results)
    process_matches(re.finditer(PQS_loop5, line), "loop5", results)
    process_matches(re.finditer(PQS_loop7, line), "loop7", results)
    process_matches(re.finditer(PQS_loop_long, line), "loop_long", results)
    process_matches(re.finditer(PQS_bulge_single_b7, line), "bulge_single_b7", results)
    process_matches(re.finditer(PQS_bulge_mutiple_b1, line), "bulge_mutiple_b1", results)
    process_matches(re.finditer(PQS_bulge_mutiple_b5, line), "bulge_mutiple_b5", results)
    process_matches(re.finditer(PQS_tetrad_G2, line), "tetrad_G2", results)

df = pd.DataFrame(results, columns=["Transcript_ID", "Strand", "PQS_Type", "Motif_Sequence", "Motif_Position"])
df.to_excel("/media/whushowmaker/简单卷/rG4/PQS.xlsx",index_label=False,index=None)