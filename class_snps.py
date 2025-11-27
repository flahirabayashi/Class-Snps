import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq

def load_gff_annotations(gff_file):
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "CDS":
                continue

            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]

            attributes = cols[8]
            gene_id = "unknown"
            for f2 in attributes.split(";"):
                if f2.startswith("ID="):
                    gene_id = f2.replace("ID=", "")

            genes.append({
                "start": start,
                "end": end,
                "strand": strand,
                "gene": gene_id
            })

    return genes


def find_gene(position, annotations):
    for gene in annotations:
        if gene["start"] <= position <= gene["end"]:
            return gene
    return None


def classify_snp(ref_base, alt_base, position, gene, genome_seq):
    if alt_base == ".":
        return "indel"

    pos0 = position - 1

    if gene["strand"] == "+":
        offset = pos0 - (gene["start"] - 1)
        codon_index = offset // 3
        base_in_codon = offset % 3
        codon_start = gene["start"] - 1 + codon_index * 3
        codon_seq = genome_seq[codon_start:codon_start+3]
    else:
        offset = (gene["end"] - 1) - pos0
        codon_index = offset // 3
        base_in_codon = offset % 3
        codon_start = gene["end"] - 1 - codon_index * 3
        codon_seq = genome_seq[codon_start-2:codon_start+1]
        codon_seq = str(Seq(codon_seq).reverse_complement())

    if len(codon_seq) != 3:
        return "codon_error"

    codon_list = list(codon_seq)

    if gene["strand"] == "+":
        codon_list[base_in_codon] = alt_base
    else:
        codon_list[base_in_codon] = str(Seq(alt_base).complement())

    new_codon = "".join(codon_list)

    aa_ref = str(Seq(codon_seq).translate())
    aa_alt = str(Seq(new_codon).translate())

    if aa_ref == aa_alt:
        return "synonymous"
    else:
        return "non-synonymous"


def parse_mummer(mummer_file, gff_file, genome_fasta):
    annotations = load_gff_annotations(gff_file)
    genome_seq = str(SeqIO.read(genome_fasta, "fasta").seq)

    results = []

    with open(mummer_file) as f:
        for line in f:
            if line.startswith("[") or line.strip() == "" or line.startswith("/"):
                continue

            cols = line.split()
            if len(cols) < 3:
                continue

            position = int(cols[0])
            ref_base = cols[1]
            alt_base = cols[2]

            gene = find_gene(position, annotations)

            if not gene:
                results.append((position, ref_base, alt_base, "intergenic", "", "intergenic"))
                continue

            effect = classify_snp(ref_base, alt_base, position, gene, genome_seq)

            results.append((position, ref_base, alt_base, gene["gene"], gene["strand"], effect))

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Classify SNPs (synonymous / non-synonymous) using MUMmer output + GFF3 annotation"
    )
    parser.add_argument("-a", "--annotation", required=True, help="GFF3 annotation file")
    parser.add_argument("-s", "--snps", required=True, help="TXT SNP file (MUMmer output)")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("-o", "--output", required=False, help="Output CSV file")

    args = parser.parse_args()

    results = parse_mummer(args.snps, args.annotation, args.fasta)

    # print to screen
    print("position\tref\talt\tgene\tstrand\teffect")
    for row in results:
        print("\t".join(map(str, row)))

    # write CSV if requested
    if args.output:
        with open(args.output, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["position", "ref", "alt", "gene", "strand", "effect"])
            for row in results:
                writer.writerow(row)

        print(f"\nCSV successfully saved to: {args.output}")


if __name__ == "__main__":
    main()
