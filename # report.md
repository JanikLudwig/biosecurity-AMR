# report
 All 33 headers say [Staphylococcus aureus strain USFL336 | 1280.4647] — same strain, same genome_id, every single one. What looks like "different genomes" is:

One bacterial genome, split into 33 fragments (contigs).

Why a genome file has multiple > entries: sequencing doesn't read a bacterial chromosome start-to-end in one pass. It reads millions of short overlapping fragments, then software (the assembler) stitches them back together. Repetitive regions, low-coverage spots, and other gaps mean the assembler usually can't reconstruct one single continuous sequence — it reconstructs the genome in pieces instead. Each > entry here is one such piece: scaffold ERS093100SCcontig000001, ...000002, etc. — literally labeled "contig" in the header. Add them up: 33 contigs → 2,918,311 bp total, which is exactly the expected size of one S. aureus genome (~2.8–2.9 Mbp). The biggest contig is 663 kb, the rest taper down to smaller fragments — classic draft-assembly shape.

This is completely normal for public genome databases — a "complete, closed" genome (one circular chromosome, zero gaps) is actually the exception; most deposited assemblies are drafts like this one. It'd only be a real problem if headers named different strains/organisms in the same file, which isn't happening here.