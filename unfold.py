import sys


class Unfold:
    def __init__(self, args):
        self.args = args

    def unfold(self):
        if self.args.scratch:
            # unfolded protein
            if self.args.fasta is not None:
                fasta = self.args.fasta
            else:
                protein_fasta = "proteins/villin/1vii.fasta"
                fasta = open(protein_fasta).read().split("\n")[1]
            print("folding %s" % fasta)
            from lib import write_unfolded
            write_unfolded(fasta, "/tmp/unfolded.pdb")
            pdb = PDBFile("/tmp/unfolded.pdb")
        else:
            # already folded protein
            pdb = PDBFile(self.args.pdb)
        return pdb
