import sys


class Parser:
    def __init__(self, argparse):
        self.argparse = argparse

    def parse(self):
        parser = self.argparse.ArgumentParser(
            description='Fold some proteins.')
        parser.add_argument('--scratch', action='store_true')
        parser.add_argument('--temp', type=int, default=300)
        parser.add_argument('--steps', type=int, default=100000,
                            help="2500000000 should fold the protein")
        parser.add_argument('--writes', type=int,
                            default=1000, help="default is 1000")
        parser.add_argument('--out', type=str, default="/tmp/output.pdb")
        parser.add_argument('--pdb', type=str,
                            default="proteins/villin/1vii.pdb")
        parser.add_argument('--fasta', type=str, default=None)
        args = parser.parse_args(sys.argv[1:])
        return args
