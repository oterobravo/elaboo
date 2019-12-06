import unittest
import emsa
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

class TestEMSA(unittest.TestCase):
    
    def test_MSA_reading_fasta(self):
        with open("test/align.fasta") as fastatest:
            input_alignment = AlignIO.read(fastatest, "fasta")
        self.assertEqual(len(input_alignment), 2)
        self.assertEqual(len(input_alignment[0]), 4)
    
    def test_MSA_reading_sample_phy(self):
        with open("test/simulated_test_alignment.phy") as phytest:
            input_alignment = AlignIO.read(phytest, "phylip-relaxed")
        self.assertEqual(len(input_alignment), 30)
        self.assertEqual(len(input_alignment[0]), 1000)
    
    def test_alphabet(self):
        emsa1 = emsa.Emsa([SeqRecord(Seq("AT", generic_dna), id="a")])
        emsa2 = emsa.Emsa([SeqRecord(Seq("AU", generic_dna), id="a")])
        emsa3 = emsa.Emsa([SeqRecord(Seq("AF", generic_dna), id="a")])
        emsa1.Alphabet()
        emsa2.Alphabet()
        emsa3.Alphabet()
        self.assertEqual(emsa1.alphabet, "DNA")
        self.assertEqual(emsa2.alphabet, "RNA")
        self.assertEqual(emsa3.alphabet, "PROT")