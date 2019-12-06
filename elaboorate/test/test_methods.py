import unittest
import emsa 
from elaboo_parser import elaboo_parser, start_logging
import elaboo_methods as em
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class TestCALCULATE(unittest.TestCase):

	def test_use_external_file(self):
		input_alignment = emsa.read_alignment("test/simulated_test_alignment.phy", "phylip-relaxed")
		array_file = "test/array.txt"
		ts = em.calculate(input_alignment, 2, "1", taxa_stats_file = array_file)
		ts2 = em.calculate(input_alignment, 2, 1)
		self.assertEqual(len(ts), 3)
		self.assertEqual(len(ts2), 30)
	
	def test_output_is_array(self):
		emsa1 = emsa.Emsa([SeqRecord(Seq("ATCG", generic_dna), id="a"),
			SeqRecord(Seq("ATGG", generic_dna), id="b"),
			SeqRecord(Seq("ATGA", generic_dna), id="c"),
			SeqRecord(Seq("ATCT", generic_dna), id="d"),
			SeqRecord(Seq("ATTT", generic_dna), id="e")])
		ts = em.calculate(emsa1, 2, 1)
		self.assertEqual(len(ts), 5)
		self.assertEqual(len(ts[1]), 9)

