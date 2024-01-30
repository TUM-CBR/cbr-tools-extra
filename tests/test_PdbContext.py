import unittest

from cbrextra_metal import PyPdbContext

from .resources import PdbFiles

class TestPdbContextr(unittest.TestCase):

    def test_load_pdb(self):
        PyPdbContext.open_file(PdbFiles.f_3auf)

if __name__ == '__main__':
    unittest.main()