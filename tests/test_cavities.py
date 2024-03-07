from cbrextra import testutils
from cbrextra.cavities.data import Points
from cbrextra.cavities.main import CavitiesInstance
import pytest

class TestCavitiesInstance:

    def test_get_cavitites(self):

        test_atoms = testutils.open_models(
            Points,
            __file__,
            "cavities",
            "2gvi.atoms.json"
        )

        test_min_volume = 2
        test_max_volume = 125
        expected_cavities = [79]

        for i, atoms in enumerate(test_atoms):
            instance = CavitiesInstance(atoms)
            cavs = instance.get_cavities(test_min_volume, test_max_volume)

            assert len(cavs) == expected_cavities[i], "Got unexpected number of cavitis for a known model"

            for cav in cavs:

                # The algorithm uses cubes instead of spheres,
                # we return speheres. However, for asserting volume
                # we must use the cube formula
                volume = sum(v**3 for v in cav.radii)

                assert volume >= test_min_volume, "Got a cavity with volume below than minimum treshold"
                assert volume <= test_max_volume, "Got a cavity with volume above maximum treshold"