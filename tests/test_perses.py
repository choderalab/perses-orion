from unittest import TestCase

from datarecord import OERecord
from floe.test import CubeTestRunner

from cubes.mycube import MyCube


class PersesTest(TestCase):

    def setUp(self):
        super().setUp()
        # Here we create the counter cube and test runner
        self.cube = PersesCube("perses")
        self.runner = CubeTestRunner(self.cube)

        # Set parameters to speed up tests
        self.runner.set_parameters(iterations=2)
        self.runner.set_parameters(nsteps_per_iteration=10)

        # This method *must* be called prior to using the cube
        self.runner.start()

    def test_success(self):
        # set_parameters must be called before calling start

        # We pass records to the process function this way
        record = OERecord()
        self.cube.process(record, "ligand_port")

        # Now check the output
        self.assertEqual(self.runner.outputs["success"].qsize(), num_to_send)
        self.assertEqual(self.runner.outputs["failure"].qsize(), 0)

        # This calls the end process on the cube
        self.runner.finalize()

    #def test_false(self):
    #    record = OERecord()
    #    self.runner.set_parameters(switch=False)
    #    self.runner.start()
    #    self.runner.cube.process(record, "intake")
    #    self.runner.finalize()
    #    self.assertEqual(self.runner.outputs["success"].qsize(), 0)
    #    self.assertEqual(self.runner.outputs["failure"].qsize(), 1)
