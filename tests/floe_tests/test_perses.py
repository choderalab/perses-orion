# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import os
import json
import pytest
from subprocess import check_output

from artemis.wrappers import WorkFloeWrapper, DatasetWrapper, OutputDatasetWrapper
from artemis.test import FloeTestCase
from artemis.decorators import package
from artemis.packaging import OrionTestPackage

from openeye.oechem import oeifstream
from datarecord import OEReadRecords
from datarecord.utils import TemporaryPath

import cubes

PACKAGE_DIR = os.path.dirname(os.path.dirname(cubes.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "test_data", "star-benchmark")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")


test_package = OrionTestPackage(manifest=dict(requirements="requirements_dev.txt"))
# Add the contents of the regular package
test_package.add_directory(PACKAGE_DIR)
# Remove the tests as have different requirements
test_package.remove_directory("tests/")
# Remove tasks.py as it requires invoke
test_package.remove_file("tasks.py")

with TemporaryPath(suffix=".txt") as path:
    results = check_output(["python", "setup.py", "--requires"], cwd=PACKAGE_DIR)
    requirements = json.loads(results.decode())
    with open(path, "w") as ofs:
        for result in requirements:
            # Create a file with orion requirements
            ofs.write("{}\n".format(result))
    # Add the orion requirements requirements
    test_package.add_file(path, dest="requirements.txt")


@pytest.mark.floetest
@package(test_package)
class TestPersesFloe(FloeTestCase):

    def test_simple_run(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "relative_binding.py"),
            run_timeout=1200
        )
        protein_file = DatasetWrapper.get_dataset(os.path.join(FILE_DIR, "Thrombin_protein.pdb"))
        reference_ligand_file = DatasetWrapper.get_dataset(os.path.join(FILE_DIR, "ligand7.sdf"))
        target_ligands_file = DatasetWrapper.get_dataset(os.path.join(FILE_DIR, "ligand0.sdf"))
        success_output_file = OutputDatasetWrapper(extension=".oedb")
        failure_output_file = OutputDatasetWrapper(extension=".oedb")
        workfloe.start(
            {
                "promoted": {
                    "n_iterations" : 20,
                    "n_steps_per_iteration" : 2,
                    "solvent_padding" : 0.0,
                    "protein": protein_file.identifier,
                    "vacuum_test": True, # test ligand in vacuum only
                    "reference_ligand": reference_ligand_file.identifier,
                    "target_ligands": target_ligands_file.identifier,
                    "success": success_output_file.identifier,
                    "failure": failure_output_file.identifier,
                }
            }
        )
        # Faked locally
        self.assertEqual(workfloe.state, "complete")
        # Also faked
        self.assertEqual(
            len(workfloe.reason),
            0,
            "Failed with reason {}".format(workfloe.reason)
        )

        ifs = oeifstream()
        with open(success_output_file.path, "rb") as ifs:
            records = list(OEReadRecords(ifs))
        count = len(records)
        self.assertEqual(count, 10)
