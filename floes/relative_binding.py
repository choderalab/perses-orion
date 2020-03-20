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

from floe.api import WorkFloe
from orionplatform.cubes import DatasetReaderCube, DatasetWriterCube
from cubes.perses import PersesCube


# Declare and document floe
job = WorkFloe("Peres Floe", title="Perses Floe")
job.description = (
    "Run a star-map relative free energy calculation"
)
job.classification = [['Molecular Dynamics']]
job.uuid = "155b90cf-90fd-4068-8558-3eac7c01c615"
job.tags = [tag for lists in job.classification for tag in lists]

# Declare Cubes
protein_input_cube = DatasetReaderCube("protein_input_cube")
reference_ligand_input_cube = DatasetReaderCube("reference_ligand_input_cube")
target_ligands_input_cube = DatasetReaderCube("target_ligands_input_cube")
perses_cube = PersesCube("perses_cube")
success_output_cube = DatasetWriterCube("success_output_cube", title='success')
failure_output_cube = DatasetWriterCube("failure_output_cube", title='failure')

# Add cubes to floe
job.add_cube(protein_input_cube)
job.add_cube(reference_ligand_input_cube)
job.add_cube(target_ligands_input_cube)
job.add_cube(perses_cube)
job.add_cube(success_output_cube)
job.add_cube(failure_output_cube)

# Promote parameters
protein_input_cube.promote_parameter(
    "data_in", promoted_name="protein", title="Protein"
)
reference_ligand_input_cube.promote_parameter(
    "data_in", promoted_name="reference_ligand", title="Reference ligand"
)
target_ligands_input_cube.promote_parameter(
    "data_in", promoted_name="target_ligands", title="Target ligands"
)
success_output_cube.promote_parameter(
    "data_out", promoted_name="success", title="Predicted relative binding free energies"
)

failure_output_cube.promote_parameter(
    "data_out", promoted_name="failure", title="Failed relative binding free energy calculations"
)

perses_cube.promote_parameter(
    "n_iterations", promoted_name="n_iterations", title="Total number of iterations"
)

perses_cube.promote_parameter(
    "n_steps_per_iteration", promoted_name="n_steps_per_iteration", title="Number of MD steps per iteration"
)

perses_cube.promote_parameter(
    "protein_forcefield", promoted_name="protein_forcefield", title='Force field parameters to be applied to the protein'
)

perses_cube.promote_parameter(
    "ligand_forcefield", promoted_name="ligand_forcefield", title='Force field to be applied to the ligand'
)
perses_cube.promote_parameter(
    "vacuum_test", promoted_name="vacuum_test", title="If True, just run a quick test in vacuum"
)

protein_input_cube.success.connect(perses_cube.protein_port)
reference_ligand_input_cube.success.connect(perses_cube.reference_ligand_port)
target_ligands_input_cube.success.connect(perses_cube.target_ligands_port)
perses_cube.success.connect(success_output_cube.intake)
perses_cube.failure.connect(failure_output_cube.intake)

if __name__ == "__main__":
    job.run()
