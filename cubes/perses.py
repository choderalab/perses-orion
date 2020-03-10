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

from floe.api import BooleanParameter, ComputeCube
from orionplatform.mixins import RecordPortsMixin
from snowball.utils.query_utils import get_receptor, get_init_mol
from datarecord import OERecord, Types, OEPrimaryMolField, OEField
from orionplatform.parameters import FieldParameter

from simtk import unit
import perses

class PersesCube(RecordPortsMixin, ComputeCube):
    # Cube documentation.  This documentation for this cube, and all other cubes in this repository, can be converted
    # to html by calling 'invoke docs' from the root directory of this repository.  This documentation will also
    # appear in the Orion Floe editor.
    title = "Perses relative free energy calculations cube"
    classification = [["Free Energy"]]
    tags = ["Ligand", "Protein", "Free Energy", "Perses", "Relative Binding Free Energy", "Alchemical"]
    description = """
    Compute a relative binding free energy using Perses.

    This cube uses Perses to perform a relative alchemical free energy calculation.

    See https://github.com/choderalab/perses for more information about perses.
    """
    uuid = "0fc3762d-35cc-4704-b4c3-820321d98040"

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda9"}, # Can we upgrade this?
        "spot_policy": {"default": "Prohibited"}, # TODO: Figure out how to allow spot policy
        #"spot_policy": {"default": "Allowed"}, # TODO:
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    AVAILABLE_PROTEIN_FORCEFIELDS = [
        'protein.amber14SB.xml',
        'ff99SBildn.xml',
        ]

    AVAILABLE_LIGAND_FORCEFIELDS = [
        'openff-1.0.0',
        'smirnoff99Frosst-1.1.0',
        'gaff-2.11',
        'gaff-1.81',
        ]

    AVAILABLE_SOLVENT_FORCEFIELDS = [
        'tip3p_standard.xml',
        'tip3pfb_standard.xml',
        ]

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default=AVAILABLE_PROTEIN_FORCEFIELDS[0],
        choices=AVAILABLE_PROTEIN_FORCEFIELDS,
        help_text='Force field parameters to be applied to the protein')

    solvent_forcefield = parameter.StringParameter(
        'solvent_forcefield',
        default=AVAILABLE_SOLVENT_FORCEFIELDS[0],
        choices=AVAILABLE_SOLVENT_FORCEFIELDS,
        help_text='Force field parameters to be applied to water and ions')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        default=AVAILABLE_LIGAND_FORCEFIELDS[0],
        choices=AVAILABLE_LIGAND_FORCEFIELDS,
        help_text='Force field to be applied to the ligand')

    temperature = parameter.StringParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=5000,
        help_text="Total number of Perses iterations for the entire floe.")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=250,
        help_text="Number of MD steps per iteration")

    checkpoint_interval = parameter.IntegerParameter(
        'checkpoint_interval',
        default=500,
        help_text="Full checkpoint interval (iterations)")

    timestep = parameter.DecimalParameter(
        'timestep',
        default=4.0,
        help_text="Timestep (fs)")

    hmr = parameter.BooleanParameter(
        'hmr',
        default=True,
        description='On enables Hydrogen Mass Repartitioning (HMR)')

    suffix = parameter.StringParameter(
        'suffix',
        default='perses',
        help_text='Filename suffix for output simulation files')

    nonbonded_method = parameter.StringParameter(
        'nonbonded_method',
        default='PME',
        choices=['PME', 'CutoffPeriodic'],
        help_text='Nonbonded method to use')

    nstates = parameter.IntParameter(
        'nstates',
        default=12,
        help_text='Number of alchemical intermediate states')

    # Ports
    protein_port = RecordInputPort("protein_port", initializer=True)
    reference_ligand_port = RecordInputPort("reference_ligand_port", initializer=True)
    ligand_port = RecordInputPort("ligand_port")

    # Fields
    log_field = OEField("log_field", Types.String)
    DDG_field = OEField("DDG_field", Types.Float)
    dDDG_field = OEField("dDDG_field", Types.Float)

    def begin(self):
        # Retrieve receptor
        self._receptor = self.protein_port.get_value(OEPrimaryMolField())

        # Retrieve reference ligand
        self._reference_ligand = self.reference_ligand_port.get_value(OEPrimaryMolField())

    def process(self, record, port):
        # Make sure we have a molecule defined
        if not record.has_value(OEPrimaryMolField()):
            record.set_value(self.args.log_field, 'Record is missing an input molecule field')
            self.failure.emit(record)
        mol = record.get_value(OEPrimaryMolField())

        # Generate arbitrary 3D coordinates for target ligand
        from openeye import oeomega
        omegaOpts = oeomega.OEOmegaOptions()
        omega = oeomega.OEOmega(omegaOpts)
        ret_code = omega.Build(mol)
        if ret_code != oeomega.OEOmegaReturnCode_Success:
            record.set_value(self.args.log_field, oeomega.OEGetOmegaError(ret_code))
            self.failure.emit(record)
            oechem.OEWriteMolecule(ofs, mol)

        # Prepare input for perses
        # TODO: Use tempfile
        from openeye import oechem
        protein_pdb_filename = 'receptor.pdb'
        with oechem.oemolostream(protein_pdb_filename) as ofs
            oechem.OEWriteMolecule(ofs, self._receptor)
        ligands_mol2_filename = 'ligands.mol2'
        with oechem.oemolostream(ligand_mol2_filename) as ofs
            oechem.OEWriteMolecule(ofs, self._reference_ligand) # molecule 0
            oechem.OEWriteMolecule(ofs, mol) # molecule 1

        # Set up perses calculation
        atom_map = None # generate atom map automatically
        old_ligand_index = 0 # first molecule in SDF
        new_ligand_index = 1 # second molecule in SDF
        forcefield_files = [self.args.protein_forcefield, self.args.solvent_forcefield]
        trajectory_directory = 'lig0to1' # required for perses.analyze.load_simulations
        trajectory_prefix = 'out' # required for perses.analyze.load_simulations
        phases = ['solvent', 'complex']
        from perses.app.relative_setup import RelativeFEPSetup
        fe_setup = RelativeFEPSetup(ligands_mol2_filename, old_ligand_index, new_ligand_index, forcefield_files,
                        phases=phases,
                        protein_pdb_filename=protein_pdb_filename,
                        temperature=self.args.temperature * unit.kelvin, pressure=self.args.pressure * unit.atmospheres,
                        small_molecule_forcefield=self.args.small_molecule_forcefield,
                        trajectory_directory=trajectory_directory,
                        trajectory_prefix=trajectory_prefix,
                        nonbonded_method=self.args.nonbonded_method)

        # Analyze the data
        from perses.analyze.load_simulations import Simulation
        simulation = Simulation(0, 1)
        simulation.load_data()

        # Set output molecule information
        record.set_value(self.DDG_field, simulation.comdg)
        record.set_field(self.dDDG_field, simulation.comddg)
        self.success.emit(record)

    # Uncomment this and implement to cleanup the cube at the end of the run
    def end(self):
        # TO DO: Clean up?
        pass

class ParallelPersesCube(ParallelMixin, PersesCube):
    title = "Parallel " + PersesCube.title
