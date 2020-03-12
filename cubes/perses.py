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

from floe.api.cubes import ComputeCube, ParallelMixin
from floe.api.parameter import IntegerParameter, StringParameter, DecimalParameter, BooleanParameter
from orionplatform.mixins import RecordPortsMixin
from datarecord import OERecord, Types, OEPrimaryMolField, OEField
from orionplatform.parameters import FieldParameter
from orionplatform.ports import RecordInputPort

import yaml

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
        'amber/protein.ff14SB.xml',
        'amber/ff99SBildn.xml',
        ]

    AVAILABLE_LIGAND_FORCEFIELDS = [
        'openff-1.0.0',
        'smirnoff99Frosst-1.1.0',
        'gaff-2.11',
        'gaff-1.81',
        ]

    AVAILABLE_SOLVENT_FORCEFIELDS = [
        'amber/tip3p_standard.xml',
        'amber/tip3pfb_standard.xml',
        ]

    protein_forcefield = StringParameter(
        'protein_forcefield',
        default=AVAILABLE_PROTEIN_FORCEFIELDS[0],
        choices=AVAILABLE_PROTEIN_FORCEFIELDS,
        help_text='Force field parameters to be applied to the protein')

    solvent_forcefield = StringParameter(
        'solvent_forcefield',
        default=AVAILABLE_SOLVENT_FORCEFIELDS[0],
        choices=AVAILABLE_SOLVENT_FORCEFIELDS,
        help_text='Force field parameters to be applied to water and ions')

    ligand_forcefield = StringParameter(
        'ligand_forcefield',
        default=AVAILABLE_LIGAND_FORCEFIELDS[0],
        choices=AVAILABLE_LIGAND_FORCEFIELDS,
        help_text='Force field to be applied to the ligand')

    temperature = StringParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    temperature = DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    n_iterations = IntegerParameter(
        'n_iterations',
        default=5000,
        help_text="Total number of iterations.")

    n_steps_per_iteration = IntegerParameter(
        'n_steps_per_iteration',
        default=250,
        help_text="Number of MD steps per iteration")

    checkpoint_interval = IntegerParameter(
        'checkpoint_interval',
        default=500,
        help_text="Full checkpoint interval (iterations)")

    timestep = DecimalParameter(
        'timestep',
        default=4.0,
        help_text="Timestep (fs)")

    hmr = BooleanParameter(
        'hmr',
        default=True,
        description='On enables Hydrogen Mass Repartitioning (HMR)')

    suffix = StringParameter(
        'suffix',
        default='perses',
        help_text='Filename suffix for output simulation files')

    nonbonded_method = StringParameter(
        'nonbonded_method',
        default='PME',
        choices=['PME', 'CutoffPeriodic'],
        help_text='Nonbonded method to use')

    solvent_padding = DecimalParameter(
        'solvent_padding',
        default=9.0,
        help_text="Solvent padding (A)")

    n_states = IntegerParameter(
        'n_states',
        default=11,
        help_text='Number of alchemical intermediate states')

    # Ports
    protein_port = RecordInputPort("protein_port", initializer=True)
    reference_ligand_port = RecordInputPort("reference_ligand_port", initializer=True)
    target_ligands_port = RecordInputPort("target_ligands_port")

    # Fields
    log_field = OEField("log_field", Types.String)
    DDG_field = OEField("DDG_field", Types.Float)
    dDDG_field = OEField("dDDG_field", Types.Float)

    def begin(self):
        # Retrieve receptor
        mols = [record.get_value(OEPrimaryMolField()) for record in self.protein_port]
        if len(mols) != 1:
            raise Exception(f'{len(mols)} molecules found on protein_port')
        self._receptor = mols[0]
        self.log.info(f'Receptor: {self._receptor.NumAtoms()}')

        # Retrieve reference ligand
        mols = [record.get_value(OEPrimaryMolField()) for record in self.reference_ligand_port]
        if len(mols) != 1:
            raise Exception(f'{len(mols)} molecules found on reference_ligand_port')
        self._reference_ligand = mols[0]
        self.log.info(f'Reference ligand: {self._reference_ligand.NumAtoms()}')

        # Create YAML file
        self.log.info('Creating YAML file...')
        setup_options = dict()
        setup_options['phases'] = ['solvent', 'complex']
        setup_options['protein_pdb'] = 'receptor.pdb'
        setup_options['ligand_file'] = 'ligands.sdf'
        setup_options['old_ligand_index'] = 0
        setup_options['new_ligand_index'] = 1
        setup_options['forcefield_files'] = [self.args.protein_forcefield, self.args.solvent_forcefield]
        setup_options['temperature'] = self.args.temperature
        setup_options['pressure'] = self.args.pressure
        setup_options['small_molecule_forcefield'] = self.args.ligand_forcefield
        setup_options['atom_expression'] = ['IntType']
        setup_options['n_steps_per_move_application'] = self.args.n_steps_per_iteration
        setup_options['fe_type'] = 'repex'
        setup_options['checkpoint_interval'] = self.args.checkpoint_interval
        setup_options['n_cycles'] = self.args.n_iterations
        setup_options['n_states'] = self.args.n_states
        setup_options['n_equilibration_iterations'] = 0
        setup_options['trajectory_directory'] = 'lig0to1'
        setup_options['trajectory_prefix'] = 'out'
        setup_options['atom_selection'] = 'not water'
        setup_options['timestep'] = self.args.timestep
        setup_options['solvent_padding'] = self.args.solvent_padding
        setup_options['save_setup_pickle_as'] = 'out.pkl'

        self.log.info('Writing YAML file...')
        self.yaml_filename = 'perses.yaml'
        with open(self.yaml_filename, 'w') as output:
            yaml.dump(setup_options, output)
            self.log.info(yaml.dump(setup_options))

    def process(self, record, port):
        # Make sure we have a molecule defined
        if not record.has_value(OEPrimaryMolField()):
            record.set_value(self.args.log_field, 'Record is missing an input molecule field')
            self.failure.emit(record)
        mol = record.get_value(OEPrimaryMolField())

        # Report which compound we are processing
        from openeye.oechem import OEMolToSmiles
        smiles = OEMolToSmiles(mol);
        self.log.info(f"Processing compound {smiles}")

        # Generate arbitrary 3D coordinates for target ligand
        from openeye import oeomega
        omegaOpts = oeomega.OEOmegaOptions()
        omega = oeomega.OEOmega(omegaOpts)
        ret_code = omega.Build(mol)
        if ret_code != oeomega.OEOmegaReturnCode_Success:
            record.set_value(self.args.log_field, oeomega.OEGetOmegaError(ret_code))
            self.failure.emit(record)
            oechem.OEWriteMolecule(ofs, mol)

        # Set up perses calculation
        from perses.app.setup_relative_calculation import getSetupOptions, run_setup, run
        self.log.info(f"Loading setup options...")
        setup_options = getSetupOptions(self.yaml_filename)
        self.log.info(str(setup_options))

        from tempfile import TemporaryDirectory
        import os
        cwd = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            self.log.info(f"Entering temporary directory {tmpdir}")
            os.chdir(tmpdir)

            # Prepare input for perses
            # TODO: Use tempdir in future for filesystem reasons
            self.log.info(f"Writing receptor...")
            from openeye import oechem
            protein_pdb_filename = 'receptor.pdb'
            with oechem.oemolostream(protein_pdb_filename) as ofs:
                oechem.OEWriteMolecule(ofs, self._receptor)
            self.log.info(f"Writing ligands...")
            ligands_sdf_filename = 'ligands.sdf'
            with oechem.oemolostream(ligands_sdf_filename) as ofs:
                oechem.OEWriteMolecule(ofs, self._reference_ligand) # molecule 0
                oechem.OEWriteMolecule(ofs, mol) # molecule 1

            self.log.info(f"Setting up perses calculation...")
            perses_setup = run_setup(setup_options)

            self.log.info(f"Running calculations...")
            run(self.yaml_filename)

            # Analyze the data
            self.log.info(f"Analyzing calculations...")
            from perses.analyze.load_simulations import Simulation
            simulation = Simulation(0, 1)
            simulation.load_data()
            os.chdir(cwd)

        # Set output molecule information
        # TODO: Store trajectory or final snapshots
        self.log.info(f"DDG = {simulation.comdg} +- {simulation.comddg} kcal/mol...")
        record.set_value(self.DDG_field, simulation.comdg)
        record.set_field(self.dDDG_field, simulation.comddg)
        self.success.emit(record)

    # Uncomment this and implement to cleanup the cube at the end of the run
    def end(self):
        # TO DO: Clean up?
        pass

class ParallelPersesCube(ParallelMixin, PersesCube):
    title = "Parallel " + PersesCube.title
