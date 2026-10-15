import json
import math
from pathlib import Path
import unittest

import marge_runtime as runtime
from reframed.core.cbmodel import CBModel, CBReaction, Gene, GPRAssociation, Protein
from reframed.core.model import Compartment, Metabolite
from reframed.cobra.transcriptomics import marge


ROOT = Path(__file__).resolve().parents[1]


def toy(isozymes=False):
    model = CBModel('allocation_toy')
    model.add_compartment(Compartment('c'))
    for name in ['s', 'p']:
        model.add_metabolite(Metabolite(name, compartment='c'))
    for name in ['G_g1', 'G_g2']:
        model.add_gene(Gene(name))
    model.add_reaction(CBReaction('R_in', reversible=False, stoichiometry={'s': 1}, lb=0, ub=10))
    model.add_reaction(CBReaction('R_biomass_SC5_notrace', reversible=False, stoichiometry={'p': -1}, lb=0, ub=1000, objective=1))
    if isozymes:
        association = GPRAssociation()
        for name in ['G_g1', 'G_g2']:
            protein = Protein(); protein.genes = [name]; association.proteins.append(protein)
        model.add_reaction(CBReaction('R_iso', reversible=False, stoichiometry={'s': -1, 'p': 1}, lb=0, ub=1000, gpr_association=association))
    else:
        for index, name in enumerate(['G_g1', 'G_g2'], 1):
            association = GPRAssociation(); protein = Protein(); protein.genes = [name]; association.proteins.append(protein)
            model.add_reaction(CBReaction('R_path' + str(index), reversible=False, stoichiometry={'s': -1, 'p': 1}, lb=0, ub=1000, gpr_association=association))
    return model


class UpstreamMargeTests(unittest.TestCase):
    def setUp(self):
        protocol = json.loads((ROOT / 'research/ale_marge_protocol.json').read_text(encoding='utf-8'))
        runtime.activate(protocol)

    def test_analytical_two_path_solution_and_ratio_direction(self):
        a, b = marge(toy(), rel_expr={'G_g1': 2, 'G_g2': 0.5})
        self.assertAlmostEqual(a['R_path1'], 10 / 3)
        self.assertAlmostEqual(b['R_path1'], 20 / 3)
        self.assertAlmostEqual(a['R_path2'], 20 / 3)
        self.assertAlmostEqual(b['R_path2'], 10 / 3)
        joint = next(s for s in runtime.INSTANCES if s.model is None)
        self.assertAlmostEqual(joint.joint_snapshots[0]['objective_value'], 0)
        # The second stage should minimise total proteome, which on this toy model is u_g1 + u_g2 = 20.
        # ReFramed 1.2.1 files the reference-proteome capacity row under the same id as the proteome
        # pseudo-metabolite balance, so the balance is replaced and proteome_synth is decoupled from
        # gene usage: the second-stage objective collapses to its lower bound of zero.
        self.assertAlmostEqual(joint.joint_snapshots[1]['objective_value'], 0)

    def test_balance_preserving_variant_restores_the_second_stage(self):
        runtime.StudyGLPK.KEEP_PROTEOME_BALANCE = True
        try:
            marge(toy(), rel_expr={'G_g1': 2, 'G_g2': 0.5})
        finally:
            runtime.StudyGLPK.KEEP_PROTEOME_BALANCE = False
        joint = next(s for s in runtime.INSTANCES if s.model is None)
        # With the capacity row filed under its own id the balance survives, and the total gene usage
        # of the two paths equals the analytical value again instead of collapsing to zero.
        self.assertAlmostEqual(joint.joint_snapshots[1]['objective_value'], 20)

    def test_isozyme_gene_usage_changes_with_constant_net_reaction_flux(self):
        a, b = marge(toy(isozymes=True), rel_expr={'G_g1': 2, 'G_g2': 0.5})
        self.assertAlmostEqual(a['R_iso'], 10)
        self.assertAlmostEqual(b['R_iso'], 10)
        joint = next(s for s in runtime.INSTANCES if s.model is None)
        values = {v['id']: v['value'] for v in joint.joint_snapshots[1]['variables']}
        self.assertAlmostEqual(values['u_g1_b'] / values['u_g1_a'], 2)
        self.assertAlmostEqual(values['u_g2_b'] / values['u_g2_a'], 0.5)

    def test_zero_and_negative_relative_expression_are_rejected(self):
        for value in [0, -1]:
            with self.assertRaises(RuntimeError):
                marge(toy(), rel_expr={'G_g1': value})

    def test_infinite_bounds_are_converted_to_solver_unbounded_values(self):
        solver = runtime.StudyGLPK()
        solver.add_variable('x', 0, math.inf)
        self.assertIsNone(solver.problem.variables['x'].ub)
        self.assertEqual(solver.problem.configuration.tolerances.feasibility, 1e-9)

    def test_every_failed_solve_is_rejected(self):
        solver = runtime.StudyGLPK()
        solver.add_variable('x', 0, 1)
        solver.add_constraint('infeasible', {'x': 1}, '>', 2)
        with self.assertRaisesRegex(RuntimeError, 'MARGE solve failed'):
            solver.solve({'x': 1}, minimize=True)

    def test_full_paper_gate_stays_closed(self):
        from ale_marge_demo import strict_status
        protocol = json.loads((ROOT / 'research/ale_marge_protocol.json').read_text(encoding='utf-8'))
        state = strict_status(protocol)
        self.assertFalse(state['strict_entry_allowed'])
        self.assertFalse(state['paper_exact_reproduction_completed'])
        self.assertEqual(len(state['missing_required_inputs']), 5)


if __name__ == '__main__':
    unittest.main()
