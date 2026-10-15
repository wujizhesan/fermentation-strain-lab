import collections
import collections.abc
import hashlib
from importlib.metadata import version
import math
from pathlib import Path

import numpy
import scipy


if not hasattr(collections, 'Iterable'):
    collections.Iterable = collections.abc.Iterable
if not hasattr(scipy, 'compress'):
    scipy.compress = numpy.compress

import reframed
import reframed.solvers
from reframed.solvers.optlang_solver import OptLangSolver
from reframed.solvers.solver import Solver, VarType
from reframed.solvers.solution import Status
from optlang.glpk_interface import Model
import swiglpk as glp


INSTANCES = []


def finite_bound(value):
    return float(value) if value is not None and math.isfinite(value) else None


def export_lp(problem, solution):
    handle = problem.problem
    variables = []
    names = {}
    for index in range(1, glp.glp_get_num_cols(handle) + 1):
        name = glp.glp_get_col_name(handle, index)
        names[index] = name
        variable = problem.variables[name]
        variables.append({'id': name, 'lb': finite_bound(variable.lb), 'ub': finite_bound(variable.ub),
            'value': float(solution.values[name]), 'objective': glp.glp_get_obj_coef(handle, index)})
    rows = []
    length = len(variables) + 1
    indices, coefficients = glp.intArray(length), glp.doubleArray(length)
    for index in range(1, glp.glp_get_num_rows(handle) + 1):
        name = glp.glp_get_row_name(handle, index)
        constraint = problem.constraints[name]
        count = glp.glp_get_mat_row(handle, index, indices, coefficients)
        rows.append({'id': name, 'lb': finite_bound(constraint.lb), 'ub': finite_bound(constraint.ub),
            'terms': [[names[indices[j]], coefficients[j]] for j in range(1, count + 1)]})
    return {'status': solution.status.value, 'objective_value': float(solution.fobj),
        'direction': problem.objective.direction, 'variables': variables, 'constraints': rows}


class StudyGLPK(OptLangSolver):
    # Diagnostic switch only. Upstream marge() removes and rebuilds any row whose id repeats, so the
    # proteome pseudo-metabolite balance is destroyed by the later capacity row of the same name.
    # With the switch on, the second row is filed under proteome_cap_<condition> and the balance survives.
    KEEP_PROTEOME_BALANCE = False

    def __init__(self, model=None):
        Solver.__init__(self)
        self.model = model
        self.problem = Model()
        self.problem.configuration.tolerances.feasibility = 1e-9
        self.problem.configuration.verbosity = 0
        self.calls = []
        self.joint_snapshots = []
        if model is not None:
            self.build_problem(model)
        INSTANCES.append(self)

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, update=True):
        if (self.KEEP_PROTEOME_BALANCE and constr_id in ('proteome_a', 'proteome_b')
                and constr_id in self.constr_ids):
            constr_id = 'proteome_cap_' + constr_id.rsplit('_', 1)[-1]
        return super().add_constraint(constr_id, lhs, sense, rhs, update)

    def add_variable(self, var_id, lb=-math.inf, ub=math.inf, vartype=VarType.CONTINUOUS, update=True):
        return super().add_variable(var_id, finite_bound(lb), finite_bound(ub), vartype, update)

    def solve(self, *args, **kwargs):
        solution = super().solve(*args, **kwargs)
        if solution.status != Status.OPTIMAL:
            raise RuntimeError(f'MARGE solve failed: {solution.status.value}')
        self.calls.append({'status': solution.status.value, 'objective_value': float(solution.fobj),
            'biomass_a': solution.values.get('R_biomass_SC5_notrace_a') if solution.values else None,
            'biomass_b': solution.values.get('R_biomass_SC5_notrace_b') if solution.values else None,
            'reference_biomass': solution.values.get('R_biomass_SC5_notrace') if solution.values else None,
            'reference_proteome': solution.values.get('proteome_synth') if solution.values else None})
        if self.model is None and len(self.joint_snapshots) < 2:
            self.joint_snapshots.append(export_lp(self.problem, solution))
        return solution


def activate(protocol, keep_proteome_balance=False):
    root = Path(__file__).resolve().parent
    package = Path(reframed.__file__).parent
    if version('reframed') != protocol['reframed_version']:
        raise ValueError('Pinned ReFramed version required')
    for relative, expected in protocol['upstream_module_sha256'].items():
        if hashlib.sha256((package / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Upstream module changed: {relative}')
    INSTANCES.clear()
    StudyGLPK.KEEP_PROTEOME_BALANCE = keep_proteome_balance
    reframed.solvers.solvers['study_glpk'] = StudyGLPK
    reframed.solvers.set_default_solver('study_glpk')
    return {'reframed': version('reframed'), 'optlang': version('optlang'), 'scipy': scipy.__version__,
        'numpy': numpy.__version__, 'python_compatibility_aliases': ['collections.Iterable', 'scipy.compress'],
        'solver': 'GLPK', 'source_large_bounds_preserved': True,
        'proteome_balance_kept': keep_proteome_balance,
        'local_row_rename': 'proteome_cap_<condition>' if keep_proteome_balance else None}
