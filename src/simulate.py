import json
import sys
import numpy as np
import argparse as ap
import scipy.stats as stats
import simpy as sim
from pathlib import Path
from scipy.special import expit
from collections import defaultdict
from pandas import DataFrame, Series


def parse_arguments():
    p = ap.ArgumentParser()
    p.add_argument('--spec',
        default=Path('./spec.json'),
        type=Path,
        help='simulation specification file',
        metavar='PATH')
    p.add_argument('-T', '--time',
        help='number of time steps to simulate',
        type=int,
        default=10,
        metavar='INT')
    p.add_argument('-o', '--output',
        help='output experiment report to this file',
        type=Path,
        default=Path('./output.csv'),
        metavar='FILE')
    args = p.parse_args()
    return args


def parse_vector(v, length=None, dtype=np.float64):
    v = np.array(v, dtype=dtype)
    if (length is not None) and (len(v) != length):
        msg = f'vector had invalid length {len(v)} when {lenght} was expected'
        raise ValueError(msg)
    return v


def parse_matrix(m, shape=None, dtype=np.float64):
    m = np.array(m, dtype=dtype)
    if isinstance(shape, tuple) and (m.shape != shape):
        msg = f'data had invalid shape {m.shape} when {shape} was expected'
        raise ValueError(msg)
    return m


def validate_probabilities(v):
    if (v > 1).any() or (v < 0).any():
        msg = f'invalid probabilities (must be between 0 and 1)'
        raise ValueError(msg)


def validate_stochastic_vector(v):
    validate_probabilities(v)
    if v.sum() != 1.0:
        msg = f'stochastic vector expected instead of {v}'
        raise ValueError(msg)


def validate_stochastic_matrix(mat, rows=True, cols=False):
    (m, n) = mat.shape
    if rows:
        for i in range(m):
            validate_stochastic_vector(mat[i, :])
    if cols:
        for i in range(n):
            validate_stochastic_vector(mat[:, i])


def validate_positivity(v):
    if (v < 0).any():
        msg = f''
        ValueError(msg)


def parse_spec(args):
    with args.spec.open() as fp:
        spec = json.load(fp)
    for (k, v) in spec['species'].items():
        spec['species'][k] = parse_species(v) 
    return spec


def parse_species(spec):
    n = len(spec['states'])
    
    # parse and validate the transition matrix
    mat = parse_matrix(spec['transition'], shape=(n, n))
    validate_stochastic_matrix(mat)
    spec['transition'] = mat
    
    # parse and validate reproduction fitness
    v = spec['reproduction']['fitness']
    v = parse_vector(v, length=n)
    validate_probabilities(v)
    spec['reproduction']['fitness'] = v

    # parse and validate fertility data
    mu = spec['reproduction']['fertility']['mu']
    sigma = spec['reproduction']['fertility']['sigma']
    mu = parse_vector(mu, length=n)
    sigma = parse_vector(sigma, length=n)
    validate_positivity(mu)
    validate_positivity(sigma)
    spec['reproduction']['fertility']['mu'] = mu
    spec['reproduction']['fertility']['sigma'] = sigma

    # parse and validate initial population
    p = parse_vector(spec['population'], length=n, dtype=int)
    validate_positivity(p)
    spec['population'] = p

    # parse and validate capture probabilities
    p = parse_vector(spec['capture'], length=n)
    validate_probabilities(p)
    spec['capture'] = p

    # compute state index
    spec['state_idx'] = dict((spec['states'][i], i) for i in range(n))

    # parse state at birth
    b = spec['birth_state']
    if not (b in spec['state_idx']):
        msg = f'invalid state at birth: {b}'
        raise ValueError(msg)
    spec['birth_state'] = spec['state_idx'][b]
    
    return spec


class Population(object):
    def __init__(self, species, pop=False):
        self._species = species
        if pop:
            self._pop = species['population']
        else:
            self._pop = np.zeros(len(species['states']), dtype=int)

    def update_population(self, delta):
        self._pop += delta

    def transition(self):
        t = self._species['transition']
        n = len(self._pop)
        dists = [stats.multinomial(self._pop[i], t[i, :]) for i in range(n)]
        pop = np.zeros(n, dtype=int)
        for d in dists:
            pop += d.rvs().reshape((3,))
        self._pop = pop
        return pop

    def capture(self):
        ns = self._pop
        ps = self._species['capture']
        cap = stats.binom(ns, ps).rvs()
        return cap

    def reproduce(self):
        p = self._species['reproduction']['fitness']
        n = stats.binom(self._pop, p).rvs()
        mu = self._species['reproduction']['fertility']['mu']
        sigma = self._species['reproduction']['fertility']['sigma']
        rep = stats.norm(n * mu, np.sqrt(n) * sigma).rvs()
        rep = np.round(np.clip(rep, 0.0, None)).astype(int)
        rep = rep.sum()
        spawn = np.zeros(len(self._species['states']), dtype=int)
        spawn[self._species['birth_state']] = rep
        return spawn

    def population(self):
        return self._pop.copy()


def zero_factory(n):
    return lambda: np.zeros(n, dtype=int)


class Experiment(object):
    _cohort_names = ('observed', 'marked', 'unknown')
    _report_cols = (
        'time',
        'species',
        'cohort',
        'state',
        'capture',
        'population'
    )

    def __init__(self, env, spec):
        self._env = env
        self._spec = spec
        self._make_cohorts() 
        self._make_empty_report()
        self._process = env.process(self.perform())

    def _make_empty_report(self):
        self._report = DataFrame(columns=self._report_cols)

    def _make_cohorts(self):
        cohorts = dict()
        for (s, species) in self._spec['species'].items():
            cohorts[s] = dict()
            for c in self._cohort_names:
                pop = Population(species, pop=(c == 'unknown'))
                cohorts[s][c] = pop
        self._cohorts = cohorts

    def _log_data(self, species, cohort, capture):
        states = self._spec['species'][species]['states']
        n = len(states) 
        now = self._env.now
        entry = {
            'time': [now] * n,
            'species': [species] * n,
            'cohort': [cohort] * n,
            'state': states,
            'capture': capture,
            'population': self._cohorts[species][cohort].population()
        }
        entry = DataFrame(data=entry)
        self._report = self._report.append(entry, ignore_index=True)

    def perform(self):
        while True:
            now = self._env.now
            print(f'time = {now}', file=sys.stderr)
            for s in self._spec['species'].keys():
                n = len(self._spec['species'][s]['states'])

                # capture fish and record data
                cap = defaultdict(zero_factory(n))
                for c in self._cohort_names:
                    tmp = self._cohorts[s][c].capture()
                    self._log_data(s, c, tmp)
                    cap[c] -= tmp
                    cap['observed'] += tmp
                
                # update cohort populations to reflect the result
                # of capture expedition
                for c in self._cohort_names:
                    self._cohorts[s][c].update_population(cap[c])

                # allow fish to reproduce
                rep = defaultdict(zero_factory(n))
                for c in self._cohort_names:
                    rep[c] = self._cohorts[s][c].reproduce()
                self._cohorts[s]['observed'].update_population(0)
                self._cohorts[s]['marked'].update_population(rep['observed'])
                self._cohorts[s]['unknown'].update_population(rep['marked'] + rep['unknown'])
            
            # advance to next time step
            yield self._env.timeout(1.0)
            for s in self._spec['species'].keys():
                for c in self._cohort_names:
                    self._cohorts[s][c].transition()

    def report(self):
        return self._report.copy()

def main():
    args = parse_arguments()
    spec = parse_spec(args)
    env = sim.Environment()
    experiment = Experiment(env, spec)
    env.run(until=args.time)
    experiment.report().to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
