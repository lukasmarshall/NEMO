# Copyright (C) 2012, 2013, 2014 Ben Elliston
# Copyright (C) 2014, 2015, 2016 The University of New South Wales
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

"""Evolutionary programming applied to NEM optimisations."""

import sys
import csv
import json
import argparse
import numpy as np

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from deap import cma
try:
    from scoop import futures
except ImportError:  # pragma: no cover
    print 'WARNING: scoop not loaded'

import nem
import generators
import scenarios
import costs
import configfile as cf
import consts
import transmission

parser = argparse.ArgumentParser(description='Bug reports to: nemo-devel@lists.ozlabs.org')
parser.add_argument("-c", "--carbon-price", type=int,
                    default=cf.get('costs', 'co2-price-per-t'),
                    help='carbon price ($/t) [default: %s]' % cf.get('costs', 'co2-price-per-t'))
parser.add_argument("-d", "--demand-modifier", type=str, action="append", help='demand modifier [default: unchanged]')
parser.add_argument("-g", "--generations", type=int, default=cf.get('optimiser', 'generations'),
                    help='generations [default: %s]' % cf.get('optimiser', 'generations'))
parser.add_argument("-r", "--discount-rate", type=float, default=cf.get('costs', 'discount-rate'),
                    help='discount rate [default: %s]' % cf.get('costs', 'discount-rate'))
parser.add_argument("-s", "--supply-scenario", type=str, default='re100', help='generation mix scenario [default: \'re100\']')
parser.add_argument("-t", "--transmission", action="store_true", help="include transmission [default: False]")
parser.add_argument("-v", "--verbose", action="store_true", help="be verbose")

parser.add_argument("--bioenergy-limit", type=float, default=cf.get('limits', 'bioenergy-twh-per-yr'),
                    help='Limit on annual energy from bioenergy (TWh/y) [default: %s]' %
                    cf.get('limits', 'bioenergy-twh-per-yr'))
parser.add_argument("--ccs-storage-costs", type=float, default=cf.get('costs', 'ccs-storage-costs-per-t'),
                    help='CCS storage costs ($/t) [default: %s]' % cf.get('costs', 'ccs-storage-costs-per-t'))
parser.add_argument("--coal-price", type=float, default=cf.get('costs', 'coal-price-per-gj'),
                    help='black coal price ($/GJ) [default: %s]' % cf.get('costs', 'coal-price-per-gj'))
parser.add_argument("--costs", type=str, default=cf.get('costs', 'technology-cost-class'),
                    help='cost scenario [default: %s]' % cf.get('costs', 'technology-cost-class'))
parser.add_argument("--emissions-limit", type=float, help='CO2 emissions limit (Mt/y) [default: None]')
parser.add_argument("--fossil-limit", type=float, help='Fraction of energy from fossil fuel [default: None]')
parser.add_argument("--gas-price", type=float, default=cf.get('costs', 'gas-price-per-gj'),
                    help='gas price ($/GJ) [default: %s]' % cf.get('costs', 'gas-price-per-gj'))
parser.add_argument("--hydro-limit", type=float, default=cf.get('limits', 'hydro-twh-per-yr'),
                    help='Limit on annual energy from hydro (TWh/y) [default: %s]' %
                    cf.get('limits', 'hydro-twh-per-yr'))
parser.add_argument("--lambda", type=int, dest='lambda_', help='override CMA-ES lambda value')
parser.add_argument("--list-scenarios", action="store_true")
parser.add_argument("--min-regional-generation", type=float,
                    help='minimum share of energy generated intra-region [default: None]')
parser.add_argument("--nsp-limit", type=float, default=cf.get('limits', 'nonsync-penetration'),
                    help='Non-synchronous penetration limit [default: %s]' %
                    cf.get('limits', 'nonsync-penetration'))
parser.add_argument("--reliability-std", type=float, help='reliability standard (%% unserved)')
parser.add_argument("--reserves", type=int, default=cf.get('limits', 'minimum-reserves-mw'),
                    help='minimum operating reserves [default: %s MW]' %
                    cf.get('limits', 'minimum-reserves-mw'))
parser.add_argument("--seed", type=int, help='seed for random number generator [default: None]')
parser.add_argument("--sigma", type=float, default=cf.get('optimiser', 'sigma'),
                    help='CMA-ES sigma value [default: %s]' % cf.get('optimiser', 'sigma'))
parser.add_argument("--trace-file", type=str, help='Filename for evaluation trace (comma separated) [default: None]')
parser.add_argument('--version', action='version', version='1.0')
args = parser.parse_args()

if __name__ == '__main__' and args.list_scenarios:
    for key in sorted(scenarios.supply_scenarios):
        descr = scenarios.supply_scenarios[key].__doc__
        print '%20s' % key, '\t', descr.split('\n')[0]
    print
    sys.exit(0)

if __name__ == '__main__':
    print vars(args)

np.set_printoptions(precision=5)
context = nem.Context()

# Set the system non-synchronous penetration limit.
context.nsp_limit = args.nsp_limit
assert 0 <= context.nsp_limit <= 1, \
    "NSP limit must be in the interval [0,1]"

# Override the reliability standard (if the user gives this option).
if args.reliability_std is not None:
    context.relstd = args.reliability_std

# Likewise for the minimum share of regional generation.
context.min_regional_generation = args.min_regional_generation
assert context.min_regional_generation is None or \
    (0 <= context.min_regional_generation <= 1), \
    "Minimum regional generation must be in the interval [0,1]"

cost_class = costs.cost_switch(args.costs)
context.costs = cost_class(args.discount_rate, args.coal_price, args.gas_price, args.ccs_storage_costs)
context.costs.carbon = args.carbon_price


def txcost(x):
    """Transmission cost expression."""
    return 0 if x == 0 else 965 if x > 5000 else 16319 * pow(x, -0.332)

context.costs.transmission = transmission.Transmission(txcost, args.discount_rate)
context.track_exchanges = args.transmission

# Set up the scenario.
scenarios.supply_switch(args.supply_scenario)(context)
# Apply each demand modifier in the order given on the command line.
if args.demand_modifier is not None:
    for arg in args.demand_modifier:
        scenarios.demand_switch(arg)(context)

if args.verbose and __name__ == '__main__':
    docstring = scenarios.supply_switch(args.supply_scenario).__doc__
    assert docstring is not None
    # Prune off any doctest test from the docstring.
    docstring = docstring.split('\n')[0]
    print "supply scenario: %s (%s)" % (args.supply_scenario, docstring)
    print context.generators

if args.trace_file is not None:
    with open(args.trace_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['# score', 'penalty', 'reasoncode', 'parameter values'])
        csvfile.close()

reasons = {'unserved': 1,
           'emissions': 2,
           'fossil': 4,
           'bioenergy': 8,
           'hydro': 16,
           'reserves': 32,
           'min-regional-gen': 64}


def cost(ctx):
    """Sum up the costs."""
    score = 0
    for g in ctx.generators:
        score += (g.capcost(ctx.costs) / ctx.costs.annuityf * ctx.years) \
            + g.opcost(ctx.costs)

    penalty = 0
    reason = 0

    ### Penalty: unserved energy
    minuse = ctx.demand.sum() * (ctx.relstd / 100)
    use = max(0, ctx.unserved_energy - minuse)
    if use > 0:
        reason |= reasons['unserved']
    penalty += pow(use, 3)

    ### Penalty: minimum reserves
    if args.reserves > 0:
        for i in range(ctx.timesteps):
            reserve, spilled = 0, 0
            for g in context.generators:
                try:
                    spilled += g.hourly_spilled[i]
                except KeyError:
                    # non-variable generators may not have spill data
                    pass

                # Calculate headroom for each generator, except pumped hydro and
                # CST -- tricky to calculate capacity
                if isinstance(g, nem.generators.Fuelled) and not \
                   isinstance(g, nem.generators.PumpedHydro) and not \
                   isinstance(g, nem.generators.CST):
                    reserve += g.capacity - g.hourly_power[i]

            if reserve + spilled < args.reserves:
                reason |= reasons['reserves']
                penalty += pow(args.reserves - reserve + spilled, 3)

    ### Penalty: minimum share of regional generation
    if ctx.min_regional_generation is not None:
        regional_generation_shortfall = 0
        for rgn in ctx.regions:
            regional_generation = 0
            for g in ctx.generators:
                if g.region() is rgn:
                    regional_generation += sum(g.hourly_power.values())
            min_regional_generation = sum(context.demand[rgn]) * ctx.min_regional_generation
            regional_generation_shortfall += max(0, min_regional_generation - regional_generation)
        if regional_generation_shortfall > 0:
            reason |= reasons['min-regional-gen']
        penalty += pow(regional_generation_shortfall, 3)

    ### Penalty: total emissions
    if args.emissions_limit is not None:
        emissions = 0
        for g in ctx.generators:
            try:
                emissions += sum(g.hourly_power.values()) * g.intensity
            except AttributeError:
                # not all generators have an intensity attribute
                pass
        # exceedance in tonnes CO2-e
        emissions_exceedance = max(0, emissions - args.emissions_limit * pow(10, 6) * ctx.years)
        if emissions_exceedance > 0:
            reason |= reasons['emissions']
        penalty += pow(emissions_exceedance, 3)

    ### Penalty: limit fossil to fraction of annual demand
    if args.fossil_limit is not None:
        fossil_energy = 0
        for g in ctx.generators:
            if isinstance(g, generators.Fossil):
                fossil_energy += sum(g.hourly_power.values())
        fossil_exceedance = max(0, fossil_energy - ctx.demand.sum() * args.fossil_limit * ctx.years)
        if fossil_exceedance > 0:
            reason |= reasons['fossil']
        penalty += pow(fossil_exceedance, 3)

    ### Penalty: limit biofuel use
    biofuel_energy = 0
    for g in ctx.generators:
        if isinstance(g, generators.Biofuel):
            biofuel_energy += sum(g.hourly_power.values())
    biofuel_exceedance = max(0, biofuel_energy - args.bioenergy_limit * consts.twh * ctx.years)
    if biofuel_exceedance > 0:
        reason |= reasons['bioenergy']
    penalty += pow(biofuel_exceedance, 3)

    ### Penalty: limit hydro use
    hydro_energy = 0
    for g in ctx.generators:
        if isinstance(g, generators.Hydro) and \
           not isinstance(g, generators.PumpedHydro):
            hydro_energy += sum(g.hourly_power.values())
    hydro_exceedance = max(0, hydro_energy - args.hydro_limit * consts.twh * ctx.years)
    if hydro_exceedance > 0:
        reason |= reasons['hydro']
    penalty += pow(hydro_exceedance, 3)

    if args.transmission:
        maxexchanges = ctx.exchanges.max(axis=0)
        np.fill_diagonal(maxexchanges, 0)
        for i in range(1, maxexchanges.shape[0]):
            # then put the max (upper, lower) into lower
            # and zero the upper entries
            for j in range(i):
                maxexchanges[i, j] = max(maxexchanges[i, j], maxexchanges[j, i])
                maxexchanges[j, i] = 0
        try:
            # existing transmission is "free"
            maxexchanges = np.maximum(0, maxexchanges - nem.polyons.existing_net)
        except AttributeError:
            # skip if not present
            pass

        costmat = ctx.costs.transmission.cost_matrix(maxexchanges) * ctx.years
        score += costmat.sum()

    # Express $/yr as an average $/MWh over the period
    return score / ctx.demand.sum(), penalty / ctx.demand.sum(), reason


def set_generators(chromosome):
    """Set the generator list from the chromosome."""
    i = 0
    for gen in context.generators:
        for (setter, min_cap, max_cap) in gen.setters:
            # keep parameters within bounds
            newval = max(min(chromosome[i], max_cap), min_cap)
            setter(newval)
            i += 1
    # Check every parameter has been set.
    assert i == len(chromosome), '%d != %d' % (i, len(chromosome))


def eval_func(chromosome):
    """Average cost of energy (in $/MWh)."""
    set_generators(chromosome)
    nem.run(context)
    score, penalty, reason = cost(context)
    if args.trace_file is not None:
        # write the score and individual to the trace file
        with open(args.trace_file, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([score, penalty, reason] + list(chromosome))
    return score + penalty,


creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)
toolbox = base.Toolbox()
try:
    toolbox.register("map", futures.map)
except NameError:  # pragma: no cover
    pass

# See:
# https://deap.readthedocs.org/en/master/api/algo.html#deap.cma.Strategy
# for additional parameters that can be passed to cma.Strategy.
numparams = sum([len(g.setters) for g in context.generators])
if args.lambda_ is None:
    # let DEAP choose
    strategy = cma.Strategy(centroid=[0] * numparams, sigma=args.sigma)
else:
    strategy = cma.Strategy(centroid=[0] * numparams, sigma=args.sigma, lambda_=args.lambda_)

toolbox.register("generate", strategy.generate, creator.Individual)
toolbox.register("update", strategy.update)
toolbox.register("evaluate", eval_func)


def run():
    """Run the evolution."""
    if args.verbose and __name__ == '__main__':
        print "objective: minimise", eval_func.__doc__

    if args.seed is not None:
        np.random.seed(args.seed)
    hof = tools.HallOfFame(1)
    stats_fit = tools.Statistics(lambda ind: ind.fitness.values)
    stats_hof = tools.Statistics(lambda ignored: hof[0].fitness.values)
    mstats = tools.MultiStatistics(fitness=stats_fit, hallfame=stats_hof)
    mstats.register("min", np.min)

    try:
        algorithms.eaGenerateUpdate(toolbox, ngen=args.generations,
                                    stats=mstats, halloffame=hof, verbose=True)
    except KeyboardInterrupt:  # pragma: no cover
        print 'user terminated early'

    (score,) = hof[0].fitness.values
    print 'List:', [max(0, param) for param in hof[0]]

    set_generators(hof[0])
    nem.run(context)
    context.verbose = True
    print context
    score, penalty, reason = cost(context)
    print 'Score: %.2f $/MWh' % score
    if reason > 0:
        print 'Penalty: %.2f $/MWh' % penalty
        print 'Constraints violated:',
        for label, code in zip(reasons, reasons.values()):
            if reason & code:
                print label,
        print
    if args.transmission:
        x = context.exchanges.max(axis=0)
        print np.array_str(x, precision=1, suppress_small=True)
        f = open('results.json', 'w')
        obj = {'exchanges': x.tolist(), 'generators': context}
        json.dump(obj, f, cls=nem.Context.JSONEncoder, indent=True)
        f.close()
    print 'Done'

if __name__ == '__main__':
    run()
