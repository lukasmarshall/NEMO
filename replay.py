# Copyright (C) 2012, 2013, 2014 Ben Elliston
# Copyright (C) 2014, 2015 The University of New South Wales
#
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

"""Replay runs from a text file of generators."""
import argparse
import json
import re
import numpy as np

import costs
import configfile as cf
import nem
import scenarios
import utils

parser = argparse.ArgumentParser(description='Bug reports to: nemo-devel@lists.ozlabs.org')
parser.add_argument("-f", type=str, help='replay file', required=True)
parser.add_argument("-d", "--demand-modifier", type=str, default=[], action="append", help='demand modifier')
parser.add_argument("--no-legend", action="store_false", help="hide legend")
parser.add_argument("-t", "--transmission", action="store_true", help="show region exchanges [default: False]")
parser.add_argument("-v", action="count", help='verbose mode')
parser.add_argument("-x", action="store_true", help='producing a balancing plot')
parser.add_argument("--nsp-limit", type=float, default=cf.get('limits', 'nonsync-penetration'),
                    help='Non-synchronous penetration limit [default: %s]' %
                    cf.get('limits', 'nonsync-penetration'))
parser.add_argument("--spills", action="store_true", help='plot spills')
args = parser.parse_args()

context = nem.Context()
context.costs = costs.NullCosts()
assert 0 <= args.nsp_limit <= 1
context.nsp_limit = args.nsp_limit
# Apply each demand modifier argument (if any) in the given order.
for arg in args.demand_modifier:
    scenarios.demand_switch(arg)(context)


def run_one(chromosome):
    """Run a single simulation."""
    context.set_capacities(chromosome)
    context.track_exchanges = args.transmission
    context.verbose = args.v > 1
    nem.run(context)
    context.verbose = args.v > 0
    print context
    if args.transmission:
        np.set_printoptions(precision=3)
        x = context.exchanges.max(axis=0)
        print np.array_str(x, precision=1, suppress_small=True)
        with open('results.json', 'w') as f:
            json.dump(x.tolist(), f)
    print


with open(args.f) as replayfile:
    for line in replayfile:
        if re.search(r'^\s*$', line):
            continue
        if re.search(r'^\s*#', line):
            print line,
            continue
        if not re.search(r'^\s*[\w\-\+]+:\s*\[.*\]\s*.?$', line):
            print 'skipping malformed input:', line
            continue
        match = re.match(r'^\s*([\w\-\+]+):\s*\[(.*)\]\s*.?$', line)
        scenario = match.group(1)
        capacitylist = match.group(2)
        scenarios.supply_switch(scenario)(context)
        capacities = [float(st) for st in capacitylist.split(',')]  # str -> float
        print 'scenario', scenario
        run_one(capacities)
        if args.x:  # pragma: no cover
            utils.plot(context, spills=args.spills, showlegend=args.no_legend)
