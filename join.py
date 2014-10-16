#!/usr/bin/env python

# Copyright 2012 Esteban Hurtado
#
# This file is part of Cutedots.
#
# Cutedots is distributed under the terms of the Reciprocal Public License 1.5.
#
# Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.
#
# You should have received a copy of the Reciprocal Public License along with
# Cutedots. If not, see <http://http://opensource.org/licenses/rpl-1.5>.

import pickle
from raw_curves import point_sq_dist, point_ave

# Max. sq. deviation
max_sq = 300


###########################################################################################

# Exception for breaking nested for
class break_nested_for(Exception):
    pass

# Parse command line
import argparse
parser = argparse.ArgumentParser(description='Display 3D motion data.')
parser.add_argument('file', type=argparse.FileType('rb'), help="Input data file")
parser.add_argument('-g', type=int, default=10, help='max number of frames to fill')
args = parser.parse_args()
max_gap = args.g

# Open data file
data = pickle.load(args.file)

print("Initial number of curves:", len(data.curves))

# Find pairs of curves that would fit
#####################################

try:

    # For frame gaps from 1 to max_gap
    for gap in range(1, max_gap):
        idx1 = (gap-1)
        idx2 = gap-1-idx1
        changed = True
        while changed:
            changed = False
            N = len(data.curves)
            try:
                # For each pair of curves
                for i in range(N):
                    for j in range(N):
                        if i != j:
                            c1 = data.curves[i]
                            c2 = data.curves[j]

                            # Check gap between curves
                            if (c1.end_frame + gap) == c2.begin_frame:

                                # for each point in gap see if predictions are close
                                for g in range(gap):
                                    # Predict missing point from both curves
                                    p1 = c1.predict(g)
                                    p2 = c2.back_predict(gap-1-g)
                                    sq_dev = point_sq_dist(p1,p2)

                                    # if predictions are close
                                    if sq_dev < max_sq:
                                        print("joining curves", c1.name, "and", c2.name)
                                        # Append missing points to first curve
                                        for i in range(gap):
                                            p1 = c1.predict(i)
                                            p2 = c2.back_predict(gap-1-i)
                                            missing_p = point_ave(p1,p2)
                                            c1.point_data.append(missing_p)
                                        # Then append second curve points to first curve
                                        c1.point_data.extend(c2.point_data)
                                        c1.end_frame = c1.begin_frame + len(c1.point_data)
                                        # Remove second curve
                                        data.curves.remove(c2)
                                        changed = True
                                        raise break_nested_for
            except break_nested_for:
                pass

except:
    print("There were errors")
    pass

print("Final number of curves:", len(data.curves))
print("Writing output...", end=' ')
pickle.dump(data, file(args.file.name, "wb"), 2)
print("[Ok]")
