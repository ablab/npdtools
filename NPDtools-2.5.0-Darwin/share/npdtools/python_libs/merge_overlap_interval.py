#!/usr/bin/env python

import sys


def __read_intervals(fpath):
    intervals = []
    with open(fpath) as f:
        for line in f:
            start, end = map(int, line.split())
            intervals.append((start, end))
    return intervals


def __write_intervals(handler, intervals):
    for (start, end) in intervals:
        handler.write('{start} {end}\n'.format(**locals()))


def __merge_intervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


# for replacing merge_overlap_interval.cpp
def main(fpath):
    intervals = __read_intervals(fpath)
    merged = __merge_intervals(intervals)
    __write_intervals(sys.stdout, merged)


if __name__ == '__main__':
    main(sys.argv[1])