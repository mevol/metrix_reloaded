from __future__ import absolute_import, division, print_function
import procrunner
import os, sys
import argparse
import json


def apply_dmin(mtz_file, dmin):

    prefix, ext = os.path.splitext(mtz_file)
    hklout = prefix + "_" + str(dmin) + "A" + ext
    keywords = "RESOLUTION {0}\nGO".format(dmin)
    keywords = keywords.encode()
    try:
        result = procrunner.run(["mtzutils", "hklin", mtz_file,
            "hklout", hklout], stdin=keywords)
    except FileNotFoundError:
        sys.exit("Could not run mtzutils")

    if result.returncode or result.stderr:
        sys.exit("mtzutils failed")

    return hklout

def merging_statistics(mtz_file):

    prefix, ext = os.path.splitext(mtz_file)
    stats = prefix + "_stats.json"
    try:
        result = procrunner.run(["xia2.merging_statistics", mtz_file,
            "anomalous=true", "json.file_name={}".format(stats), "json.indent=2"])
    except FileNotFoundError:
        sys.exit("Could not run xia2.merging_statistics")

    if result.returncode or result.stderr:
        sys.exit("xia2.merging_statistics failed")

    return stats

if __name__ == "__main__":

    description = "Calculate data processing stats at a chosen resolution cutoff"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("hklin", help="Scaled, unmerged MTZ file")
    parser.add_argument("dmin", type=float, help="High resolution limit")
    options = parser.parse_args()

    hklin = options.hklin
    dmin = options.dmin

    if not os.path.exists(options.hklin):
        sys.exit("Cannot find MTZ file {}".format(options.hklin))

    new_mtz = apply_dmin(hklin, dmin)
    stats = merging_statistics(new_mtz)

    with open(stats, "r") as f:
        stats = json.load(f)

    # Example of extracting overall statistics
    overall = stats['overall']
    print("Anomalous slope: {0:.5f}".format(overall['anom_probability_plot_expected_delta']['slope']))
    print("Anomalous completeness: {0:.5f}".format(overall['anom_completeness']))
    print("dF/F: {0:.5f}".format(overall['anom_signal']))
    print("dI/s(dI): {0:.5f}".format(overall['delta_i_mean_over_sig_delta_i_mean']))
