#!/usr/bin/env python
"""Calculate and save jet flavor fractions from the Coffea producer outputs."""

import argparse
import json
from pathlib import Path

import numpy as np
from coffea import util
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

import helpers as h
from common_binning import JERC_Constants
from JetEtaBins import JetEtaBins, PtBins
from fileNames.available_datasets import dataset_dictionary
from plotters.pltStyle import pltStyle
from uncertainty_helpers import (
    FlavorFractions,
    combine_flavors,
    get_spline_smooth,
    sum_output,
)


SCHEMA_VERSION = 1
DEFAULT_HERWIG_SAMPLES = ("QCD-MG-Her", "DY-MG-Her", "Herwig-TTBAR")
DEFAULT_PYTHIA_SAMPLES = ("QCD-MG-Py", "DY-MG-Py", "Pythia-TTBAR")


def parse_sample_list(value):
    return tuple(sample.strip() for sample in value.split(",") if sample.strip())


def _sample_paths(samples, input_dir):
    return {
        sample: Path(input_dir) / f"CoffeaJERCOutputs_L5_{sample}.coffea"
        for sample in samples
    }


def fraction_metadata(
    herwig_samples,
    pythia_samples,
    eta_binning,
    include_unmatched,
    combine_antiflavour,
    input_dir="out",
):
    samples = tuple(herwig_samples) + tuple(pythia_samples)
    unknown_samples = [
        sample
        for sample in dict.fromkeys(samples)
        if not any(
            sample == known or sample.startswith(known + "_")
            for known in dataset_dictionary
        )
    ]
    if unknown_samples:
        raise ValueError(
            "Samples are not registered in dataset_dictionary: "
            + ", ".join(unknown_samples)
        )
    paths = _sample_paths(samples, input_dir)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Coffea inputs:\n" + "\n".join(missing))

    input_files = {}
    dataset_config = {}
    for sample, path in paths.items():
        stat = path.stat()
        input_files[sample] = {
            "path": str(path.resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
        registry_key = next(
            known
            for known in dataset_dictionary
            if sample == known or sample.startswith(known + "_")
        )
        _, xsec_config, _ = dataset_dictionary[registry_key]
        dataset_config[sample] = {
            "registry_key": registry_key,
            "xsec_config": xsec_config,
        }
        if isinstance(xsec_config, str):
            xsec_path = Path(__file__).resolve().parent / xsec_config
            if xsec_path.is_file():
                xsec_stat = xsec_path.stat()
                dataset_config[sample]["xsec_file"] = {
                    "path": str(xsec_path.resolve()),
                    "size": xsec_stat.st_size,
                    "mtime_ns": xsec_stat.st_mtime_ns,
                }

    return {
        "schema_version": SCHEMA_VERSION,
        "herwig_samples": list(herwig_samples),
        "pythia_samples": list(pythia_samples),
        "eta_binning": eta_binning,
        "pt_binning": "MC_truth",
        "include_unmatched": bool(include_unmatched),
        "combine_antiflavour": bool(combine_antiflavour),
        "input_files": input_files,
        "dataset_config": dataset_config,
    }


def default_cache_path(metadata, output_dir="out_txt"):
    sample_tags = "-".join(
        metadata["herwig_samples"] + metadata["pythia_samples"]
    )
    return Path(output_dir) / (
        f"flavor_fractions_{sample_tags}_{metadata['eta_binning']}.json"
    )


def _fraction_flavors(include_unmatched, combine_antiflavour):
    flavors = (
        ["g", "cs", "b", "ud", "unmatched"]
        if include_unmatched
        else ["g", "c", "b", "ud", "s"]
    )
    if not combine_antiflavour:
        flavors = h.get_flavor_antiflavor_list(flavors)
        flavors.extend(
            flavor
            for flavor in _original_flavors(include_unmatched)
            if flavor not in h.barable_flavors
        )
    return flavors


def _original_flavors(include_unmatched):
    flavors = (
        ["g", "cs", "b", "ud", "unmatched"]
        if include_unmatched
        else ["g", "c", "b", "ud", "s"]
    )
    return [flav for flav in flavors if flav != "g"]


def _stored_flavors(include_unmatched, combine_antiflavour):
    if combine_antiflavour:
        return _fraction_flavors(include_unmatched, True)
    flavors = []
    for flavor in _original_flavors(include_unmatched):
        flavors.append(flavor)
        if flavor in h.barable_flavors:
            flavors.append(flavor + "bar")
    return flavors


def _make_fraction_objects(sample_data, pt_bins, jeteta_bins, eta_binning):
    qfracs = {
        flav: np.asarray(values, dtype=float)
        for flav, values in sample_data["mean"].items()
    }
    qfrac_var = {
        flav: np.asarray(values, dtype=float)
        for flav, values in sample_data["variance"].items()
    }
    for flavor, values in qfracs.items():
        for eta_idx in range(jeteta_bins.nbins):
            valid = np.isfinite(values[:, eta_idx]) & (values[:, eta_idx] != 0)
            if np.count_nonzero(valid) < 4:
                raise ValueError(
                    f"Insufficient nonzero fraction bins to build a spline for "
                    f"{flavor} in eta bin {eta_idx}."
                )
    qfrac_splines = {
        flav: np.array(
            [
                get_spline_smooth(
                    qfracs[flav][:, eta_idx],
                    qfracs[flav][:, eta_idx],
                    pt_bins,
                    s=0,
                )
                for eta_idx in range(jeteta_bins.nbins)
            ]
        )
        for flav in qfracs
    }
    qfrac_splines_2d = {
        flav: RegularGridInterpolator(
            (np.log10(pt_bins.centres), jeteta_bins.centres),
            qfracs[flav],
            fill_value=None,
        )
        for flav in qfracs
    }
    return (
        qfracs,
        qfrac_var,
        qfrac_splines,
        qfrac_splines_2d,
        FlavorFractions(qfrac_splines, eta_binning),
    )


def load_fraction_cache(cache_path, expected_metadata):
    """Load a fraction cache only when its setup and input files match."""
    cache_path = Path(cache_path)
    with cache_path.open() as cache_file:
        payload = json.load(cache_file)

    if not isinstance(payload, dict):
        raise ValueError("Fraction cache must contain a JSON object.")
    if payload.get("metadata") != expected_metadata:
        raise ValueError("Fraction cache metadata does not match this setup.")
    if "samples" not in payload or "binning" not in payload:
        raise ValueError("Fraction cache is missing required data.")
    if not isinstance(payload["binning"], dict) or not isinstance(
        payload["samples"], dict
    ):
        raise ValueError("Fraction cache has invalid binning or sample data.")

    eta_binning = expected_metadata["eta_binning"]
    jeteta_bins = JetEtaBins(eta_binning, absolute=True)
    pt_bins = PtBins("MC_truth")
    binning = payload["binning"]
    if not np.array_equal(binning.get("pt_edges"), pt_bins.edges):
        raise ValueError("Fraction cache has incompatible pt binning.")
    if not np.array_equal(binning.get("eta_edges"), jeteta_bins.edges):
        raise ValueError("Fraction cache has incompatible eta binning.")
    expected_samples = set(
        expected_metadata["herwig_samples"] + expected_metadata["pythia_samples"]
    )
    if set(payload["samples"]) != expected_samples:
        raise ValueError("Fraction cache does not contain the expected samples.")
    expected_flavors = set(
        _stored_flavors(
            expected_metadata["include_unmatched"],
            expected_metadata["combine_antiflavour"],
        )
    )
    sample_plot = {}
    qfrac_spline_dict = {}
    for sample, sample_data in payload["samples"].items():
        if not isinstance(sample_data, dict):
            raise ValueError(f"Fraction cache data for {sample} is invalid.")
        if set(sample_data.get("mean", {})) != expected_flavors or set(
            sample_data.get("variance", {})
        ) != expected_flavors:
            raise ValueError(
                f"Fraction cache flavors for {sample} do not match this setup."
            )
        for field in ("mean", "variance"):
            for flavor, values in sample_data[field].items():
                if np.asarray(values).shape != (pt_bins.nbins, jeteta_bins.nbins):
                    raise ValueError(
                        f"Fraction cache {field} for {sample}/{flavor} has "
                        "incompatible dimensions."
                    )
        qfracs, qfrac_var, splines, splines_2d, fraction_obj = (
            _make_fraction_objects(sample_data, pt_bins, jeteta_bins, eta_binning)
        )
        sample_plot[sample] = [qfracs, qfrac_var, splines, splines_2d]
        qfrac_spline_dict[sample] = fraction_obj

    return {
        "samples": payload["samples"],
        "sample_plot": sample_plot,
        "qfrac_spline_dict": qfrac_spline_dict,
        "pt_bins": pt_bins,
        "jeteta_bins": jeteta_bins,
        "labels": payload.get("labels", {}),
    }


def _calculate_sample_fractions(
    sample,
    input_path,
    flavors,
    include_unmatched,
    combine_antiflavour,
    jeteta_bins,
    pt_bins,
):
    output = util.load(str(input_path))
    xsec_dict, label = h.get_xsec_dict(sample, dataset_dictionary)
    output = sum_output(output, sample, xsec_dict)
    flavors_to_obtain = flavors + ["all"] if include_unmatched else flavors
    histograms = combine_flavors(
        output,
        flavors_to_obtain,
        sumeta=False,
        combine_antiflavour=combine_antiflavour,
    )
    rebinned = {
        flav: h.rebin_hist(
            h.sum_neg_pos_eta(histogram), "jeteta", jeteta_bins.edges
        )
        for flav, histogram in histograms.items()
    }

    if not include_unmatched:
        rebinned["all"] = sum(rebinned[flav] for flav in flavors)

    values = {}
    for flav, histogram in rebinned.items():
        vals = histogram.values().copy()
        vals[vals == 0] = np.nan
        values[flav] = vals

    if combine_antiflavour:
        fractions = {flav: values[flav] / values["all"] for flav in flavors}
        total_variance = rebinned["all"].variances() / values["all"] ** 2
        fraction_variances = {
            flav: fractions[flav] ** 2
            * (
                rebinned[flav].variances() / values[flav] ** 2
                + total_variance**2
            )
            for flav in flavors
        }
    else:
        fractions = {}
        fraction_variances = {}
        for flav in _original_flavors(include_unmatched):
            if flav not in h.barable_flavors:
                fractions[flav] = values[flav] / values["all"]
                total_variance = rebinned["all"].variances() / values["all"] ** 2
                fraction_variances[flav] = fractions[flav] ** 2 * (
                    rebinned[flav].variances() / values[flav] ** 2
                    + total_variance**2
                )
                continue
            antiflav = flav + "bar"
            total = values[flav] + values[antiflav]
            total_variance = (
                (rebinned[flav] + rebinned[antiflav]).variances() / total**2
            )
            for component in (flav, antiflav):
                fractions[component] = values[component] / total
                fraction_variances[component] = fractions[component] ** 2 * (
                    rebinned[component].variances() / values[component] ** 2
                    + total_variance**2
                )

    data = {
        "mean": {flav: vals.tolist() for flav, vals in fractions.items()},
        "variance": {
            flav: vals.tolist() for flav, vals in fraction_variances.items()
        },
    }
    objects = _make_fraction_objects(data, pt_bins, jeteta_bins, jeteta_bins.bin_type)
    return label, data, objects


def calculate_flavor_fractions(
    herwig_samples=DEFAULT_HERWIG_SAMPLES,
    pythia_samples=DEFAULT_PYTHIA_SAMPLES,
    eta_binning="Summer20Flavor",
    include_unmatched=False,
    combine_antiflavour=False,
    input_dir="out",
    output_file=None,
    save_plots=True,
    save_root=False,
    ratio_title="Her7/Py8",
):
    herwig_samples = tuple(herwig_samples)
    pythia_samples = tuple(pythia_samples)
    if include_unmatched:
        if "QCD-Py" not in herwig_samples:
            herwig_samples += ("QCD-Py",)
        if "QCD-MG-Py" not in pythia_samples:
            pythia_samples += ("QCD-MG-Py",)
    if len(herwig_samples) != len(pythia_samples):
        raise ValueError("Herwig and Pythia sample lists must have equal lengths.")
    if not herwig_samples:
        raise ValueError("At least one Herwig/Pythia sample pair is required.")

    metadata = fraction_metadata(
        herwig_samples,
        pythia_samples,
        eta_binning,
        include_unmatched,
        combine_antiflavour,
        input_dir,
    )
    if output_file is None:
        output_file = default_cache_path(metadata)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    jeteta_bins = JetEtaBins(eta_binning, absolute=True)
    pt_bins = PtBins("MC_truth")
    flavors = _fraction_flavors(include_unmatched, combine_antiflavour)
    paths = _sample_paths(
        set(herwig_samples).union(pythia_samples), input_dir
    )
    samples_data = {}
    sample_plot = {}
    qfrac_spline_dict = {}
    labels = {}
    for sample in dict.fromkeys((*herwig_samples, *pythia_samples)):
        print(f"Calculating flavor fractions for {sample}")
        label, data, objects = _calculate_sample_fractions(
            sample,
            paths[sample],
            flavors,
            include_unmatched,
            combine_antiflavour,
            jeteta_bins,
            pt_bins,
        )
        samples_data[sample] = data
        labels[sample] = label
        qfracs, qfrac_var, splines, splines_2d, fraction_obj = objects
        sample_plot[sample] = [qfracs, qfrac_var, splines, splines_2d]
        qfrac_spline_dict[sample] = fraction_obj

    payload = {
        "metadata": metadata,
        "binning": {
            "pt_edges": pt_bins.edges.tolist(),
            "pt_centres": pt_bins.centres.tolist(),
            "eta_edges": jeteta_bins.edges.tolist(),
            "eta_centres": jeteta_bins.centres.tolist(),
        },
        "labels": labels,
        "samples": samples_data,
    }
    with output_file.open("w") as cache_file:
        json.dump(payload, cache_file, allow_nan=True)
    print(f"Saved flavor fraction cache: {output_file}")

    if save_plots:
        pltStyle(style="hep", font_frac=1.05)
        from uncertainty_plotters import (
            legend_dict_q_vs_qbar,
            plot_Efractions_ratio,
        )

        (Path("fig") / "fractions").mkdir(parents=True, exist_ok=True)
        ndraw = (
            jeteta_bins.nbins
            if include_unmatched
            or (not combine_antiflavour and eta_binning == "HCalPart")
            else 1
        )
        for her_sample, py_sample in zip(herwig_samples, pythia_samples):
            legenddict = {
                py_sample: labels[py_sample],
                her_sample: labels[her_sample],
            }
            for eta_idx in range(ndraw):
                plot_data = sample_plot
                if not combine_antiflavour:
                    plot_data = {
                        sample: [
                            {
                                flav: item[flav]
                                for flav in _original_flavors(include_unmatched)
                            }
                            if idx in (0, 1)
                            else item
                            for idx, item in enumerate(items)
                        ]
                        for sample, items in sample_plot.items()
                    }
                    plot_Efractions_ratio(
                        plot_data,
                        [py_sample, her_sample],
                        eta_idx,
                        jeteta_bins,
                        pt_bins,
                        legenddict=legenddict,
                        legenddict2=legend_dict_q_vs_qbar,
                        saveplot=True,
                        ratio_lim=(0.8, 1.2),
                        draw_spline=True,
                        ratio_title=ratio_title,
                    )
                else:
                    plot_Efractions_ratio(
                        plot_data,
                        [py_sample, her_sample],
                        eta_idx,
                        jeteta_bins,
                        pt_bins,
                        legenddict=legenddict,
                        saveplot=True,
                        ratio_title=ratio_title,
                    )
                plt.close("all")

    if save_root:
        from save_json import save_fractions_root

        root_data = {
            sample: (
                {
                    flav: np.asarray(values)
                    for flav, values in data["mean"].items()
                },
                {
                    flav: np.asarray(values)
                    for flav, values in data["variance"].items()
                },
            )
            for sample, data in samples_data.items()
        }
        root_path = output_file.with_suffix(".root")
        save_fractions_root(root_data, pt_bins, jeteta_bins, str(root_path))

    return load_fraction_cache(output_file, metadata)


def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--herwig-samples",
        type=parse_sample_list,
        help="Comma-separated sample tags. Defaults to: "
        + ", ".join(DEFAULT_HERWIG_SAMPLES),
    )
    parser.add_argument(
        "--pythia-samples",
        type=parse_sample_list,
        help="Comma-separated sample tags. Defaults to: "
        + ", ".join(DEFAULT_PYTHIA_SAMPLES),
    )
    parser.add_argument(
        "--eta-binning",
        choices=tuple(JERC_Constants.StrToBinsDict(absolute=True)),
        default="Summer20Flavor",
        help="Absolute-eta binning scheme.",
    )
    parser.add_argument(
        "--include-unmatched",
        action="store_true",
        help="Include unmatched jets in the flavor fractions.",
    )
    parser.add_argument(
        "--combine-antiflavour",
        action="store_true",
        help="Combine particle and antiparticle flavors (disabled by default).",
    )
    parser.add_argument("--input-dir", default="out")
    parser.add_argument("--output", dest="output_file")
    parser.add_argument(
        "--ratio-title",
        default="Her7/Py8",
        help="Y-axis title for the flavor-fraction ratio subplot.",
    )
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--save-root", action="store_true")
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    herwig_samples = (
        args.herwig_samples
        if args.herwig_samples is not None
        else DEFAULT_HERWIG_SAMPLES
    )
    pythia_samples = (
        args.pythia_samples
        if args.pythia_samples is not None
        else DEFAULT_PYTHIA_SAMPLES
    )
    calculate_flavor_fractions(
        herwig_samples=herwig_samples,
        pythia_samples=pythia_samples,
        eta_binning=args.eta_binning,
        include_unmatched=args.include_unmatched,
        combine_antiflavour=args.combine_antiflavour,
        input_dir=args.input_dir,
        output_file=args.output_file,
        save_plots=not args.no_plots,
        save_root=args.save_root,
        ratio_title=args.ratio_title,
    )


if __name__ == "__main__":
    main()
