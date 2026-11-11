#!/usr/bin/env python
"""Derive flavor-dependent JEC uncertainties from response and fraction inputs."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common_binning import JERC_Constants
from RatioPlotFitRes import CorrectionEvaluator, RatioPlotFitRes, get_uncertainty
from correction_fitter import correction_fitter
from flavour_fractions import (
    DEFAULT_HERWIG_SAMPLES,
    DEFAULT_PYTHIA_SAMPLES,
    calculate_flavor_fractions,
    default_cache_path,
    fraction_metadata,
    load_fraction_cache,
    parse_sample_list,
)
from uncertainty_helpers import ptmax_global, ptmin_global


FLAVOR_UNCERTAINTY_NAMES = {
    "g": "PureGluon",
    "c": "PureCharm",
    "b": "PureBottom",
    "ud": "PureUpDown",
    "q": "PureQuark",
    "s": "PureStrange",
    "u": "PureUp",
    "d": "PureDown",
    "cbar": "PureAntiCharm",
    "bbar": "PureAntiBottom",
    "udbar": "PureAntiUpDown",
    "qbar": "PureAntiQuark",
    "sbar": "PureAntiStrange",
}


def _uncertainty_flavors(combine_antiflavour):
    flavors = ["g", "c", "b", "ud", "s"]
    if combine_antiflavour:
        return flavors
    return [
        item
        for flav in flavors
        if flav != "g"
        for item in (flav, flav + "bar")
    ]


def _correction_configuration(combine_antiflavour):
    correction_dir = Path("Summer20UL18_V2_MC")
    correction_prefix = "Summer20UL18_V2_MC_L5Flavor_AK4PFchs"
    if combine_antiflavour:
        correction_base = correction_prefix + "_MGQCD"
    else:
        correction_base = (
            correction_prefix + "_split_antiflav_Summer20Flavor_old"
        )
    return correction_dir, correction_base


def _load_correction_evaluator(
    eta_binning,
    flavors,
    correction_base,
    correction_dir,
    correction_type,
    combine_antiflavour,
    rerun_fits,
    save_fit_plots,
):
    if rerun_fits:
        if not combine_antiflavour:
            raise ValueError(
                "Refitting split flavor/antiflavor corrections is not supported "
                "by correction_fitter.py."
            )
        correction_fitter(
            saveplots=save_fit_plots,
            do_simfit=True,
            do_Mikkofit=False,
            correction_for="MGQCD_Py",
            eta_binning=eta_binning,
            combine_antiflavour=True,
        )
        correction_fitter(
            saveplots=save_fit_plots,
            do_simfit=True,
            do_Mikkofit=False,
            correction_for="MGQCD_Her",
            eta_binning=eta_binning,
            combine_antiflavour=True,
        )

    eta_suffix = "_" + eta_binning if eta_binning != "Summer20Flavor" else ""
    correction_files = (
        correction_dir / f"{correction_base}{eta_suffix}.txt",
        correction_dir / f"{correction_base}_Her{eta_suffix}.txt",
    )
    missing_corrections = [
        str(path) for path in correction_files if not path.is_file()
    ]
    if missing_corrections:
        message = "Missing correction payloads:\n" + "\n".join(
            missing_corrections
        )
        if combine_antiflavour:
            message += "\nUse --rerun-fits to produce the standard Run 2 payloads."
        raise FileNotFoundError(
            message
        )

    evaluator = CorrectionEvaluator(
        binning=eta_binning,
        flavors=flavors,
        correction_txt_dir=str(correction_dir) + "/",
        correction_txt=correction_base,
        inverse=False,
        ptmin_global=ptmin_global,
        ptmax_global=ptmax_global,
        use_corrections="S" if correction_type == "simfit" else correction_type,
    )
    return evaluator


def _load_ratio_evaluator(
    evaluator,
    eta_binning,
    flavors,
    jeteta_bins,
    pt_bins,
):
    from uncertainty_plotters import plot_ratio_comparisons_samples

    eta_binning_str = "_" + eta_binning if eta_binning != "HCalPart" else ""
    ratio_fit = RatioPlotFitRes(
        binning=eta_binning,
        ptmin_global=ptmin_global,
        ptmax_global=ptmax_global,
    )
    correction_name = evaluator.correction_txt + evaluator.eta_binning_str
    evaluator_names = {
        "Sum20Py": correction_name,
        "Sum20Her": correction_name + "_Her",
    }
    for flavor in flavors:
        fit_results = [
            plot_ratio_comparisons_samples(
                flavor,
                eta_idx,
                jeteta_bins,
                pt_bins.centres,
                eta_binning_str,
                evaluator=evaluator,
                evaluator_names=evaluator_names,
                plotnewfit=True,
                plotcorrectionratios=True,
                show_original_uncertainties=True,
                max_point_fit_idx=5,
            )
            for eta_idx in range(jeteta_bins.nbins)
        ]
        ratio_fit.add_flavor(flavor, np.asarray(fit_results))
    return ratio_fit


def calculate_flavor_uncertainties(
    herwig_samples=DEFAULT_HERWIG_SAMPLES,
    pythia_samples=DEFAULT_PYTHIA_SAMPLES,
    eta_binning="Summer20Flavor",
    combine_antiflavour=False,
    input_dir="out",
    fraction_cache=None,
    recompute_fractions=False,
    correction_type="T",
    correction_dir=None,
    correction_base=None,
    rerun_fits=False,
    save_plots=True,
    output_file=None,
    uncertainty_files=(),
    output_dir=None,
    overwrite_output=False,
    divide_herwig_pythia=False,
):
    if len(herwig_samples) != len(pythia_samples):
        raise ValueError("Herwig and Pythia sample lists must have equal lengths.")
    if not herwig_samples:
        raise ValueError("At least one Herwig/Pythia sample pair is required.")
    if correction_type not in {"J", "T", "ratio", "simfit"}:
        raise ValueError("correction_type must be one of: J, T, ratio, simfit.")

    metadata = fraction_metadata(
        herwig_samples,
        pythia_samples,
        eta_binning,
        include_unmatched=False,
        combine_antiflavour=combine_antiflavour,
        input_dir=input_dir,
    )
    if fraction_cache is None:
        fraction_cache = default_cache_path(metadata)
    fraction_cache = Path(fraction_cache)

    fractions = None
    if not recompute_fractions:
        try:
            fractions = load_fraction_cache(fraction_cache, metadata)
            print(f"Using matching flavor fraction cache: {fraction_cache}")
        except FileNotFoundError:
            print(f"Flavor fraction cache not found: {fraction_cache}; recomputing.")
        except json.JSONDecodeError as error:
            print(
                f"Flavor fraction cache is invalid ({error}); recomputing: "
                f"{fraction_cache}"
            )
        except ValueError as error:
            print(f"{error} Recomputing fractions.")

    if fractions is None:
        fractions = calculate_flavor_fractions(
            herwig_samples=herwig_samples,
            pythia_samples=pythia_samples,
            eta_binning=eta_binning,
            include_unmatched=False,
            combine_antiflavour=combine_antiflavour,
            input_dir=input_dir,
            output_file=fraction_cache,
            save_plots=save_plots,
        )

    qfrac_spline_dict = fractions["qfrac_spline_dict"]
    required_samples = {
        "QCD-MG-Her",
        "QCD-MG-Py",
        "DY-MG-Her",
        "DY-MG-Py",
    }
    missing_samples = sorted(required_samples - qfrac_spline_dict.keys())
    if missing_samples:
        raise ValueError(
            "The uncertainty calculation requires these cached fraction "
            f"samples: {', '.join(missing_samples)}."
        )

    flavors = _uncertainty_flavors(combine_antiflavour)
    jeteta_bins = fractions["jeteta_bins"]
    pt_bins = fractions["pt_bins"]

    default_dir, default_base = _correction_configuration(combine_antiflavour)
    if rerun_fits and (correction_dir is not None or correction_base is not None):
        raise ValueError(
            "--rerun-fits writes the standard Run 2 correction payload paths; "
            "do not combine it with custom correction directory or base values."
        )
    correction_dir = Path(correction_dir) if correction_dir else default_dir
    correction_base = correction_base or default_base
    evaluator = _load_correction_evaluator(
        eta_binning,
        flavors,
        correction_base,
        correction_dir,
        correction_type,
        combine_antiflavour,
        rerun_fits,
        save_plots,
    )
    if correction_type == "ratio":
        if not combine_antiflavour:
            raise ValueError(
                "Ratio correction fits are not supported for split "
                "flavor/antiflavor uncertainties."
            )
        difference_evaluator = _load_ratio_evaluator(
            evaluator,
            eta_binning,
            flavors,
            jeteta_bins,
            pt_bins,
        )
    else:
        difference_evaluator = evaluator

    qcd_pair_idx = next(
        (
            idx
            for idx, (her_sample, py_sample) in enumerate(
                zip(herwig_samples, pythia_samples)
            )
            if "QCD" in her_sample and "QCD" in py_sample
        ),
        None,
    )
    if qcd_pair_idx is None:
        raise ValueError(
            "Herwig/Pythia samples must include a matching QCD pair to form "
            "the reference difference."
        )

    qcd_her = herwig_samples[qcd_pair_idx]
    qcd_py = pythia_samples[qcd_pair_idx]
    uncertainty_values = {flavor: np.empty((jeteta_bins.nbins, pt_bins.nbins)) for flavor in flavors}
    for eta_idx, eta_center in enumerate(jeteta_bins.centres):
        eta_values = np.full(pt_bins.nbins, eta_center)
        eta_values_reference = np.full(pt_bins.nbins, jeteta_bins.edges[0])
        qcd_difference = difference_evaluator.resum_ratio_to_mix(
            eta_values,
            pt_bins.centres,
            qfrac_spline_dict[qcd_her],
            qfrac_spline_dict[qcd_py],
            divide_herwig_pythia,
        )
        curves = get_uncertainty(
            eta_values,
            pt_bins.centres,
            eta_values_reference,
            qcd_difference,
            evaluator=difference_evaluator,
            qfrac_spline_dict=qfrac_spline_dict,
            divideHerPy=divide_herwig_pythia,
        )
        for flavor in flavors:
            if combine_antiflavour:
                uncertainty_values[flavor][eta_idx] = (
                    curves[flavor + "100"] - curves["Rref"]
                )
            else:
                uncertainty_values[flavor][eta_idx] = curves[flavor + "100"]

    uncertainty_payload = {
        "metadata": {
            "eta_binning": eta_binning,
            "combine_antiflavour": combine_antiflavour,
            "correction_type": correction_type,
            "divide_herwig_pythia": divide_herwig_pythia,
            "herwig_samples": list(herwig_samples),
            "pythia_samples": list(pythia_samples),
            "fraction_cache": str(fraction_cache),
            "correction_dir": str(correction_dir),
            "correction_base": correction_base,
        },
        "binning": {
            "pt_edges": pt_bins.edges.tolist(),
            "pt_centres": pt_bins.centres.tolist(),
            "eta_edges": jeteta_bins.edges.tolist(),
            "eta_centres": jeteta_bins.centres.tolist(),
        },
        "uncertainties": {
            flavor: values.tolist()
            for flavor, values in uncertainty_values.items()
        },
    }
    if output_file is None:
        sample_tags = "-".join(
            (*herwig_samples, *pythia_samples)
        )
        split_suffix = "" if combine_antiflavour else "_split_antiflav"
        output_file = Path("out_txt") / (
            f"flavour_uncertainties_{correction_type}_{sample_tags}_"
            f"{eta_binning}{split_suffix}.json"
        )
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as result_file:
        json.dump(uncertainty_payload, result_file, allow_nan=True)
    print(f"Saved flavor uncertainties: {output_file}")

    if save_plots:
        plot_dir = Path("fig") / "uncertainty" / correction_type
        plot_dir.mkdir(parents=True, exist_ok=True)
        for eta_idx, eta_center in enumerate(jeteta_bins.centres):
            fig, ax = plt.subplots()
            for flavor in flavors:
                ax.plot(
                    pt_bins.centres,
                    uncertainty_values[flavor][eta_idx] * 100,
                    label=flavor,
                )
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_xscale("log")
            ax.set_xlabel(r"$p_{T}$ (GeV)")
            ax.set_ylabel("JEC uncertainty (%)")
            ax.set_title(jeteta_bins.idx2plot_str(eta_idx))
            ax.legend()
            fig.tight_layout()
            plot_name = plot_dir / f"uncertainty_eta{eta_idx}.png"
            fig.savefig(plot_name)
            fig.savefig(plot_name.with_suffix(".pdf"))
            plt.close(fig)
        print(f"Saved uncertainty plots under {plot_dir}")

    if uncertainty_files:
        new_uncertainties = _make_payload_sections(
            uncertainty_values,
            pt_bins,
            jeteta_bins,
        )
        _write_uncertainty_files(
            uncertainty_files,
            new_uncertainties,
            combine_antiflavour,
            output_dir,
            overwrite_output,
        )

    return uncertainty_payload


def _make_payload_sections(uncertainty_values, pt_bins, jeteta_bins):
    sections = {}
    for flavor, values in uncertainty_values.items():
        if flavor not in FLAVOR_UNCERTAINTY_NAMES:
            raise ValueError(
                f"No uncertainty-source section mapping is defined for {flavor!r}."
            )
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"Cannot write non-finite uncertainty values for {flavor}."
            )
        name = "Flavor" + FLAVOR_UNCERTAINTY_NAMES[flavor]
        negative_eta_rows = [
            [-jeteta_bins.edges[idx + 1], -jeteta_bins.edges[idx], pt_bins.nbins * 3]
            + list(
                np.column_stack(
                    (
                        pt_bins.centres,
                        values[idx],
                        values[idx],
                    )
                ).reshape(-1)
            )
            for idx in range(jeteta_bins.nbins - 1, -1, -1)
        ]
        positive_eta_rows = [
            [jeteta_bins.edges[idx], jeteta_bins.edges[idx + 1], pt_bins.nbins * 3]
            + list(
                np.column_stack(
                    (
                        pt_bins.centres,
                        values[idx],
                        values[idx],
                    )
                ).reshape(-1)
            )
            for idx in range(jeteta_bins.nbins)
        ]
        sections[name] = np.asarray(negative_eta_rows + positive_eta_rows)
    return sections


def _render_sections(sections, digits):
    lines = []
    for section, rows in sections.items():
        lines.extend([f"[{section}]\n", '{1 JetEta 1 JetPt "" Correction JECSource}\n'])
        lines.extend(
            " ".join(str(value) for value in np.round(row, digits)) + "\n"
            for row in rows
        )
    return lines


def _replace_sections(source_path, sections, output_path, digits, overwrite):
    source_path = Path(source_path)
    output_path = Path(output_path)
    if not source_path.is_file():
        raise FileNotFoundError(f"Uncertainty source file not found: {source_path}")
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Pass --overwrite-output to replace it."
        )

    source_lines = source_path.read_text().splitlines(keepends=True)
    output_lines = []
    seen = set()
    idx = 0
    while idx < len(source_lines):
        line = source_lines[idx]
        if line.startswith("["):
            section = line.strip()[1:-1]
            if section in sections:
                output_lines.extend(_render_sections({section: sections[section]}, digits))
                seen.add(section)
                idx += 1
                while idx < len(source_lines) and not source_lines[idx].startswith("["):
                    idx += 1
                continue
        output_lines.append(line)
        idx += 1

    missing = set(sections) - seen
    if missing:
        raise ValueError(
            f"{source_path} does not contain expected sections: "
            + ", ".join(sorted(missing))
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(output_lines))
    print(f"Saved uncertainty payload: {output_path}")


def _write_uncertainty_files(
    uncertainty_files,
    sections,
    combine_antiflavour,
    output_dir,
    overwrite,
):
    for filename in uncertainty_files:
        source_path = Path(filename)
        if output_dir:
            output_path = Path(output_dir) / (
                source_path.stem
                + ("_run2flavor" if combine_antiflavour else "_run2flavor_antiflavor")
                + source_path.suffix
            )
        else:
            output_path = source_path.with_name(
                source_path.stem
                + ("_run2flavor" if combine_antiflavour else "_run2flavor_antiflavor")
                + source_path.suffix
            )
        if combine_antiflavour:
            _replace_sections(
                source_path,
                sections,
                output_path,
                digits=4,
                overwrite=overwrite,
            )
        else:
            if output_path.exists() and not overwrite:
                raise FileExistsError(
                    f"Output already exists: {output_path}. "
                    "Pass --overwrite-output to replace it."
                )
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text("".join(_render_sections(sections, digits=6)))
            print(f"Saved uncertainty payload: {output_path}")


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
        "--combine-antiflavour",
        action="store_true",
        help="Combine particle and antiparticle flavors (disabled by default).",
    )
    parser.add_argument("--input-dir", default="out")
    parser.add_argument("--fraction-cache")
    parser.add_argument("--recompute-fractions", action="store_true")
    parser.add_argument("--correction-type", choices=("J", "T", "ratio", "simfit"), default="T")
    parser.add_argument("--correction-dir")
    parser.add_argument(
        "--correction-base",
        help="Correction filename base, excluding the Her suffix and eta-binning suffix.",
    )
    parser.add_argument("--rerun-fits", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--output")
    parser.add_argument(
        "--uncertainty-files",
        nargs="*",
        default=(),
        help="Existing uncertainty payloads to copy and update. No payloads are written unless specified.",
    )
    parser.add_argument("--output-dir")
    parser.add_argument("--overwrite-output", action="store_true")
    parser.add_argument("--divide-herwig-pythia", action="store_true")
    return parser


def main(argv=None):
    args = _build_parser().parse_args(argv)
    calculate_flavor_uncertainties(
        herwig_samples=(
            args.herwig_samples
            if args.herwig_samples is not None
            else DEFAULT_HERWIG_SAMPLES
        ),
        pythia_samples=(
            args.pythia_samples
            if args.pythia_samples is not None
            else DEFAULT_PYTHIA_SAMPLES
        ),
        eta_binning=args.eta_binning,
        combine_antiflavour=args.combine_antiflavour,
        input_dir=args.input_dir,
        fraction_cache=args.fraction_cache,
        recompute_fractions=args.recompute_fractions,
        correction_type=args.correction_type,
        correction_dir=args.correction_dir,
        correction_base=args.correction_base,
        rerun_fits=args.rerun_fits,
        save_plots=not args.no_plots,
        output_file=args.output,
        uncertainty_files=args.uncertainty_files,
        output_dir=args.output_dir,
        overwrite_output=args.overwrite_output,
        divide_herwig_pythia=args.divide_herwig_pythia,
    )


if __name__ == "__main__":
    main()
