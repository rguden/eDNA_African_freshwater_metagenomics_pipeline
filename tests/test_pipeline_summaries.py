"""Run summary steps against small fixtures without bioinformatics dependencies.

Run with: python -m unittest discover -s tests -v
These tests need Bash plus the standard Unix tools used by the pipeline.
"""

import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "eDNA_African_freshwater_metagenomics_pipeline.sh"
)
BASH = shutil.which("bash")


@unittest.skipUnless(BASH, "Bash is required to exercise the pipeline dispatcher")
class PipelineSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.output = self.root / "output"
        self.taxonomy = self.root / "taxonomy"
        self.environment = os.environ.copy()
        self.environment.update(
            OUTPUT_DIR=str(self.output),
            READS_DIR=str(self.root / "reads"),
            METADATA=str(self.root / "metadata.tsv"),
            CAT_TAX=str(self.taxonomy),
            THREADS="1",
        )

    def write(self, relative_path, text):
        path = self.root / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")
        return path

    def run_pipeline(self, *arguments, expected_returncode=0):
        # A timeout also bounds the historical BAT merged-file self-read loop.
        try:
            result = subprocess.run(
                [BASH, str(SCRIPT), *arguments],
                cwd=self.root,
                env=self.environment,
                capture_output=True,
                text=True,
                timeout=8,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            self.fail(f"Pipeline did not finish within 8 seconds: {error.cmd}")
        self.assertEqual(
            result.returncode,
            expected_returncode,
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def read_table(self, relative_path):
        with (self.output / relative_path).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            self.assertTrue(all(None not in row for row in rows), rows)
            return rows

    def write_gtdb_classifications(self, sample, entries):
        rows = ["user_genome\tclassification\tclosest_genome_reference"]
        for genome, species in entries:
            classification = (
                "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
                "o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;"
                f"s__{species}"
            )
            rows.append(f"{genome}\t{classification}\tGCF_000000001.1")
        self.write(
            f"output/gtdbtk/{sample}/classify/gtdbtk.bac120.summary.tsv",
            "\n".join(rows) + "\n",
        )

    def test_no_steps_requires_explicit_selection(self):
        result = self.run_pipeline(expected_returncode=2)
        self.assertIn("Usage:", result.stdout + result.stderr)

    def test_gtdb_summary_preserves_species_with_spaces(self):
        self.write_gtdb_classifications(
            "lake_alpha", [("lake_alpha_bin.1", "Pseudomonas fluorescens")]
        )

        self.run_pipeline("-s", "gtdbtk_summary")

        rows = self.read_table("gtdbtk/gtdbtk_summary.tsv")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["MAG"], "lake_alpha")
        self.assertEqual(rows[0]["Genome"], "lake_alpha_bin.1")
        self.assertEqual(rows[0]["Species"], "s__Pseudomonas fluorescens")
        self.assertEqual(rows[0]["Genus"], "g__Pseudomonas")

    def test_checkm_merge_matches_sample_bin_for_every_quality(self):
        # Local bin identifiers repeat across samples; classified FASTA IDs
        # contain the sample and bin, with no quality suffix.
        cases = [
            ("lake_alpha", "HQ", "Pseudomonas fluorescens"),
            ("lake_beta", "MQ", "Pseudomonas putida"),
            ("lake_gamma", "LQ", "Pseudomonas aeruginosa"),
        ]
        checkm = ["Sample\tBin\tCompleteness\tContamination\tCategory"]
        for sample, quality, species in cases:
            self.write_gtdb_classifications(sample, [(f"{sample}_bin.1", species)])
            checkm.append(f"{sample}\tbin.1\t95.0\t1.0\t{quality}")
        checkm.append("lake_alpha\tbin.2\t40.0\t1.0\tLQ")
        self.write(
            "output/checkm_classified_prok/checkm_summary.tsv",
            "\n".join(checkm) + "\n",
        )

        self.run_pipeline("-s", "gtdbtk_summary,gtdbtk_merge_checkm")

        rows = self.read_table("gtdbtk/gtdbtk_checkm_merged.tsv")
        self.assertEqual(len(rows), 4)
        by_bin = {(row["Sample"], row["Bin"]): row for row in rows}
        for sample, quality, species in cases:
            with self.subTest(sample=sample, quality=quality):
                row = by_bin[(sample, "bin.1")]
                self.assertEqual(row["Category"], quality)
                self.assertEqual(row["Domain"], "d__Bacteria")
                self.assertEqual(row["Species"], f"s__{species}")
        missing = by_bin[("lake_alpha", "bin.2")]
        for rank in ("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"):
            self.assertEqual(missing[rank], "NA")

    def test_checkm_merge_accepts_legacy_ids_and_prefers_current_ids(self):
        self.write_gtdb_classifications(
            "lake_alpha",
            [
                ("lake_alpha_bin.1", "Pseudomonas fluorescens"),
                ("lake_alpha_bin.1_HQ", "Pseudomonas putida"),
                ("lake_alpha_bin.2_MQ", "Pseudomonas aeruginosa"),
            ],
        )
        self.write(
            "output/checkm_classified_prok/checkm_summary.tsv",
            "Sample\tBin\tCompleteness\tContamination\tCategory\n"
            "lake_alpha\tbin.1\t95.0\t1.0\tHQ\n"
            "lake_alpha\tbin.2\t70.0\t1.0\tMQ\n",
        )

        self.run_pipeline("-s", "gtdbtk_summary,gtdbtk_merge_checkm")

        rows = self.read_table("gtdbtk/gtdbtk_checkm_merged.tsv")
        self.assertEqual(len(rows), 2)
        by_bin = {row["Bin"]: row for row in rows}
        self.assertEqual(by_bin["bin.1"]["Species"], "s__Pseudomonas fluorescens")
        self.assertEqual(by_bin["bin.2"]["Species"], "s__Pseudomonas aeruginosa")

    def test_bat_summary_is_stable_when_run_twice(self):
        taxa = [
            (2, "Bacteria", "superkingdom"),
            (10, "Synthetic kingdom", "kingdom"),
            (11, "Pseudomonadota", "phylum"),
            (12, "Gammaproteobacteria", "class"),
            (13, "Pseudomonadales", "order"),
            (14, "Pseudomonadaceae", "family"),
            (15, "Pseudomonas", "genus"),
            (16, "Pseudomonas fluorescens", "species"),
        ]
        self.write(
            "taxonomy/names.dmp",
            "".join(
                f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n"
                for taxid, name, _ in taxa
            ),
        )
        self.write(
            "taxonomy/nodes.dmp",
            "".join(f"{taxid}\t|\t1\t|\t{rank}\t|\n" for taxid, _, rank in taxa),
        )
        bat_header = "# bin\tclassification\treason\tlineage\n"
        self.write(
            "output/bat_prok/lake_alpha/lake_alpha.bin2classification.txt",
            bat_header + "bin.1.fa\ttaxid assigned\tbased on ORFs\t2;10;11;12;13;14;15;16\n",
        )
        self.write(
            "output/bat_prok/lake_beta/lake_beta.bin2classification.txt",
            bat_header + "bin.1.fa\ttaxid assigned\tbased on ORFs\t2;11;12;13;14;15\n",
        )

        self.run_pipeline("-s", "bat_prok_summary_structured")
        merged = self.output / "bat_prok_summary/all_samples.taxonomy_structured.tsv"
        first_output = merged.read_bytes()
        self.run_pipeline("-s", "bat_prok_summary_structured")

        self.assertEqual(merged.read_bytes(), first_output)
        rows = self.read_table("bat_prok_summary/all_samples.taxonomy_structured.tsv")
        self.assertEqual(len(rows), 2)
        by_sample = {row["sample"]: row for row in rows}
        self.assertEqual(set(by_sample), {"lake_alpha", "lake_beta"})
        self.assertEqual(by_sample["lake_alpha"]["species"], "Pseudomonas fluorescens")
        self.assertEqual(by_sample["lake_alpha"]["kingdom"], "Synthetic kingdom")
        self.assertEqual(by_sample["lake_beta"]["species"], "NA")
        self.assertEqual(by_sample["lake_beta"]["kingdom"], "NA")
        self.assertTrue(all(row["bin"] == "bin.1.fa" for row in rows))
        self.assertTrue(all(row["superkingdom"] == "Bacteria" for row in rows))

    def prepare_remag(self, mode):
        self.write("output/tiara/lake_alpha/euk.fa", ">contig_1\nACGT\n")
        self.write("output/mapping/lake_alpha.bam", "synthetic BAM placeholder\n")
        self.write_gtdb_classifications(
            "lake_alpha", [("lake_alpha_bin.1", "Pseudomonas fluorescens")]
        )
        mock = self.write(
            "bin/remag",
            '#!/bin/bash\n'
            'printf "called\\n" >> "$MOCK_REMAG_CALLS"\n'
            'while (( $# )); do\n'
            '    case "$1" in\n'
            '        -o) out="$2"; shift 2 ;;\n'
            '        *) shift ;;\n'
            '    esac\n'
            'done\n'
            'case "$MOCK_REMAG_MODE" in\n'
            '    success)\n'
            '        printf "bin\\n1\\n" > "$out/bins.csv"\n'
            '        mkdir -p "$out/bins"\n'
            '        printf ">contig_1\\nACGT\\n" > "$out/bins/bin.1.fa" ;;\n'
            '    no_bins) exit 0 ;;\n'
            '    failure) printf "synthetic reMAG failure\\n" >&2; exit 23 ;;\n'
            'esac\n',
        )
        mock.chmod(0o755)
        self.environment.update(
            SAMPLE="lake_alpha",
            FORCE="0",
            MOCK_REMAG_MODE=mode,
            MOCK_REMAG_CALLS=str(self.root / "remag_calls.txt"),
            PATH=str(mock.parent) + os.pathsep + self.environment.get("PATH", ""),
        )

    def test_remag_success_continues_to_next_selected_step(self):
        self.prepare_remag("success")

        self.run_pipeline("-s", "remag,gtdbtk_summary")

        self.assertTrue((self.output / "remag/lake_alpha/REMAG_COMPLETE").is_file())
        self.assertEqual(len(self.read_table("gtdbtk/gtdbtk_summary.tsv")), 1)

    def test_remag_valid_no_bins_continues_to_next_selected_step(self):
        self.prepare_remag("no_bins")

        self.run_pipeline("-s", "remag,gtdbtk_summary")

        self.assertTrue((self.output / "remag/lake_alpha/NO_BINS").is_file())
        self.assertEqual(len(self.read_table("gtdbtk/gtdbtk_summary.tsv")), 1)

    def test_remag_completed_sample_skips_tool_and_continues(self):
        self.prepare_remag("failure")
        self.write("output/remag/lake_alpha/bins.csv", "bin\n1\n")

        self.run_pipeline("-s", "remag,gtdbtk_summary")

        self.assertFalse((self.root / "remag_calls.txt").exists())
        self.assertEqual(len(self.read_table("gtdbtk/gtdbtk_summary.tsv")), 1)

    def test_remag_failure_reports_tool_diagnostics_and_exit_status(self):
        self.prepare_remag("failure")

        result = self.run_pipeline(
            "-s", "remag,gtdbtk_summary", expected_returncode=23
        )

        self.assertIn("ERROR: reMAG failed for lake_alpha", result.stdout)
        self.assertIn("Exit code: 23", result.stdout)
        self.assertIn("synthetic reMAG failure", result.stdout)
        self.assertFalse((self.output / "gtdbtk/gtdbtk_summary.tsv").exists())


if __name__ == "__main__":
    unittest.main()
