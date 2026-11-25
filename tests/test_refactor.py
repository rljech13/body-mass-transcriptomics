"""Regression checks for moved existing calculations; no model fitting or new analyses."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from bodymass.runtime import Context, PACKAGE_ROOT


class RefactorTests(unittest.TestCase):
    def context(self):
        return Context(
            {
                "paths": {
                    k: str(Path(tempfile.gettempdir()) / "bodymass-test-inputs" / k)
                    for k in ["avian", "mammals", "resources", "tools_home", "workstation"]
                }
            }
        )

    def test_help_and_plan_do_not_read_scientific_inputs(self):
        proc = subprocess.run(
            [sys.executable, "-m", "bodymass", "--list"],
            cwd=PACKAGE_ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        stages = json.loads(proc.stdout)
        self.assertEqual(set(stages), {"avian", "groups", "networks", "scenic", "figures"})
        with tempfile.TemporaryDirectory() as folder:
            config = Path(folder) / "config.json"
            config.write_text(json.dumps(self.context().settings))
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "bodymass",
                    "scenic",
                    "grn",
                    "--config",
                    str(config),
                    "--plan",
                    "--seed",
                    "20261023",
                ],
                cwd=PACKAGE_ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertFalse(json.loads(proc.stdout)["executes_analysis"])
            self.assertEqual(sorted(p.name for p in Path(folder).iterdir()), ["config.json"])

    def test_existing_bird_fit_functions_preserve_all_eighteen_frozen_results(self):
        poster = json.loads((PACKAGE_ROOT / "data/frozen/poster_inputs.json").read_text())
        for unit in ["avian.53_three_checks", "avian.78_holding_cores"]:
            api = self.context().unit(unit)
            for expected in poster["tests"]:
                if expected["version"] != "corrected":
                    continue
                rows = {
                    r["species"]: r
                    for r in poster["birds"]
                    if r["tissue"] == expected["tissue"] and r["set"] == expected["set"]
                }
                values = [rows[s]["log10_sum"] for s in api.SPECIES_ORDER]
                r, p = api.fit(values)
                self.assertAlmostEqual(r, expected["r"], places=10)
                self.assertAlmostEqual(p, expected["p"], places=10)

    def test_grn_worker_preserves_original_arguments(self):
        from bodymass import scenic

        received = {}

        def original_backend(**kwargs):
            received.update(kwargs)
            return "called"

        class Matrix:
            def __getitem__(self, index):
                return ("column", index[1])

        scenic._GRN_STATE.update(
            infer_partial_network=original_backend,
            SGBM_KWARGS={"loss": "squared_error"},
            EARLY_STOP_WINDOW_LENGTH=25,
            TFX="TF matrix",
            TFN=["TF1"],
            GN=["G0", "G1"],
            X=Matrix(),
            SEED=20261023,
        )
        self.assertEqual(scenic._infer_target(1), "called")
        self.assertEqual(received["regressor_type"], "GBM")
        self.assertEqual(received["target_gene_name"], "G1")
        self.assertEqual(received["target_gene_expression"], ("column", 1))
        self.assertEqual(received["seed"], 20261023)
        self.assertFalse(received["include_meta"])
        self.assertEqual(received["early_stop_window_length"], 25)
        import multiprocessing as mp

        if "fork" in mp.get_all_start_methods():
            with mp.get_context("fork").Pool(1) as pool:
                self.assertEqual(pool.map(scenic._infer_target, [1]), ["called"])


if __name__ == "__main__":
    unittest.main()
