"""Path configuration and lazy initialization; no scientific calculations on import."""

from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
import importlib
import json

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Context:
    settings: dict
    argv: list = field(default_factory=list)
    _cache: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_file(cls, filename, argv=()):
        path = Path(filename).resolve()
        settings = json.loads(path.read_text())
        settings["_config_file"] = str(path)
        return cls(settings, list(argv))

    def root(self, name):
        if name == "repository":
            return PACKAGE_ROOT
        raw = self.settings.get("paths", {}).get(name)
        if not raw:
            raise ValueError(f"Configure paths.{name} in your configuration file.")
        return Path(raw).expanduser().resolve()

    def run_dir(self, name):
        explicit = self.settings.get("runs", {}).get(name)
        if explicit:
            return Path(explicit).expanduser().resolve()
        if name == "network_summary":
            return PACKAGE_ROOT / "data/snapshots/networks"
        if name in {"figure_bank", "regulon_scatter", "target_heatmap"}:
            return PACKAGE_ROOT / "figures" / name
        names = {
            "wgcna": "growth_regulators_2026-10-21_v1",
            "scenic": "ma2016_pyscenic_2026-10-23_v1",
            "maps": "ruslan_mammals_2026-10-22_v1",
        }
        if name == "avian":
            return self.root("avian")
        return self.root("mammals") / "results" / names[name]

    def resolve_path(self, original):
        prefixes = [
            ("/path/to/user-home/chicken/kleiber_birds_gemcat", "avian"),
            ("/path/to/storage/kleiber_birds_gemcat", "mammals"),
            ("/path/to/storage", "resources"),
            ("/path/to/user-home", "tools_home"),
            ("/path/to/workstation", "workstation"),
        ]
        for prefix, key in prefixes:
            if original == prefix or original.startswith(prefix + "/"):
                return self.root(key) / original[len(prefix) :].lstrip("/")
        return Path(original)

    def unit(self, name):
        """Initialize a lazily bound group of existing functions, without running its main."""
        if name not in self._cache:
            from .registry import UNITS

            module_name, builder_name = UNITS[name]
            module = importlib.import_module("." + module_name, __package__)
            self._cache[name] = getattr(module, builder_name)(self)
        return self._cache[name]

    def describe(self, family, stage, argv):
        from .registry import STAGES

        return {
            "workflow": family,
            "stage": stage,
            "backend": STAGES[family][stage],
            "arguments": list(argv),
            "configuration": self.settings.get("_config_file"),
            "paths": self.settings.get("paths", {}),
            "runs": self.settings.get("runs", {}),
            "executes_analysis": False,
        }


def verify_data():
    """Verify packaged bytes; this is integrity checking, not biological replication."""
    import hashlib

    manifest = json.loads((PACKAGE_ROOT / "provenance/data_manifest.json").read_text())
    failures = []
    for item in manifest["files"]:
        path = PACKAGE_ROOT / item["path"]
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            failures.append(item["path"])
    if failures:
        raise ValueError("Missing or changed packaged data: " + ", ".join(failures))
    return {
        "status": "passed",
        "checked_files": len(manifest["files"]),
        "biological_replication": False,
    }
