"""One entry point for the existing analyses; every stage is explicit."""

import argparse
import json
import os
import subprocess
from .runtime import Context, PACKAGE_ROOT
from .registry import STAGES


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Existing body-mass analyses, consolidated by workflow."
    )
    parser.add_argument("workflow", nargs="?", choices=sorted(STAGES))
    parser.add_argument("stage", nargs="?")
    parser.add_argument("--config", default="config/local.json")
    parser.add_argument(
        "--list", action="store_true", help="List workflows and stages without reading data."
    )
    parser.add_argument(
        "--plan", action="store_true", help="Show the configured operation without running it."
    )
    args, remaining = parser.parse_known_args(argv)
    if args.list or args.workflow is None:
        print(json.dumps({k: list(v) for k, v in STAGES.items()}, indent=2))
        return 0
    if args.stage not in STAGES[args.workflow]:
        parser.error("Choose a stage: " + ", ".join(STAGES[args.workflow]))
    context = Context.from_file(args.config, remaining)
    if args.plan:
        print(json.dumps(context.describe(args.workflow, args.stage, remaining), indent=2))
        return 0
    backend = STAGES[args.workflow][args.stage]
    if isinstance(backend, str):
        context.unit(backend).main()
    else:
        script, operation, run = backend
        env = os.environ.copy()
        if context.settings.get("paths", {}).get("workstation"):
            env["BODYMASS_WORKSTATION"] = str(context.root("workstation"))
        subprocess.run(
            [
                context.settings.get("rscript", "Rscript"),
                "--vanilla",
                str(PACKAGE_ROOT / "R" / script),
                operation,
                str(context.run_dir(run)),
                *remaining,
            ],
            env=env,
            check=True,
        )
    return 0
