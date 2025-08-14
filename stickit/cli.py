import argparse, json
from .pipeline import stic_generation
from .utils import load_config

def main():
    ap = argparse.ArgumentParser(prog="stickit")
    sp = ap.add_subparsers(dest="cmd", required=True)

    runp = sp.add_parser("run", help="Run STIC generation from YAML config")
    runp.add_argument("--config", required=True, help="YAML config path")

    args = ap.parse_args()
    if args.cmd == "run":
        cfg = load_config(args.config)
        sticsets = stic_generation(config=cfg)  # will read input from cfg['io']['input']
        summary = {s.parent_key: sum(len(st.conformers) for st in s.stics) for s in sticsets}
        print(json.dumps(summary, indent=2))

