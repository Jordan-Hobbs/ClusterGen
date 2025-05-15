import argparse
import tomllib
import sys
import logging

import utils

class ProgramController:
    def __init__(self, commands: dict) -> None:
        """Initialize the ProgramController with available commands."""
        self.commands = commands

        # Create parent parser for arguments shared across all subparsers
        self.top_parent_parser = argparse.ArgumentParser(add_help=False)

        self.top_parent_parser.add_argument(
            "-c", "--Config",
            type=str,
            help="Name of the config file used to input required settings."
        )
        self.top_parent_parser.add_argument(
            "-v", "--Verbose",
            action="store_true",
            help="Enable verbose output for debugging."
        )

        # Create parser and subparsers for input arguments
        self.parser = argparse.ArgumentParser(
            prog="ClusterGen",
            parents=[self.top_parent_parser],
            description=(
                "Generate hex clusters of molcules from some input smiles "
                "code. These can then be submitted to some HPC cluster for " 
                "optimisation or some specified calculation. Here by default "
                "the programe is optimised for the SLURM based HPC at the "
                "University of Leeds"
            )
        )

        self.top_subparsers = self.parser.add_subparsers(
            required=True,
            dest="Command",
            help=(
                "Select operating mode: \"cluster_gen\" or "
                "\"cluster_analysis\"."
            )
        )

        self.parse_gen_input()

    def parse_gen_input(self) -> None:
        self.cluster_gen = self.top_subparsers.add_parser(
            "cluster_gen",
            parents=[self.top_parent_parser],
            help="Generate input files for a conformer search."
        )

        # Positional arguments (required, no dashes)
        self.cluster_gen.add_argument(
            "SmilesString",
            type=str,
            help="SMILES string of the compound to be optimised."
        )
        self.cluster_gen.add_argument(
            "NumClusters", 
            type=int,
            help="Number of clusters to generate."
        )
        self.cluster_gen.add_argument(
            "NumRings", 
            type=int,
            help="Number of rings per cluster layer."
        )
        # Optional arguments (with both short and long flags)
        self.cluster_gen.add_argument(
            "-nl", "--NumLayers", 
            type=int,
            help="Number of stacked layers of hex rings."
        )
        self.cluster_gen.add_argument(
            "-o", "--Output", 
            type=str, 
            default="clusters.zip",
            help="Name of the output zip file (default: clusters.zip)."
        )
        self.cluster_gen.add_argument(
            "-nc", "--NumCPUs",
            type=int,
            default=8,
            help="Number of CPUs allocated to the HPC calculation."
        )
        self.cluster_gen.add_argument(
            "-t", "--RunTime",
            type=str,
            default="24:00:00",
            help="Maximum runtime allowed on the HPC system."
        )
        self.cluster_gen.add_argument(
            "-e", "--Email",
            type=str,
            help=(
                "Email address to receive updates on job start, finish, "
                "and failure."
            )
        )
        self.cluster_gen.add_argument(
            "-a", "--NumAP", 
            type=int, 
            default=0,
            help="Number of inverted molecules per cluster (default: 0)."
        )
        self.cluster_gen.add_argument(
            "-d", "--MolSep", 
            type=float,
            default=6.0,
            help="Radial separation between molecules (default: 6.0 Å)."
        )
        self.cluster_gen.add_argument(
            "-m", "--MinSep", 
            type=float, 
            default=2.5,
            help="Minimum allowed interatomic distance (default: 2.5 Å)."
        )
        self.cluster_gen.add_argument(
            "-cl", "--CRESTOptLevel",
            type=str,
            default="extreme",
            choices=[
                "crude", "vloose", "loose", "normal",
                "tight", "vtight", "extreme"
            ],
            help="Convergence conditions for CREST optimisation steps."
        )
        self.cluster_gen.add_argument(
            "-cm", "--CRESTMethod",
            type=str,
            default="gfnff",
            choices=["gfnff", "gfn0", "gfn1", "gfn2"],
            help=(
                "Tight binding method used by CREST for optimisation and "
                "conformer searching."
            )
        )
        self.cluster_gen.add_argument(
            "-rn", "--RDNumConf",
            type=int,
            default=100,
            help="Number of conformers for the initial RDKit search (default: 100)."
        )
        self.cluster_gen.set_defaults(func=self.commands["cluster_gen"])

    def load_from_config(self) -> None:
        """Load configuration parameters from a TOML file if provided."""
        try:
            with open(self.args.Config, "rb") as file:
                config = tomllib.load(file)
            logging.info("Config file found.")
        except FileNotFoundError:
            logging.warning(
                f"Config file not found. Is \"{self.args.Config}\" a valid "
                "path? Using default parameters instead."
            )
            return

        flat_config = utils.flatten_dict(config.get(self.args.Command))

        # Get command-line arguments passed (excluding --Config)
        cmdline_args = [
            arg.lstrip("-") for arg in sys.argv[2:]
            if arg.startswith(("-", "--")) and arg not in ("-c", "--Config")
        ]
        
        subparser = self.top_subparsers.choices[self.args.Command]

        # Convert short arg names to long form
        short_to_long = self._build_short_to_long_map(subparser)
        cmdline_args = [
            short_to_long.get(arg.lstrip("-"), arg.lstrip("-")) 
            for arg in cmdline_args
        ]

        # Filter out config keys that were already supplied via CLI
        for key in cmdline_args:
            flat_config.pop(key, None)

        # Inject remaining config values into the args Namespace
        for key, value in flat_config.items():
            setattr(self.args, key, value)

        logging.info("Config parameters loaded successfully.")

    def _build_short_to_long_map(self, subparser) -> dict:
        """Build a mapping of short flags to long flags."""
        mapping = {}
        for action in subparser._actions:
            short = long = None
            for option in action.option_strings:
                if option.startswith("--"):
                    long = option.lstrip("-")
                elif option.startswith("-"):
                    short = option.lstrip("-")
            if short and long:
                mapping[short] = long
        return mapping

    def get_args(self) -> argparse.Namespace:
        """Parse and return command-line arguments."""
        self.args = self.parser.parse_args()

        # Set logging level
        if getattr(self.args, "Verbose", False):
            logging.basicConfig(level=logging.DEBUG, format='%(message)s')
        else:
            logging.basicConfig(level=logging.INFO, format='%(message)s')

        return self.args
