import logging
import cluster_gen
import cmd_parser

def main():
    # Setup logging to file and console
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filename='cluster_gen.log',
        filemode='w'
    )

    # Also log to console at INFO level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logging.getLogger().addHandler(console_handler)

    # Sets command functions
    commands={
        "cluster_gen": cluster_gen.cluster_gen
    }
    
    # Initialize parser
    parser = cmd_parser.ProgramController(commands)
    args = parser.get_args()

    # Loads config file if config file provided
    if args.Config:
        logging.info(
            "Config input detected. Attempting to load parameters from "
            f"\"{args.Config}\""
        )
        parser.load_from_config()
    else:
        logging.info(
            "No config file provided. Using command line and defaults "
            "instead."
        )

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
