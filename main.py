import cluster_gen
import cmd_parser

def main():
    # Sets command functions
    commands={
        "cluster_gen": cluster_gen.cluster_gen
    }
    
    # Initialize parser
    parser = cmd_parser.ProgramController(commands)
    args = parser.get_args()

    # Loads config file if config file provided
    if args.Config:
        print(
            "Config input detected. Attempting to load parameters from "
            f"\"{args.Config}\""
        )
        parser.load_from_config()
    else:
        print(
            "No config file provided. Using command line and defaults "
            "instead."
        )

    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()