import argparse
from DEplots.dashboard.dashboard import _cli_wrapper


def _dash_parser(subparsers, name):
    parser = subparsers.add_parser(
        name,
        description="Runs the RAPDOR GUI"
    )
    parser.add_argument(
        '--config',
        type=str,
        help="Path to the config file containing configuration for dash app",
        required=True,
        default=None
    )
    parser.add_argument(
        '--port',
        type=str,
        help="Port to run the Dash server (Default: 8080)",
        default="8080"
    )
    parser.add_argument(
        '--host',
        type=str,
        help="Host IP used by the dash server to serve the application (Default:127.0.0.1)",
        default="127.0.0.1"
    )
    parser.add_argument(
        '--debug',
        action="store_true",
        help="Runs dashboard in debug mode",
    )
    parser.add_argument(
        '--processes',
        type=int,
        help="Number of used cpu cores (Default: 1)",
        default=1
    )
    return parser


class DEPlots:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            "DEPlots suite",
            usage="DEPlots <command> [<args>]"

        )
        self.methods = {
            "Dash": (_dash_parser, _cli_wrapper),
        }
        self.subparsers = self.parser.add_subparsers()
        self.__addparsers()

    def __addparsers(self):
        for name, (parser_add, func) in self.methods.items():
            subp = parser_add(self.subparsers, name)
            subp.set_defaults(func=func)

    def parse_args(self):
        args = self.parser.parse_args()
        return args

    def run(self):
        args = self.parse_args()
        args.func(args)


def main():
    DEPlots().run()


if __name__ == '__main__':
    main()

