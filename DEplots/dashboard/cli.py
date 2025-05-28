import DEplots.dashboard
from DEplots.dashboard import get_data


def cli_wrapper(
        config_file: str = None,
        run_dir: str = None,
        debug: bool = False,
        port: int = 8080,
        host: str = "127.0.0.1",
        processes: int = 1
):
    DEplots.dashboard.DASH_DATA, DEplots.dashboard.CONFIG = get_data(config_file, run_dir)
    from DEplots.dashboard.app import app, get_layout

    app.layout = get_layout()
    app.run(debug=debug, port=port, host=host, processes=processes, threaded=False)


def _cli_wrapper(args):
    print(args)
    cli_wrapper(args.config, args.run_dir, args.debug, args.port, args.host, args.processes)





if __name__ == '__main__':
    config_file = "/home/rabsch/PythonProjects/DEPlots/testData/config.yaml"
    rd = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/"
    #rd = "/home/rabsch/PythonProjects/SynDESeq2/RUNS/"
    cli_wrapper(config_file=config_file, run_dir=rd, debug=True, processes=3)
