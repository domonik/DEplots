import DEplots.dashboard
from DEplots.dashboard import get_data, get_coverage_data


def cli_wrapper(
        config_file: str = None,
        run_dir: str = None,
        debug: bool = False,
        port: int = 8080,
        host: str = "127.0.0.1",
        processes: int = 1
):
    DEplots.dashboard.DASH_DATA = get_data(config_file, run_dir)
    cov_design, cov_data, gff, line_mapping = get_coverage_data(config_file)
    if cov_data is not None:
        DEplots.dashboard.COVERAGE_DATA = cov_data
        DEplots.dashboard.COVERAGE_DESIGN = cov_design
        DEplots.dashboard.GFF = gff
        DEplots.dashboard.LINE_MAPPING = line_mapping
    from DEplots.dashboard.app import app, get_layout

    app.layout = get_layout()
    app.run(debug=debug, port=port, host=host, processes=processes, threaded=False)


def _cli_wrapper(args):
    cli_wrapper(args.config, args.run_dir, args.debug, args.port, args.host, args.processes)





if __name__ == '__main__':
    config_file = "/home/rabsch/PythonProjects/DEPlots/testData/config.yaml"
    rd = "/home/rabsch/PythonProjects/RlocSeq/Pipeline/RUNS/"
    cli_wrapper(config_file=config_file, run_dir=rd, debug=True, processes=3)
