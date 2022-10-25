import multiprocessing as mp
import os, sys
import argparse
import importlib
from src.utils import LatexReport

"""
instructions:

python3 submit.py launch_local
python3 submit.py --collect

python3 submit.py launch_condor --priority group_u_CMST3.all --queue 1nh
python3 submit.py --collect
"""

# ________________________________________________________________________________
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--config", help="path to config file", default="config/cfg_fcchh_I.py")
    parser.add_argument("--outdir", help="path output directory", default="validation_fcchh_I")
    parser.add_argument("--collect", help="collect jobs and produce report", action="store_true")

    subparsers = parser.add_subparsers(dest="command")

    launch_local = subparsers.add_parser(
        "launch_local", help="launch locally using multi-processing"
    )
    launch_local.add_argument("--dry", help="check local submission", action="store_true")

    launch_condor = subparsers.add_parser("launch_condor", help="launching on condor")
    launch_condor.add_argument("--dry", help="check condor submission", action="store_true")
    launch_condor.add_argument(
        "--priority", help="priority list for condor", default="group_u_CMST3.all"
    )
    launch_condor.add_argument(
        "--queue",
        help="queue for condor",
        choices=[
            "espresso",
            "microcentury",
            "longlunch",
            "workday",
            "tomorrow",
            "testmatch",
            "nextweek",
        ],
        default="longlunch",
    )

    args = parser.parse_args()

    config = args.config
    outdir = os.path.abspath(args.outdir)
    homedir = os.getcwd()

    ## to avoid matplotlib warnings
    os.system("mkdir -p .matplotlib")
    os.system('export MPLCONFIGDIR=".matplotlib"')

    # import config file
    config = os.path.abspath(args.config)
    os.system("cp {} .".format(config))
    config_basename = os.path.basename(config)
    config_module = config_basename.strip(".py")
    print(config_module)
    cfg = importlib.import_module(config_module)

    particles = []
    for plot in cfg.eff_plots:
        for hist in plot.eff_histos:
            part = hist.particle
            if part not in particles:
                particles.append(part)

    for plot in cfg.eff_tag_plots:
        for hist in plot.eff_histos:
            part = hist.particle
            if part not in particles:
                particles.append(part)

    for plot in cfg.reso_plots:
        for hist in plot.res_histos:
            part = hist.particle
            if part not in particles:
                particles.append(part)

    validation_files = dict()
    for p in particles:
        validation_files[p] = "{}/particle_gun_{}/val_{}.root".format(outdir, p.pid, p.pid)
    print("")
    print("=================================================")
    print("")

    print(" CONFIGURATION    = {}".format(config))
    print(" CARD             = {}".format(cfg.card))
    print(" HOME             = {}".format(homedir))
    print(" DELPHES PATH     = {}".format(cfg.delphes_path))
    print(" NJOBS            = {}".format(cfg.njobs))
    print(" NEVTS_PER_JOB    = {}".format(cfg.nevts_per_job))
    print(" PMIN             = {}".format(cfg.pmin))
    print(" PMAX             = {}".format(cfg.pmax))
    print(" ETAMIN           = {}".format(cfg.etamin))
    print(" ETAMAX           = {}".format(cfg.etamax))
    print(" ECM              = {}".format(cfg.ecm))
    part_gun_str = " GUNS PID         ="
    for p in particles:
        part_gun_str = "{} {}".format(part_gun_str, p.pid)
    print(part_gun_str)
    print(" OUTDIR           = {}".format(outdir))
    print("")

    threads = []
    cmdfile = ""
    ## -------- launch local submission
    if args.command == "launch_local":
        print("LOCAL SUBMISSION: ")
        print("")
        print("   NCPUS (available): {}".format(mp.cpu_count()))
        print("")
        ## clean old stuff
        # os.system("rm -rf logs job* {}".format(outdir))

    elif args.command == "launch_condor":

        print("CONDOR SUBMISSION: ")
        print("")
        print("   QUEUE     : {}".format(args.queue))
        print("   PRIORITY  : {}".format(args.priority))
        print("")
        ## clean old stuff
        # os.system("rm -rf logs job* {}".format(outdir))
        os.system("mkdir -p logs")

        cmdfile = """# here goes your shell script
        executable    = src/validation.sh

        # here you specify where to put .log, .out and .err files
        output                = logs/condor.$(ClusterId).$(ProcId).out
        error                 = logs/condor.$(ClusterId).$(ProcId).err
        log                   = logs/condor.$(ClusterId).log

        +AccountingGroup = "{}"
        +JobFlavour    = "{}"
        """.format(
            args.priority, args.queue
        )

    for p in particles:

        jobs_outdir = "{}/particle_gun_{}".format(outdir, p.pid)
        os.system("mkdir -p {}".format(jobs_outdir))

        for i in range(int(cfg.njobs)):
            seed = i

            cmd_args = "{} {} {} {} {} {} {} {} {} {} {} {}".format(
                p.pid,
                cfg.nevts_per_job,
                cfg.pmin,
                cfg.pmax,
                cfg.etamin,
                cfg.etamax,
                seed,
                jobs_outdir,
                homedir,
                cfg.delphes_path,
                cfg.card,
                config,
            )

            if cfg.log:
                cmd_args += " --log"
            # print(args)

            ## -------- launch local submission
            if args.command == "launch_local":

                cmd = "src/validation.sh {}".format(cmd_args)

                if args.dry:
                    print(cmd)
                else:
                    thread = mp.Process(target=run_cmd, args=(cmd,))
                    thread.start()
                    threads.append(thread)

            ## -------- launch condor submission
            elif args.command == "launch_condor":
                cmdfile += 'arguments="{}"\n'.format(cmd_args)
                cmdfile += "queue\n"

    ## run commands in multithreading mode
    if args.command == "launch_local":
        if not args.dry:
            for proc in threads:
                proc.join()

    elif args.command == "launch_condor":
        with open("condor_validation.sub", "w") as f:
            f.write(cmdfile)
        if not args.dry:
            pass

    ## -------- collect jobs
    if args.collect:

        ## first hadd files
        threads = []
        print("")
        print(" collecting jobs ... ")
        print("")

        for p in particles:

            jobs_outdir = "{}/particle_gun_{}".format(outdir, p.pid)
            os.chdir("{}".format(jobs_outdir))
            cmd_hadd = "cd {}; hadd -f val_{}.root".format(jobs_outdir, p.pid)
            for i in range(int(cfg.njobs)):
                cmd_hadd += " val_{}_{}.root".format(p.pid, i)

            # os.system(cmd_hadd)
            os.chdir(homedir)
            thread = mp.Process(target=run_cmd, args=(cmd_hadd,))
            thread.start()
            threads.append(thread)

        ## run hadd in multithreading mode
        for proc in threads:
            proc.join()

        validation_dir = outdir
        report_name = cfg.name
        report_latex = LatexReport(report_name, validation_dir)

        print("")
        print(" producing report ... ")
        print("")

        for section_title, subsections in cfg.report.items():
            report_latex.section(section_title)
            for subsection_title, subsubsections in subsections.items():
                report_latex.subsection(subsection_title)
                for title, plots in subsubsections.items():
                    print(title)
                    report_latex.begin_frame(title)
                    figs = []
                    for plot in plots:
                        fig = plot.plot(validation_files, outdir)
                        figs.append(fig)
                    report_latex.add_figures(figs)

        report_latex.compile()

        print(" ")
        print("{}.pdf report file has been produced".format(report_name))
        print(" ")

    ## do some cleaning
    os.system("rm {}".format(config_basename))
    os.system("rm -rf .matplotlib")


# _____________________________________________________________________________________________________________
def run_cmd(plot_cmd):
    os.system(plot_cmd)


# _______________________________________________________________________________________
if __name__ == "__main__":
    main()
