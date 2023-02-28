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
    parser.add_argument("--outdir", help="path output directory", default="validation_samples")

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

    collect = subparsers.add_parser("collect", help="collect jobs and plot")
    collect.add_argument(
        "--no_hadd", help="collect without hadding", action="store_true", default=False
    )

    args = parser.parse_args()

    config = args.config
    outdir = os.path.abspath(args.outdir)
    homedir = os.getcwd()
    delphes_path = os.path.abspath("{}/..".format(homedir))

    ## to avoid matplotlib warnings
    os.system("mkdir -p .matplotlib")
    os.system('export MPLCONFIGDIR=".matplotlib"')

    # import config file
    config = os.path.abspath(args.config)
    os.system("cp {} .".format(config))
    config_basename = os.path.basename(config)
    config_module = config_basename.strip(".py")
    proc_dir = config_module
    sample_dir = "{}/{}".format(outdir, proc_dir)
    cfg = importlib.import_module(config_module)

    particles = []

    plots = cfg.eff_plots + cfg.eff_tag_plots + cfg.reso_plots

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
        validation_files[p] = "{}/particle_gun_{}/val_{}.root".format(sample_dir, p.pid, p.pid)
    print("")
    print("=================================================")
    print("")

    print(" CONFIGURATION    = {}".format(config))
    print(" CARD             = {}".format(cfg.card))
    print(" HOME             = {}".format(homedir))
    print(" DELPHES PATH     = {}".format(delphes_path))
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
    print(" OUTDIR           = {}".format(sample_dir))
    print("")

    pool = mp.Pool()
    cmdfile = ""
    ## -------- launch local submission
    if args.command == "launch_local":
        print("LOCAL SUBMISSION: ")
        print("")
        print("   NCPUS (available): {}".format(mp.cpu_count()))
        print("")

        ## clean previous runs and logs
        os.system("rm -rf logs/* job* {} {}".format(sample_dir, proc_dir))

    elif args.command == "launch_condor":

        print("CONDOR SUBMISSION: ")
        print("")
        print("   QUEUE     : {}".format(args.queue))
        print("   PRIORITY  : {}".format(args.priority))
        print("")

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

        ## clean previous runs and logs
        os.system("rm -rf logs/* job* {} {}".format(sample_dir, proc_dir))

    for p in particles:

        jobs_outdir = "{}/particle_gun_{}".format(sample_dir, p.pid)
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
                delphes_path,
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
                    pool.apply_async(run_cmd, args=(cmd,))

            ## -------- launch condor submission
            elif args.command == "launch_condor":
                cmdfile += 'arguments="{}"\n'.format(cmd_args)
                cmdfile += "queue\n"

    ## run commands in multithreading mode
    if args.command == "launch_local":
        if not args.dry:
            pool.close()
            pool.join()

    elif args.command == "launch_condor":
        with open("condor_validation.sub", "w") as f:
            f.write(cmdfile)
        if not args.dry:
            os.system("condor_submit condor_validation.sub")

    ## -------- collect jobs
    # if args.collect:
    if args.command == "collect":

        if not args.no_hadd:
            ## first hadd files
            pool = mp.Pool()
            print("")
            print(" collecting jobs ... ")
            print("")

            for p in particles:

                jobs_outdir = "{}/particle_gun_{}".format(sample_dir, p.pid)
                os.chdir("{}".format(jobs_outdir))
                cmd_hadd = "cd {}; hadd -f val_{}.root val_{}_*.root".format(
                    jobs_outdir, p.pid, p.pid
                )
                os.chdir(homedir)
                pool.apply_async(run_cmd, args=(cmd_hadd,))

            ## run hadd in multithreading mode
            pool.close()
            pool.join()

        print("")
        print(" producing plots ... ")
        print("")

        ## produce plots in multithreading mode
        os.system("rm -rf  {}".format(proc_dir))
        os.system("mkdir -p {}".format(proc_dir))

        name_path_dict = dict()

        pool = mp.Pool()
        manager = mp.Manager()
        name_path_dict = manager.dict()

        for plot in plots:
            pool.apply_async(
                run_plot,
                args=(
                    plot,
                    validation_files,
                    proc_dir,
                    name_path_dict,
                ),
            )

        pool.close()
        pool.join()

        print("")
        print(" producing report ... ")
        print("")

        ## produce latex report
        report_name = cfg.name
        report_latex = LatexReport(report_name, sample_dir)
        for section_title, subsections in cfg.report.items():
            report_latex.section(section_title)
            for subsection_title, subsubsections in subsections.items():
                report_latex.subsection(subsection_title)
                for title, plots in subsubsections.items():
                    report_latex.begin_frame(title)
                    figs = []
                    for plot in plots:
                        figs.append(name_path_dict[plot.name])
                    report_latex.add_figures(figs)

        report_latex.compile()

        print(" ")
        print("{}.pdf report file has been produced".format(report_name))
        print(" ")

    ## do some cleaning
    os.system("rm -rf .matplotlib")


# _____________________________________________________________________________________________________________
def run_cmd(cmd):
    os.system(cmd)


# _____________________________________________________________________________________________________________
def run_plot(plot, validation_files, dir, name_path_dict):
    name = plot.name
    path = plot.plot(validation_files, dir)
    name_path_dict[name] = path


# _______________________________________________________________________________________
if __name__ == "__main__":
    main()
