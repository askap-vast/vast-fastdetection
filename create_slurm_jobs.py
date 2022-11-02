# from vastfast.setting import N_BEAM, MEM, N_CPU
import configparser

## ========= read in config file =========
config = configparser.ConfigParser()		
config.read("run_vastfast.ini")

n_beam = int(config["INSTRUMENT"]["n_beam"])

work_dir = config["PATH"]["work_dir"]
setup_file = config["EXE"]["setup_file"]
python_file = config["EXE"]["python_file"]

mem = config["RESOURCE"]["mem"]
n_cpu = config["RESOURCE"]["n_cpu"]
time = config["RESOURCE"]["time"]



## ========= write job script of a single beam ========= 
def write_job_single(beam):
    filename = "vastfast_beam{:02}.sh".format(beam)
    out = open(filename, 'w')
    out.write("#!/bin/bash \n")
    out.write("# \n")
    out.write("#SBATCH --job-name=vastfast_beam{:02} \n".format(beam))
    out.write("#SBATCH --output=vastfast_beam{:02}.log \n\n".format(beam))
    out.write("#SBATCH --time={} \n".format(time))
    out.write("#SBATCH --mem={} \n".format(mem))
    out.write("#SBATCH --cpus-per-task={} \n\n".format(n_cpu))
    out.write("export MPLBACKEND=agg \n")
    out.write("source {} \n".format(setup_file))
    out.write("cd {} \n\n".format(work_dir))
    out.write("python {} {} \n\n".format(python_file, beam))
    out.close()

## ========= create list of job scripts =========
def create_jobs():
    filename = "submit_vastfast_jobs.sh"
    out = open(filename, 'w')
    out.write("#!/bin/bash \n\n")
    for i in range(n_beam):
        write_job_single(i)
        job_file = "vastfast_beam{:02}.sh".format(i)
        out.write("sbatch " + job_file + "\n")
    out.close()




if __name__ == "__main__":
    create_jobs()
    