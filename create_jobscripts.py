N_BEAM = 36
MEM = 10
N_CPU = 4

module_file = "/fred/oz999/tiger/SS2022B-TMurphy/ozstar-dev.sh"
work_dir = "/fred/oz999/tiger/SS2022B-TMurphy/tmurphy_2022b"
python_file = "run_all_single.py"

def write_job_single(beam):
    filename = "vastfast_beam{:02}.sh".format(beam)
    out = open(filename, 'w')
    out.write("#!/bin/bash \n")
    out.write("# \n")
    out.write("#SBATCH --job-name=vastfast_beam{:02} \n".format(beam))
    out.write("#SBATCH --output=vastfast_beam{:02}.log \n\n".format(beam))
    out.write("#SBATCH --time=02:00:00\n")
    out.write("#SBATCH --mem={}g \n".format(MEM))
    out.write("#SBATCH --cpus-per-task={}\n\n".format(N_CPU))
    out.write("source " + module_file + "\n")
    out.write("cd " + work_dir + "\n\n")
    out.write("python " + "run_all_single.py " + str(beam) + "\n\n")
    out.close()

def create_jobs():
    filename = "submit_vastfast_jobs.sh"
    out = open(filename, 'w')
    out.write("#!/bin/bash \n\n")
    for i in range(N_BEAM):
        write_job_single(i)
        job_file = "vastfast_beam{:02}.sh".format(i)
        out.write("sbatch " + job_file + "\n")
    out.close()




if __name__ == "__main__":
    create_jobs()
    