import os
import psutil


def get_minenv(envstr, base):
    # retrieve environment variable, return minimum relative to base
    # cast comparison is to catch instance when envstr does not exist, returns None
    return min(int(os.getenv(envstr) or base), base)


def get_resources():
    # NOTE: memory should be reported in bytes

    # node totals are the MAXIMUM
    avail_cpu = psutil.cpu_count()
    avail_mem = psutil.virtual_memory()[0] / 1e6 # in MB since SLURM is in MB

    # check various environment variables
    avail_cpu = get_minenv("AVAIL_CPU", avail_cpu)
    avail_cpu = get_minenv("SLURM_CPUS_ON_NODE", avail_cpu)
    avail_mem = get_minenv("AVAIL_MEM", avail_mem)
    avail_mem = get_minenv("SLURM_MEM_PER_NODE", avail_mem)

    assert avail_cpu >= 1
    assert avail_mem > 0

    # convert back to bytes
    avail_mem *= 1e6

    return avail_cpu, avail_mem