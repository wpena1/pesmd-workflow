import parsl
from parsl.app.app import python_app
import numpy as np
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.optimize import curve_fit
from scipy.stats import gamma
from scipy import stats

def read_trajectory(src_dir, force_list,
                    file_template="dwinfr_KJp_F_FORCE_/RUN/dwinfr_KJp_F_FORCE__sRUN.temp1.0.colvar.out",
                    num_runs=20, f_stop=-2.00, f_start=0.0, n_force=11):
    data = {}
    # force_list = np.linspace(f_stop, f_start, n_force)
    for i in force_list:
        key = "%3.2f" % i
        data[key] = []
    for k in data:
        file_name = file_template.replace('_FORCE_', k)
        runs = []
        for i in range(1, num_runs + 1):
            final = file_name.replace('RUN', str(i))
            fullpath = os.path.join(src_dir, final)
            temp = np.genfromtxt(fullpath, comments='#', dtype=float)
            runs.append(temp)
        data[k] = np.array(runs, dtype=object)
    return data


#  @python_app
def plot_distance_distribution(runs, force, out_dir):

    import matplotlib
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(6, 5))
    for i,run in enumerate(runs):
        h = sns.kdeplot(run[:, 1], label=f'run: {i+1}')
    plt.xlabel(r"$Distance\ \AA$")
    plt.legend(ncol=3, loc='best')
    out_path = os.path.join(out_dir, f"rDistribution_{force}.png")
    plt.savefig(out_path)
    plt.close()
    return out_path


#  @python_app
def plot_distance_trajectory(runs, force, out_dir):

    import matplotlib
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure(figsize=(6, 5))
    for i,run in enumerate(runs):
        lab = np.abs(np.float64(force))
        plt.plot(run[:, 0], run[:, 1], label=F'run: {i+1}')
        plt.xlabel(r"$Time$")
        plt.ylabel(r"$Distance\ \AA$")
    plt.legend(ncol=3, loc='best')
    out_path = os.path.join(out_dir, f"distTrajectory_{force}.png")
    plt.savefig(out_path)
    plt.close()
    return out_path


def read_metad(src_dir, force_list, file_template="dwinfr_KJp_F_FORCE_/RUN/dwinfr_KJp_F_FORCE__sRUN.metad.out",
               num_runs=20):
    acc = {}
    for i in force_list:
        key = "%3.2f" % i
        acc[key] = []
    for k in acc:
        file_name = file_template.replace('_FORCE_', k)
        info = []
        for i in range(1, num_runs + 1):
            final = file_name.replace('RUN', str(i))
            fullpath = os.path.join(src_dir, final)
            with open(fullpath) as f:
                for line in f:
                    pass
            info.append([np.float64(f) for f in line.split()])
        acc[k] = np.array(info)
    return acc


def get_escape_times(data):
    mean_product = []
    escape_times = {}
    for k in data:
        runs = data[k]
        product = []
        last = []
        alpha = []
        for run in runs:
            scaled = run[0] * run[2] * 1e-6
            product.append(scaled)
            alpha.append(run[2])
            last.append(run[0])
            # print("%s %f %f %f" % (k, run[0], run[2], scaled))
            # out_file.write("%s %f %f %f \n" % (k, run[0], run[2], scaled))
        mean_product.append(np.mean(product))
        escape_times[k] = product
    return mean_product, escape_times


#  @python_app
def plot_kde_hist_cumulative(times_vector, out_dir, force):
    import matplotlib
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.figure(figsize=(5, 4))
    sns.kdeplot(times_vector, common_norm=True, cumulative=True)
    a = plt.hist(times_vector, density=True, bins=20, cumulative=True)
    plt.xscale('log')
    plt.savefig(f"{out_dir}/kde_hist_cumulative{force}.png")
    plt.close()
    return 1




# Function to calculate the TCDF with constants rate
def func(x, a):
    return 1 - np.exp(-a * x)


# Function to fit line
def line(x, a, b):
    return a * x + b


def get_ecdf(time_list, min_time, max_time, n_bins):
    time_domain = np.logspace(np.log10(min_time), np.log10(max_time), n_bins)
    N = len(time_list)
    a, b = np.histogram(time_list, time_domain)
    # mod = np.append(a,[0])
    ECDF = np.cumsum(a) / N
    # new = np.append(ECDF,[ECDF[-1]])
    if len(time_domain[:-1]) == len(ECDF):
        # print(len(mod), len(a[:-1]), len(ECDF))
        return time_domain[:-1], ECDF
    else:
        return "error"


def compute_cdfs(time_list):
    min_t = np.min(time_list)/10  # maybe another way of getting these
    max_t = np.max(time_list)*10
    nbins = 500
    x, y = get_ecdf(time_list, min_t, max_t, nbins)
    test_mean = np.mean(time_list)
    test_sigma = np.std(time_list, ddof=1)
    test_m = np.median(time_list)
    means = test_mean
    medians = test_m
    uncertainty = test_sigma / np.sqrt(len(time_list))
    stdvs = test_sigma
    mean_sigma_ratio = test_mean / test_sigma
    log2_median_ratio = np.log(2) * test_mean / test_m
    guess = 1.0 / test_mean
    # print(np.mean(time_list),min_t,max_t)
    pars, cov = curve_fit(f=func, xdata=x, ydata=y, p0=[guess], maxfev=1000000)
    stdevs = np.sqrt(np.diag(cov))
    f_stdvs = stdevs[0]
    tau = np.power(pars[0], -1)
    tau_mean_ratio = (tau / test_mean)
    tcdf = gamma.cdf(x, 1, scale=tau)
    cdfs = [x, y, tcdf]
    data = [tau, stdvs, f_stdvs, means, medians,
            uncertainty, mean_sigma_ratio, log2_median_ratio, tau_mean_ratio]
    return cdfs, data


#  @python_app
def plot_cdfs(x, y, tcdf, force, model_name, out_dir):
    import matplotlib
    matplotlib.rcParams['savefig.dpi'] = 300
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.step(x, y, c='k', label=F'ECDF: {np.abs(float(force))} pN {model_name}')
    ax.plot(x, tcdf, label=F'TCDF: {np.abs(float(force))} pN {model_name}', c='r')
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
    ax.set_xlim(min(x), max(x))
    ax.tick_params(axis='both', labelsize=14, width=2)
    plt.legend(fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(out_dir,f'F{force}_{model_name}.png')
    plt.savefig(out_path)
    plt.close()
    return 1


def do_KS2test(time_list, tau):
    size = np.int64(len(time_list) * 1e6)
    rvs1 = gamma.rvs(1, scale=tau, size=size)  # loc=min(time_list))
    ks_stat, p = stats.ks_2samp(rvs1, time_list)
    return ks_stat, p


def plot_escape_time_log(force_list, tau_list, model_name, out_dir):
    fig, ax = plt.subplots(figsize=(6, 8))
    lab = F"{model_name}"
    ax.plot(force_list, tau_list, c='k', marker='.', label=lab)
    ax.tick_params(axis='both', labelsize=14, width=2)
    ax.legend(fontsize=13, loc='lower left')
    ax.set_ylabel(r"${\tau}\ ({\mu}s)$", fontsize=14)
    ax.set_yscale('log')
    plt.xlabel(r"$Force\ (pN)$", fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=.92)
    out_path = os.path.join(out_dir,f'{model_name}_logOfFittedTaus.png')
    plt.savefig(out_path)
    plt.close()
    return 1


def plot_escape_time(force_list, tau_list, model_name, out_dir, log=False):
    force_list = np.array(force_list) * -1
    tau_list = np.array(tau_list)
    if log:
        plot_escape_time_log(force_list, tau_list, model_name, out_dir)
    else:
        fig, ax = plt.subplots(figsize=(6, 8))
        lab = F"{model_name}"
        ax.plot(force_list, tau_list, c='k', marker='.', label=lab)
        ax.tick_params(axis='both', labelsize=14, width=2)
        ax.legend(fontsize=13, loc='lower left')
        ax.set_ylabel(r"${\tau}\ ({\mu}s)$", fontsize=14)
        plt.xlabel(r"$Force\ (pN)$", fontsize=14)
        fig.tight_layout()
        fig.subplots_adjust(top=.92)
        out_path = os.path.join(out_dir, f'{model_name}_FittedTaus.png')
        plt.savefig(out_path)
        plt.close()
    return 1


def plot_rates_log(force_list, tau_list, model_name, out_dir):
    fig, ax = plt.subplots(figsize=(6, 8))
    lab = F"{model_name}"
    ax.plot(force_list, 1 / tau_list, c='k', marker='.', label=lab)
    ax.tick_params(axis='both', labelsize=14, width=2)
    ax.legend(fontsize=13, loc='lower left')
    ax.set_ylabel(r"$rate\ ({\mu}s^{-1})$", fontsize=14)
    ax.set_yscale('log')
    plt.xlabel(r"$Force\ (pN)$", fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=.92)
    out_path = os.path.join(out_dir, f'{model_name}_logOfRatesvsForce.png')
    plt.savefig(out_path)
    plt.close()
    return 1


def plot_rates(force_list, tau_list, model_name, out_dir, log=False):
    force_list = np.array(force_list) * -1
    tau_list = np.array(tau_list)
    if log:
        plot_rates_log(force_list, tau_list, model_name, out_dir)
    else:
        fig, ax = plt.subplots(figsize=(6, 8))
        lab = F"{model_name}"
        ax.plot(force_list, 1 / tau_list, c='k', marker='.', label=lab)
        ax.tick_params(axis='both', labelsize=14, width=2)
        ax.legend(fontsize=13, loc='lower left')
        ax.set_ylabel(r"$rate\ ({\mu}s^{-1})$", fontsize=14)
        plt.xlabel(r"$Force\ (pN)$", fontsize=14)
        fig.tight_layout()
        fig.subplots_adjust(top=.92)
        out_path = os.path.join(out_dir, f'{model_name}_ratesvsForce.png')
        plt.savefig(out_path)
        plt.close()
    return 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parameters for analysis')
    parser.add_argument('--srcdir', type=str, default='/scratch/wpc252/parsl', help='source directory')
    args = parser.parse_args()
    source_dir = args.srcdir
    print(source_dir)
    time_list = [389267.8260916,
                 280615.1612349,
                 352991.755287,
                 254401.5943329,
                 290813.855694,
                 127531.5566096,
                 225731.5920222,
                 350577.3045656,
                 269457.00861759996,
                 314554.8179661,
                 472361.6314188,
                 299550.19261675,
                 346090.13722035,
                 313239.37392205006,
                 218508.85789239997,
                 434531.22000145,
                 224990.1612792,
                 156855.9787566,
                 168118.58454685,
                 301677.0178401]
    plot_kde_hist_cumulative(time_list, out_dir="./", force=0.0)
    cdfs, data = compute_cdfs(time_list)
    time_domain = cdfs[0]
    ecdf = cdfs[1]
    tcdf = cdfs[2]
    tau = data[0]
    print(tau)
    plot_cdfs(x=time_domain, y=ecdf, tcdf=tcdf, force=0.0, model_name='pesmd', out_dir='./')
    s, p = do_KS2test(time_list, tau)
    print(p)

