import parsl
import os
from parsl.app.app import python_app, bash_app
from parsl.data_provider.files import File
from path import Path
from parslpw import pwconfig, pwargs

parsl.clear()
parsl.load(pwconfig)

def replace_file(sub_dict,input_file,output_file=None):
    with open(input_file, 'r') as f:
        s = f.read()
        for key in sub_dict.keys():
            s = s.replace(key, sub_dict[key])

    if output_file is None:
        print(s)
    else:
        with open(output_file,'w') as f:
            f.write(s)
        return output_file

@bash_app
def run_pesmd(inputs=[], outputs=[], stderr="run.err",pesmd_script="/home/wpc252/parsl/infrMetad/pesmd.sh" ):
    import os
    #input 0 is pesmd file
    pesmd_input_file = inputs[0].filepath
    pesmd_file_dir = os.path.dirname(pesmd_input_file)
    pesmd_file_name = os.path.basename(pesmd_input_file)
    return "cd %s; %s %s"%(pesmd_file_dir,pesmd_script,pesmd_file_name)

@bash_app
def remove_files(rundir):
    import os
    path = os.path.join(rundir, 'runinfo')
    if path == 'runinfo':
        return 0
    else:
        return f"rm -rf {path}"


def perform_analysis(source_dir, force_list, file_template, num_runs, model_name, output):
    id_run = 1
    metad_data = metaD.read_metad(source_dir, force_list, file_template=file_template, num_runs=20)
    mean_accelerations, escape_times = metaD.get_escape_times(metad_data)
    tau_list = []
    p_values = []
    pfile = open(F"{output}/pvalues{id_run}.dat", 'w')

    for force in escape_times:
        time_list = escape_times[force]
        print(force)
        metaD.plot_kde_hist_cumulative(time_list, id_run, out_dir=output,force=force)
        cdfs, data = metaD.compute_cdfs(time_list)
        time_domain = cdfs[0]
        ecdf = cdfs[1]
        tcdf = cdfs[2]
        tau = data[0]
        tau_list.append(tau)
        metaD.plot_cdfs(x=time_domain, y=ecdf, tcdf=tcdf, id=id_run, force=force, model_name=model_name, out_dir=output)
        s, p = metaD.do_KS2test(time_list, tau)
        p_values.append(p)
        pfile.write(F"tau: {tau} pvalue: {p} force: {force} \n")
    pfile.close()
    error_bars = metaD.calculate_error_bars(tau_list, num_runs)
    metaD.plot_rates(force_list, tau_list, model_name, output, id_run, error_bars)
    metaD.plot_escape_time(force_list, tau_list, model_name, output, id_run, error_bars)

if __name__ == "__main__":
    import sys
    import AnalysisMetad as metaD

    n_repeats = int(pwargs.nruns)
    outdir = pwargs.outdir
    fstop = float(pwargs.fstop)*-1
    fstart = float(pwargs.fstart)
    nforce = int(pwargs.nforce)
    force_list = metaD.np.linspace(fstop,fstart,nforce)
    run_dir = metaD.os.getcwd()
    source_dir = metaD.os.path.join('./', outdir)
    output = metaD.os.path.join(source_dir, "figures")
    calc_type = bool(pwargs.calc)
    model_type = bool(pwargs.model)
    calc_type_str = 'infr'
    if calc_type == False:
        calc_type_str ='fes'
        n_repeats=1

    print('output at: ', source_dir)
    print('figures at: ', output)
    print('current directory: ', run_dir)

    plumed_template=F'./dw{calc_type_str}.plumed.dat'
    input_template=F'./dw{calc_type_str}.pesmd.input'

    plumed_label = metaD.os.path.basename(plumed_template).split('.')[0]
    result_list = []

    print("......Going in main loop to set up......")
    for force in force_list:
        for seed in range(1,n_repeats+1):
            #put each job in separate directories, because pesmd writes other random files
            output_directory=metaD.os.path.join(source_dir,"%s_KJp_F%3.2f"%(plumed_label,force),"%i"%seed)

            metaD.os.makedirs(output_directory, exist_ok=True)

            out_label = metaD.os.path.basename(metaD.os.path.dirname(output_directory))+"_s%i"%seed
            output_prefix = metaD.os.path.join(output_directory,out_label)
            sub_dict = { "_SEED_": "%i"%seed, "_OUTPREFIX_": out_label, "_FORCE_": "%3.2f"%force,}
            if model_type==False:
                sub_dict = { "_SEED_": "%i"%seed, "_OUTPREFIX_": out_label, "_FORCE_": "%3.2f"%force,}
            output_dir = Path(output_directory)
            plumed_file = Path( replace_file(sub_dict, plumed_template,output_prefix+".plumed.dat") )
            pesmd_input_file = Path( replace_file(sub_dict,input_template,output_prefix+".pesmd.input") )

            r = run_pesmd(inputs=[pesmd_input_file, plumed_file], outputs=[output_dir])
            print("queued %d %3.2f"%(seed, force))
            result_list.append(r)
    print("......Set up done......")
    [r.result() for r in result_list]
    print("......Runs finished......")


    #  Analysis is executed below
    plt_kde_force = []
    plt_trj_time = []
    metaD.os.makedirs(output, exist_ok=True)
    distance_traj = metaD.read_trajectory(source_dir, force_list, num_runs=n_repeats)
    metad_data = metaD.read_metad(source_dir, force_list, num_runs=n_repeats)
    mean_accelerations, escape_times = metaD.get_escape_times(metad_data)

    print('......first plots.......')
    for i,force in enumerate(distance_traj):
        runs = distance_traj[force]
        plt_a = metaD.plot_distance_distribution(runs, force, output)
        plt_b = metaD.plot_distance_trajectory(runs, force, output)
        plt_kde_force.append(plt_a)
        plt_trj_time.append(plt_b)


    #  [r.result() for r in plt_kde]
    #  [r.result() for r in plt_trj_time]


    tau_list = []
    p_values = []
    plt_distribution = []
    plt_cdfs = []

    for force in escape_times:
        time_list = escape_times[force]
        plt_kde = metaD.plot_kde_hist_cumulative(time_list, out_dir=output, force=force)
        plt_distribution.append(plt_kde)
        cdfs, data = metaD.compute_cdfs(time_list)
        time_domain = cdfs[0]
        ecdf = cdfs[1]
        tcdf = cdfs[2]
        tau = data[0]
       # print(tau)
        tau_list.append(tau)
        plt_cdf = metaD.plot_cdfs(x=time_domain, y=ecdf, tcdf=tcdf, force=force, model_name='pesmd_d', out_dir=output)
        plt_cdfs.append(plt_cdf)
        s, p = metaD.do_KS2test(time_list, tau)
        p_values.append(p)
        # print(p, "fromDriver")

    #[r.result() for r in plt_distribution]
    #[r.result() for r in plt_cdfs]

    metaD.plot_rates(-1*force_list, tau_list, 'pesmd_d', output, True)
    metaD.plot_rates(-1*force_list, tau_list, 'pesmd_d', output)
    metaD.plot_escape_time(-1*force_list, tau_list, 'pesmd_d', output, True)
    metaD.plot_escape_time(-1*force_list, tau_list, 'pesmd_d', output)

    rows = ['in:force, img:traj, img:dist, img:fit, img:kde, out:tau']
    for i,force in enumerate(escape_times):
        tau_i = tau_list[i]
        plt_b = os.path.join(run_dir,plt_trj_time[i])
        plt_a = os.path.join(run_dir,plt_kde_force[i])
        plt_c = os.path.join(run_dir,output,F"F{force}_pesmd_d.png")
        plt_d = os.path.join(run_dir,output,F"kde_hist_cumulative{force}.png")
        rows.append(F"{-1*float(force)},{plt_b},{plt_a},{plt_c},{plt_d},{tau_i}")

    with open(pwargs.outcsv,'w') as file:
        file.write("\n".join(rows))
    with open(pwargs.outhtml, 'w') as html:
        html.write("""\
        <html style="overflow-y:hidden;background:white">
    <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px" href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>
    <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>
</html>
        """.format(pwargs.outcsv,'tau',pwargs.outcsv, 'tau'))

    sys.exit()
