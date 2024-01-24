import parsl
from parsl.app.app import python_app, bash_app
print(parsl.__version__, flush = True)

# import parsl utils stuff below
import parsl_utils
from parsl_utils.config import config, resource_labels, form_inputs #, executor_dict
from parsl_utils.data_provider import PWFile

print("MAIN.py:...Configuring Parsl...")
parsl.load(config)
print("MAIN.py:...Parsl config loaded...")

def replace_file(sub_dict, input_file, output_file=None):
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

@parsl_utils.parsl_wrappers.log_app
@bash_app
def run_pesmd(inputs=[], outputs=[], stdout="run.out", stderr="run.err", pesmd_script="./pesmd.sh", calctype='FES', sum_hills=1):
    import os
    #input 0 is pesmd file
    if calctype == 'RATE': #  rate case
        calctype = 1
    else: #  fes case
        calctype = 0
    pesmd_script=os.path.basename(pesmd_script)
    pesmd_input_file = inputs[0].filepath
    outputprefix = os.path.basename(pesmd_script).split("_")[0]
    pesmd_file_dir = os.path.dirname(pesmd_input_file)
    pesmd_file_name = os.path.basename(pesmd_input_file)
    return "cd %s;bash %s %s %s %s %s"%(pesmd_file_dir, pesmd_script, pesmd_file_name, outputprefix, calctype, sum_hills)

@parsl_utils.parsl_wrappers.log_app
@bash_app
def remove_files(rundir):
    import os
    path = os.path.join(rundir, 'runinfo')
    if path == 'runinfo':
        return 0
    else:
        return f"rm -rf {path}"


def is_str_numeric(txt, type):
    if type == 'float':
        try:
            float(txt)
            return True
        except ValueError:
            return False
    elif type == 'int':
        try:
            int(txt)
            return True
        except ValueError:
            return False
    else:
        print("func: is_str_numeric: type not supported")
        return 0

if __name__ == "__main__":
    import os
    import sys
    import numpy as np

    print(F"MAIN.py:...starting workflow setup...\n Resource list: {resource_labels}")
    
    model_type = str(form_inputs['model_inputs']['modeltype'])
    calc_type = str(form_inputs['model_inputs']['calctype'])

    if model_type == 'slip':
        fstop = form_inputs['model_inputs']['fstop']
        fstart = form_inputs['model_inputs']['fstart']
        nforce = form_inputs['model_inputs']['nforce']
        potential_rgx = form_inputs['model_inputs']['func']
    elif model_type == 'catch':
        fstop = form_inputs['model_inputs']['fstopc']
        fstart = form_inputs['model_inputs']['fstartc']
        nforce = form_inputs['model_inputs']['nforcec']
        potential_rgx = form_inputs['model_inputs']['funcc']
    else:
        print("MAIN.py ERROR: model type is not slip or catch")
        sys.exit()
    
    if is_str_numeric(fstop, 'float'):
        fstop = float(fstop)*-1
    else:
        if model_type == 'slip':
            fstop = 2.0*-1
        elif model_type == 'catch':
            fstop = 8.0*-1
    
    if is_str_numeric(fstart, 'float'):
        fstart = float(fstart)
    else:
        fstart = 0.0
     
    if is_str_numeric(nforce, 'int'):
        nforce = int(nforce)
    else:
        nforce = 11
    
    force_list = np.linspace(fstop,fstart,nforce)

    if calc_type == 'FES':
        simlength = form_inputs['model_inputs']['simlength']
        n_repeats = form_inputs['model_inputs']['nruns']
        outdir = form_inputs['model_inputs']['outdir']
    elif calc_type =='RATE':
        simlength = form_inputs['model_inputs']['simlengthr']
        n_repeats = form_inputs['model_inputs']['nrunsr']
        outdir = form_inputs['model_inputs']['outdirr']
    else:
        print("MAIN.py ERROR: calculation type is not FES or Rate")
        sys.exit()

    if is_str_numeric(simlength, 'float'):
        simlength = float(simlength)
    else:
        if calc_type == 'FES':
            simlength = 15.0
        else:
            simlength = 7.5

    if is_str_numeric(n_repeats, 'int'):
        n_repeats = int(n_repeats)
    else:
        if calc_type == 'FES':
            n_repeats = 1
        else:
            n_repeats = 14


    run_dir = os.getcwd()
    # for key in executor_dict:
    #    remote_dir = executor_dict[key]["resource"]["jobdir"]
    remote_dir = config.executors[0].working_dir + F"/{outdir}"
    
    source_dir = os.path.join('.', outdir)
    output = os.path.join(source_dir, "figures")
    n_steps = int(1e6*simlength/2)
    kbT = 2.249 # 300K in KJ/mol
    rid="r1"

    print('MAIN.py: output at: ', source_dir)
    print('MAIN.py: figures at: ', output)
    print('MAIN.py: current directory: ', run_dir)

    if model_type == 'slip': # Slip case
        print("MAIN.py: Slip model requested")
        if calc_type == 'FES': # FES calculation
            calc_type_str = 'fes'
            n_repeats=1
            output = os.path.join(source_dir, "slip", "fes", "figures" )
            output_data = os.path.join(source_dir, "slip", "fes", "data" )
            print("MAIN.py: FES calc requested", F'runs: {n_repeats}')
        else: # Rates calculation
            calc_type_str = 'infr'
            if n_repeats == 1:
                n_repeats = 20
            output = os.path.join(source_dir, "slip", "rates", "figures" )
            output_data = os.path.join(source_dir, "slip", "rates", "data" )
            print("MAIN.py: Rates calc requested", F'runs: {n_repeats}')
        plumed_template=F'./plumed_inputs/dw{calc_type_str}.plumed.dat'
        input_template=F'./plumed_inputs/dw{calc_type_str}.pesmd.input'
        model_name='slip'

    else: # Catch case
        print("MAIN.py: Catch model requested")
        if calc_type == 'FES':
            calc_type_str = 'fes'
            n_repeats=1
            output = os.path.join(source_dir, "catch", "fes", "figures" )
            output_data = os.path.join(source_dir, "catch", "fes", "data" )
            print("MAIN.py: FES calc requested", F'runs: {n_repeats}')
        else:
            lenfl = len(force_list)
            deltax = (22.65 - 13.5)/(lenfl-1)
            deltay = (8.65 - 7.5)/(lenfl-1)
            mins=13.25
            maxs=14.0
            minys=7.25
            maxys=8.0
            force_list_c = np.array([round(k,2) for k in np.flip(force_list)])*-1
            ln_c = len(np.where(force_list_c<3.2)[0])
            init_deltax_c = (7.5 - 3.5)/ln_c
            init_deltay_c = (-1.7 + 2.5)/ln_c

            ln_s = len(np.where(force_list_c > 3.6)[0])
            init_deltax_s = (15.0 - 12.8)/ln_s
            c_count = 0
            s_count = 0

            calc_type_str = 'infr'
            if n_repeats == 1:
                n_repeats = 20
            output = os.path.join(source_dir, "catch", "rates", "figures" )
            output_data = os.path.join(source_dir, "catch", "rates", "data" )
            positions = np.genfromtxt(F"./plumed_inputs/positions.dat", dtype=float)
            lenpos = len(positions)
            print("MAIN.py: Rates calc requested", F'runs: {n_repeats}')

        plumed_template=F'./plumed_inputs/tw{calc_type_str}.plumed.dat'
        input_template=F'./plumed_inputs/tw{calc_type_str}.pesmd.input'
        model_name='catch'

    plumed_label = os.path.basename(plumed_template).split('.')[0]
    result_list = []

    print("MAIN:.......Going in main loop to set up......")
    sumhills = 0
    if model_type == 'catch' and calc_type == 'FES': # catch, fes
        sumhills = 2
    elif model_type == 'slip' and calc_type == 'FES': # slip, fes
        sumhills = 1

    for i,force in enumerate(np.flip(force_list)):
        if model_type == 'catch' and calc_type == 'RATE':
            if lenfl == lenpos:
                pair=positions[i]
                init_x=pair[1]
                init_y=pair[2]
                stopx = pair[5]
                stopy = pair[6]

                min_x = stopx - 0.25
                max_x = stopx + 0.5
                min_y = stopy - 0.25
                max_y = stopy + 0.5
            else:
                if (force*-1) <= 3.2:
                    init_x = 3.5 + (c_count*init_deltax_c)
                    init_y = -2.5 + (c_count*init_deltay_c)
                    c_count += 1
                else:
                    init_x = 12.8 + (s_count*init_deltax_s)
                    init_y = -0.5
                    s_count += 1
                min_x = mins + (i*deltax)
                max_x = maxs + (i*deltax)
                min_y = minys + (i*deltay)
                max_y = maxys + (i*deltay)

        for seed in range(1,n_repeats+1):
            # put each job in separate directories, because pesmd writes other random files
            output_directory = os.path.join(source_dir,"%s_KJp_F%3.2f"%(plumed_label,force),"%i"%seed)
            output_dir_remote = "%s_KJp_F%3.2f"%(plumed_label,force) + "/%i"%seed
            output_dir = PWFile(url='file://usercontainer/' + output_directory, local_path=remote_dir+output_dir_remote)
            # os.makedirs(output_dir, exist_ok=True)
            out_label = f"{seed}_pesmd"
            output_prefix = os.path.join(output_directory,out_label)
            sub_dict = { "_SEED_": "%i"%seed, "_OUTPREFIX_": out_label, "_FORCE_": "%3.2f"%force, "_NSTEP_": "%d"%n_steps}
            if model_type == 'catch' and calc_type == 'RATE':
                sub_dict = { "_SEED_": "%i"%seed, "_OUTPREFIX_": out_label, "_FORCE_": "%3.2f"%force, "_INX_": "%3.2f"%init_x,\
                 "_INY_": "%3.2f"%init_y, "_NSTEP_": "%d"%n_steps, "_MIN_": "%f"%min_x, "_MAX_":"%f"%max_x, "_MINY_":"%f"%min_y,"_MAXY_":"%f"%max_y}

           
            plumed_file_path = replace_file(sub_dict, plumed_template,output_prefix+".plumed.dat")
            plumed_file = PWFile(url='file://usercontainer/'+ plumed_file_path, local_path=remote_dir+'/src')
            plumed_input_path = replace_file(sub_dict,input_template,output_prefix+".pesmd.input")
            pesmd_input_file = PWFile(url='file://usercontainer/'+ plumed_input_path, local_path=remote_dir+'/src')
            pesmd_script_path = replace_file({}, "./utils/pesmd.sh", output_prefix+".pesmd.sh")
            pesmd_script = PWFile('', pesmd_script_path)
            
            r = run_pesmd(inputs=[pesmd_input_file, plumed_file], outputs=[output_dir], pesmd_script=pesmd_script, calctype=calc_type, sum_hills=sumhills)
            print("MAIN.py: queued %d %3.2f"%(seed, force))
            result_list.append(r)


    print("MAIN.py: ......lenght of runs......", len(result_list))
    print("MAIN.py: ......Set up done......")
    [r.result() for r in result_list]
    print("MAIN.py: ......Runs finished......")
    # sys.exit()

    print("MAIN.py: ......Excuting Analysis......")
    import utils.AnalysisMetad2 as metaD

    os.makedirs(output, exist_ok=True)
    os.makedirs(output_data, exist_ok=True)

    if model_type == 'slip': # Slip case
        plt_kde_force = []
        plt_trj_time = []
        print(os.getcwd(), "current dir")
        ftemplate = F"dw{calc_type_str}_KJp_F_FORCE_/RUN/RUN_pesmd.temp_1.0.colvar.out"
        distance_traj = metaD.read_trajectory(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
        print('MAIN.py: ......first plots.......')
        for i,force in enumerate(distance_traj):
            runs = distance_traj[force]
            plt_a = metaD.plot_distance_distribution(runs, rid, force, output)
            plt_b = metaD.plot_distance_trajectory(runs, rid, force, output)
            plt_kde_force.append(plt_a)
            plt_trj_time.append(plt_b)

        if calc_type == 'FES': # Fes calc
            fes_plots = []
            ftemplate = F"dw{calc_type_str}_KJp_F_FORCE_/RUN/RUN.fes.dat"
            fes_data = metaD.read_trajectory(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
            dim=1
            for i, force in enumerate(fes_data):
                runs = fes_data[force]
                fesplot = metaD.plot_fes(runs, force, output, dim)
                fes_plots.append(fesplot)

            rows = ['in:force, img:traj, img:dist, img:fes']
            for i,force in enumerate(fes_data):
                plt_b = os.path.join(run_dir,plt_trj_time[i])
                plt_a = os.path.join(run_dir,plt_kde_force[i])
                plt_c = os.path.join(run_dir,fes_plots[i])

                rows.append(F"{-1*float(force)},{plt_b},{plt_a},{plt_c}")

            with open(F'{output_data}/{model_name}.csv','w') as csvfile:
                csvfile.write("\n".join(rows))
            with open(F'{output_data}/{model_name}.html', 'w') as html:
                csvpath=metaD.os.path.join(run_dir,output_data,F'{model_name}.csv')
                html.write("""\
                <html style="overflow-y:hidden;background:white">\
                <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px"\
                href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>\
                <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>\
                </html>""".format(csvpath,'force', csvpath, 'force'))


        else: #Rates Calc
            ftemplate = F"dw{calc_type_str}_KJp_F_FORCE_/RUN/RUN_pesmd.metad.out"
            metad_data = metaD.read_metad(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
            mean_accelerations, escape_times = metaD.get_escape_times(metad_data, output_data, rid)

            tau_list = []
            p_values = []
            plt_distribution = []
            plt_cdfs = []

            pfile = open(F"{output_data}/pvalues{rid}.dat", 'w')
            pfile.write("#force\t p-value\n")
            force_list = force_list*(1.66)*-1
            for force in escape_times:
                time_list = escape_times[force]
                escape_file = open(F"{output_data}/escape_times_{force}{rid}.dat", 'w')
                escape_file.write("#run\t escape_time\n")
                print(force)
                for i in range(len(time_list)):
                    escape_file.write("%d\t %f\n"%(i+1, time_list[i]))
                escape_file.close()
                plt_kde = metaD.plot_kde_hist_cumulative(time_list, rid, out_dir=output,force=force)
                plt_distribution.append(plt_kde)
                cdfs, data = metaD.compute_cdfs(time_list)
                time_domain = cdfs[0]
                ecdf = cdfs[1]
                tcdf = cdfs[2]
                tau = data[0]
                tau_list.append(tau)
                plt_cdf = metaD.plot_cdfs(x=time_domain, y=ecdf, tcdf=tcdf, id=rid, force=force, model_name=model_name, out_dir=output)
                plt_cdfs.append(plt_cdf)
                s, p = metaD.do_KS2test(time_list, tau)
                p_values.append(p)
                pfile.write("%f\t%f\n"%(float(force),p))
            pfile.close()

            rates = 1/np.array(tau_list)

            error_bar = metaD.calculate_error_bars(tau_list, n_repeats)

            taus_b = error_bar[:,0]
            bars = error_bar[:,1]
            rates_b = error_bar[:,2]
            r_bars = error_bar[:,3]

            taus_bars = open(F"{output_data}/taus_bars{rid}.dat", 'w')
            rates_bars = open(F"{output_data}/rates_bars{rid}.dat", 'w')
            rates_bars.write("#%s\t%s\t%s\t%s\n"%("forces","rates","rates_bootstrap", "bars"))
            taus_bars.write("#%s\t%s\t%s\t%s\n"%("forces","taus","taus_bootstrap", "bars"))

            for i in range(len(force_list)):
                taus_bars.write("%f\t%f\t%f\t%f\n"%(force_list[i],tau_list[i],taus_b[i],bars[i]))
                rates_bars.write("%f\t%f\t%f\t%f\n"%(force_list[i],rates[i],rates_b[i],r_bars[i]))
            rates_bars.close()
            taus_bars.close()

            rates_plot_path = metaD.plot_rates(force_list, tau_list, model_name, output, rid, error_bar)

            b_params = metaD.bell_params(force_list,tau_list)
            dx = round(-1*b_params[0]*kbT,3)
            k_0 = round(np.exp(b_params[1]),3)
            rows_rate = ['in:nruns,img:rates,out:k_0,out:dx']
            plt_g = os.path.join(run_dir, rates_plot_path)
            rows_rate.append(F'{n_repeats},{plt_g},{k_0},{dx}')

            with open(F'{output_data}/{model_name}_rates.csv','w') as csvfile:
                csvfile.write("\n".join(rows_rate))
            with open(F'{output_data}/{model_name}_rates.html', 'w') as html:
                csvpath=os.path.join(run_dir,output_data,F'{model_name}_rates.csv')
                html.write("""\
                <html style="overflow-y:hidden;background:white">\
                <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px"\
                href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>\
                <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>\
                </html>""".format(csvpath,'nruns', csvpath,'nruns'))


            rows = ['in:force, img:kde, img:traj, img:fit, img:dist, out:tau, out:pval']
            for i,force in enumerate(escape_times):
                plt_a = os.path.join(run_dir,plt_kde_force[i])
                plt_b = os.path.join(run_dir,plt_trj_time[i])
                plt_c = os.path.join(run_dir,plt_cdfs[i])
                plt_d = os.path.join(run_dir,plt_distribution[i])
                plt_e = tau_list[i]
                plt_f = p_values[i]

                rows.append(F"{-1*float(force)},{plt_a},{plt_b},{plt_c},{plt_d},{plt_e},{plt_f}")

            with open(F'{output_data}/{model_name}.csv','w') as csvfile:
                csvfile.write("\n".join(rows))
            with open(F'{output_data}/{model_name}.html', 'w') as html:
                csvpath=os.path.join(run_dir,output_data,F'{model_name}.csv')
                html.write("""\
                <html style="overflow-y:hidden;background:white">\
                <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px"\
                href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>\
                <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>\
                </html>""".format(csvpath,'force', csvpath, 'force'))


    else: # Catch case
        plt_kde_force = []
        plt_trj_time = []
        ftemplate = F"tw{calc_type_str}_KJp_F_FORCE_/RUN/RUN_pesmd.temp_1.0.colvar.out"
        distance_traj = metaD.read_trajectory(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
        print('MAIN.py: ......first plots.......')
        for i,force in enumerate(distance_traj):
            runs = distance_traj[force]
            plt_a = metaD.plot_distance_distribution(runs, rid, force, output)
            plt_b = metaD.plot_distance_trajectory(runs, rid, force, output)
            plt_kde_force.append(plt_a)
            plt_trj_time.append(plt_b)

        if calc_type == 'FES': # Fes Calc
            fes_plots2d = []
            fes_plotsx = []
            fes_plotsy = []

            ftemplate = F"tw{calc_type_str}_KJp_F_FORCE_/RUN/RUN.fes.dat"
            fes_data = metaD.read_trajectory(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
            dim=2
            for i, force in enumerate(fes_data):
                runs = fes_data[force]
                fesplot = metaD.plot_fes(runs, force, output, dim)
                fes_plots2d.append(fesplot)

            ftemplate = F"tw{calc_type_str}_KJp_F_FORCE_/RUN/RUN.fesx.dat"
            fes_data = metaD.read_trajectory(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
            dim=1
            for i, force in enumerate(fes_data):
                runs = fes_data[force]
                fesplot = metaD.plot_fes(runs, force, output, dim, 'x')
                fes_plotsx.append(fesplot)

            ftemplate = F"tw{calc_type_str}_KJp_F_FORCE_/RUN/RUN.fesy.dat"
            fes_data = metaD.read_trajectory(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
            dim=1
            for i, force in enumerate(fes_data):
                runs = fes_data[force]
                fesplot = metaD.plot_fes(runs, force, output, dim, 'y')
                fes_plotsy.append(fesplot)


            rows = ['in:force, img:traj, img:dist, img:fes, img:fesx, img:fesy']
            for i,force in enumerate(fes_data):
                plt_b = os.path.join(run_dir,plt_trj_time[i])
                plt_a = os.path.join(run_dir,plt_kde_force[i])
                plt_c = os.path.join(run_dir,fes_plots2d[i])
                plt_d = os.path.join(run_dir,fes_plotsx[i])
                plt_e = os.path.join(run_dir,fes_plotsy[i])

                rows.append(F"{-1*float(force)},{plt_b},{plt_a},{plt_c},{plt_d},{plt_e}")


            with open(F'{output_data}/{model_name}.csv','w') as csvfile:
                csvfile.write("\n".join(rows))
            with open(F'{output_data}/{model_name}.html', 'w') as html:
                csvpath=os.path.join(run_dir,output_data,F'{model_name}.csv')
                html.write("""\
                <html style="overflow-y:hidden;background:white">\
                <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px"\
                href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>\
                <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>\
                </html>""".format(csvpath, 'force', csvpath, 'force'))


        else: # Rate Calc
            ftemplate = F"tw{calc_type_str}_KJp_F_FORCE_/RUN/RUN_pesmd.metad.out"
            metad_data = metaD.read_metad(source_dir, force_list, file_template=ftemplate, num_runs=n_repeats)
            mean_accelerations, escape_times = metaD.get_escape_times(metad_data, output_data, rid)
            tau_list = []
            p_values = []
            plt_distribution = []
            plt_cdfs = []

            pfile = open(F"{output_data}/pvalues{rid}.dat", 'w')
            pfile.write("#force\t p-value\n")
            force_list = force_list*(1.66)*-1
            for force in escape_times:
                time_list = escape_times[force]
                escape_file = open(F"{output_data}/escape_times_{force}{rid}.dat", 'w')
                escape_file.write("#run\t escape_time\n")
                print(force)
                for i in range(len(time_list)):
                    escape_file.write("%d\t %f\n"%(i+1, time_list[i]))
                escape_file.close()
                plt_distribution.append(metaD.plot_kde_hist_cumulative(time_list, rid, out_dir=output,force=force))
                cdfs, data = metaD.compute_cdfs(time_list)
                time_domain = cdfs[0]
                ecdf = cdfs[1]
                tcdf = cdfs[2]
                tau = data[0]
                tau_list.append(tau)
                plt_cdfs.append(metaD.plot_cdfs(x=time_domain, y=ecdf, tcdf=tcdf, id=rid, force=force, model_name=model_name, out_dir=output))
                s, p = metaD.do_KS2test(time_list, tau)
                p_values.append(p)
                pfile.write("%f\t%f\n"%(float(force),p))
            pfile.close()

            rates = 1/np.array(tau_list)

            error_bar = metaD.calculate_error_bars(tau_list, n_repeats)

            taus_b = error_bar[:,0]
            bars = error_bar[:,1]
            rates_b = error_bar[:,2]
            r_bars = error_bar[:,3]

            taus_bars = open(F"{output_data}/taus_bars{rid}.dat", 'w')
            rates_bars = open(F"{output_data}/rates_bars{rid}.dat", 'w')
            rates_bars.write("#%s\t%s\t%s\t%s\n"%("forces","rates","rates_bootstrap", "bars"))
            taus_bars.write("#%s\t%s\t%s\t%s\n"%("forces","taus","taus_bootstrap", "bars"))

            for i in range(len(force_list)):
                taus_bars.write("%f\t%f\t%f\t%f\n"%(force_list[i],tau_list[i],taus_b[i],bars[i]))
                rates_bars.write("%f\t%f\t%f\t%f\n"%(force_list[i],rates[i],rates_b[i],r_bars[i]))
            rates_bars.close()
            taus_bars.close()

            rates_plot_path = metaD.plot_catch_tau_rates(force_list, tau_list, model_name, output, rid, error_bar)

            cb_params = metaD.catch_params(force_list,tau_list)
            k1c_0 = round(cb_params[0],3)
            x1c = round(cb_params[1],3)
            k1s_0 = round(cb_params[2],3)
            x1s = round(cb_params[3],3)

            rows_rate = ['in:nruns,img:rates,out:k1c_0,out:x1c,out:k1s_0,out:x1s']
            plt_g = os.path.join(run_dir, rates_plot_path)
            rows_rate.append(F'{n_repeats},{plt_g},{k1c_0},{x1c},{k1s_0},{x1s}')

            with open(F'{output_data}/{model_name}_rates.csv','w') as csvfile:
                csvfile.write("\n".join(rows_rate))
            with open(F'{output_data}/{model_name}_rates.html', 'w') as html:
                csvpath=os.path.join(run_dir,output_data,F'{model_name}_rates.csv')
                html.write("""\
                <html style="overflow-y:hidden;background:white">\
                <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px"\
                href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>\
                <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>\
                </html>""".format(csvpath,'nruns', csvpath,'nruns'))


            rows = ['in:force, img:kde, img:traj, img:fit, img:dist, out:tau, out:pval']
            for i,force in enumerate(escape_times):
                plt_a = os.path.join(run_dir,plt_kde_force[i])
                plt_b = os.path.join(run_dir,plt_trj_time[i])
                plt_c = os.path.join(run_dir,plt_cdfs[i])
                plt_d = os.path.join(run_dir,plt_distribution[i])
                plt_e = tau_list[i]
                plt_f = p_values[i]

                rows.append(F"{-1*float(force)},{plt_a},{plt_b},{plt_c},{plt_d},{plt_e},{plt_f}")

            with open(F'{output_data}/{model_name}.csv','w') as csvfile:
                csvfile.write("\n".join(rows))
            with open(F'{output_data}/{model_name}.html', 'w') as html:
                csvpath=os.path.join(run_dir,output_data,F'{model_name}.csv')
                html.write("""\
                <html style="overflow-y:hidden;background:white">\
                <a style="font-family:sans-serif;z-index:1000;position:absolute;top:15px;right:0px;margin-right:20px;font-style:italic;font-size:10px"\
                href="/preview/DesignExplorer/index.html?datafile={}&colorby={}" target="_blank">Open in New Window</a>\
                <iframe width="100%" height="100%" src="/preview/DesignExplorer/index.html?datafile={}&colorby={}" frameborder="0"></iframe>\
                </html>""".format(csvpath,'force', csvpath, 'tau'))
