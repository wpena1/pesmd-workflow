
import subprocess
import parsl
import os,sys,logging

from parsl.config import Config

from parsl.providers import CoasterProvider
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.monitoring.monitoring import MonitoringHub

from parsl.data_provider import pwrsync
from parsl.data_provider.pwrsync import PWRsyncStaging

'''
   Expected in the current working dir when this is imported:
     parsl.swift.conf
     parsl.pool.properties

  For testing:
    pool.properties_pwcoaster -> parsl-start-coaster-service.sh -> pool.properties
      (inserts session #)
'''

def bash(script,*args,checkrc=False,debug=False):

    cmd = ['bash', '-s', *args]

    result = subprocess.run(cmd, input=script, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8')
    if debug: print('DB: bash(): subprocess.run returned: ' + str(result))

    rc = result.returncode
    if debug:  print('DB: bash(): command exit code: ' + str(rc))

    cmd_output = result.stdout
    if debug: print('DB: bash(): command output: ' + str(cmd_output))

    return (rc, cmd_output)

try:
    swift_conf      = 'pw.conf'
    pool_properties = 'pool.properties'  # Link this to the coaster pool used

    try:
        script = r'''
        set -o pipefail
        sed -e 's/{{//g;s/}}//g;s/"//g;s/: / /;s/^ +//' < {0} |
            awk '
                /^site\./   {{ gsub(/^site./,"",$1); site=$1 }}
                /  URL/     {{ print(site " " $2) }}
                '
        '''.format(swift_conf)   # FIXME: Convert to python???
    except:
        swift_conf = 'swift.conf'
        script = r'''
        set -o pipefail
        sed -e 's/{{//g;s/}}//g;s/"//g;s/: / /;s/^ +//' < {0} |
            awk '
                /^site\./   {{ gsub(/^site./,"",$1); site=$1 }}
                /  URL/     {{ print(site " " $2) }}
                '
        '''.format(swift_conf)   # FIXME: Convert to python???

    (rc,output) = bash(script)

    if rc != 0:
        #print('parslpw: error: can not parse parsl.swift.conf. rc={} output={}'.format(rc,output))
        # try:
        #     exit(1)
        # except:
        raise ValueError('parslpw: error: can not parse pw.conf. rc={} output={}'.format(rc,output))

    pools=[]
    poolinfo={}
    for line in output.splitlines():
        fields = line.split()
        poolname = fields[0]
        service = fields[1].split(':')
        host = service[1][2:] # skip over // from protocol
        port = service[2]
        pools.append(poolname)
        poolinfo[poolname]={'host':host,'port':port}

    #print(pools)
    #print(poolinfo)

    pool=pools[0]
    coaster_host = poolinfo[pool]['host']
    coaster_port = poolinfo[pool]['port']

    #print('pool: {} coaster_host: {} coaster_port: {}'.format(pool,coaster_host,coaster_port))

    #print(os.path.basename(__file__) + " entered")

    def get_coaster_properties(pfile):
        props = {}
        with open(pfile) as f:
            for line in f:
                line = line.rstrip().split('#',1)[0]
                toks = line.split(':',maxsplit=1)
                if len(toks) == 2:
                    name = toks[0].lstrip().rstrip()
                    val = toks[1].lstrip().rstrip()
                    if len(val) >= 2 and ((val[0] == "'" and val[-1] == "'") or (val[0] == '"' and val[-1] == '"')):
                        val = val[1:-1]
                    props[name] = val
        return props

    #print("Running parsl version " + parsl.__version__)

    #print("test of PW Coasters starting")

    # get the props from the coaster service

    pwprops = get_coaster_properties(swift_conf)
    #print(pwprops)

    os.environ['PW_PARSL_WORKER_HOME'] = pwprops['workDirectory']
    os.environ['keepSiteDir'] = pwprops['keepSiteDir']

    pwconfig = Config(
        executors=[
            HighThroughputExecutor(
                label='coaster_single_node',
                worker_debug=False,
                cores_per_worker=int(pwprops['CORES_PER_WORKER']),
                working_dir=pwprops['workDirectory'],
                worker_logdir_root=pwprops['workDirectory'] + '/parsllogs',
                # storage_access=[PWRsyncStaging(hostname='localhost')], # FIXME: Remove hostname?
                storage_access=[PWRsyncStaging(hostname=None)], # FIXME: Remove hostname?
                address=coaster_host.replace('http://','').split(':')[0],   ### # FIXME?
                provider = CoasterProvider(
                    coaster_host = coaster_host,
                    coaster_port = coaster_port,
                    init_blocks = 1,
                    min_blocks = int(pwprops['BLOCK_MIN']),
                    max_blocks = int(pwprops['BLOCK_MAX']),
                    parallelism = 1 # Was 0.80 ### # FIXME?
                )
            )
        ],
        monitoring=MonitoringHub(
           hub_address=address_by_hostname(),
           resource_monitoring_interval=1,
       ),
    )
    #print("Config completed. pwconfig:", pwconfig)

except Exception as e:
    #print(e)

    print("loading local pwconfig")
    pwconfig = Config(
        executors=[
            HighThroughputExecutor(
                storage_access=[PWRsyncStaging(hostname=None)],
                label='coaster_single_node', # FIXME: Allow multinode
                worker_debug=False,
                cores_per_worker=1,
                working_dir='./',
                provider = LocalProvider(
                    init_blocks = 1,
                    min_blocks = 1,
                    max_blocks = 1,
                    parallelism = 1 # Was 0.80 ### # FIXME?
                )
            )
        ],
        monitoring=MonitoringHub(
           hub_address=address_by_hostname(),
           resource_monitoring_interval=1,
       ),
    )

# GET COMMAND LINE ARGS FROM PW FORM
# import argparse
# parser=argparse.ArgumentParser()
# parsed, unknown = parser.parse_known_args()
# for arg in unknown:
#     if arg.startswith(("-", "--")):
#         parser.add_argument(arg)
# pwargs=parser.parse_args()
# print("pwargs:",pwargs)
