'''Functions to set up a kubernetes cluster for downloading data.
'''

from dask_kubernetes.operator import KubeCluster, make_cluster_spec

def make_k8s_download_cluster(name, out_dir, n_workers=3, threads_per_worker=6,
                              worker_memory='2G', local_domain=None,
                              port_forward_cluster_ip=None):
    '''Set up a kubernetes cluster for downloading NWP data. Each worker has 1
    CPU, and will write data to `out_dir`. Requires the dask kubernetes operator
    to be installed on the cluster.
    '''
    worker_resources = {'limits': {'cpu': '1', 'memory': worker_memory}}
    if local_domain is None:
        scheduler_address = f'tcp://{name}-scheduler:8786'
    else:
        scheduler_address = f'tcp://{name}-scheduler.svc.{local_domain}:8786'
    env = {"DASK_SCHEDULER_ADDRESS": scheduler_address,
           'EXTRA_APT_PACKAGES': 'curl libeccodes0',
           'EXTRA_CONDA_PACKAGES': 'wgrib2',
           'EXTRA_PIP_PACKAGES': 'herbie-data humanize'}
    spec = make_cluster_spec(name=name, n_workers=n_workers,
                             resources=worker_resources, env=env)
    # give the scheduler enough time to install the conda packages
    scheduler = spec['spec']['scheduler']['spec']['containers'][0]
    scheduler['livenessProbe']['initialDelaySeconds'] = 300
    # add the volume for writing data
    volume = {'hostPath': {'path': out_dir},
              'name': 'nwpout'}
    mount = {'name': 'nwpout', 'mountPath': '/mnt/nwp'}
    spec['spec']['worker']['spec']['volumes'] = [volume]
    worker = spec['spec']['worker']['spec']['containers'][0]
    worker['volumeMounts'] = [mount]
    worker['args'].extend(['--nthreads', str(threads_per_worker)])
    return KubeCluster(custom_cluster_spec=spec,
                       port_forward_cluster_ip=port_forward_cluster_ip)
