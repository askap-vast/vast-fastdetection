import os
import numpy as np


class DataBasic:
    '''
    Create VASTER folder saving system
    '''

    def __init__(self, sbid, parent_dir):
        self.sbid = sbid
        self.parent_dir = parent_dir 
        self.parent_dir_abs = os.path.abspath(parent_dir)

        self.nbeam = 36
        self.steps = [
            'GETDATA', 'UNTAR', 'FIXDATA', 'MODELING', 'IMGFAST', 'SELCAND', 'CLNDATA', 
        ]
        self.function = [
            'run','output', 'error', 'usage', 
        ]


    def create_folder_tree(self, path):
        paths = {}
        paths['path'] = path
        paths['path_data'] = os.path.join(path, 'data') # saving visibilities, selavy catalogues 
        paths['path_models'] = os.path.join(path, 'models') 
        paths['path_images'] = os.path.join(path, 'images')
        paths['path_cand'] = os.path.join(path, 'candidates')
        paths['path_scripts'] = os.path.join(path, 'scripts')
        paths['path_logs'] = os.path.join(path, 'logfiles')

        return paths
    
    @property
    def paths(self, ):
        path = os.path.join(self.parent_dir_abs, f'SB{self.sbid}')
        return self.create_folder_tree(path)
    

    @property
    def files(self, ):

        dtype = [
            ('beam', np.int32), 
            ('step', 'U10'), 
            ('function', 'U10'), 
            ('fname', 'U256')
        ]

        shape = (self.nbeam, len(self.steps), len(self.function), 1)

        array = np.empty(shape, dtype=dtype)

        for idx in range(self.nbeam):
            for j, step in enumerate(self.steps):
                if step == 'GETDATA' or step == 'UNTAR':
                    prefix = 'bash'
                else:
                    prefix = 'slurm'

                for k, function in enumerate(self.function):
                    if function == 'run':
                        name = f'{step}_beam{idx:02d}'
                        fname = os.path.join(self.paths['path_scripts'], f'{prefix}_{name}.sh')
                    else:
                        name = f'{step}_SB{self.sbid}_beam{idx:02d}'
                        fname = os.path.join(self.paths['path_logs'], f'{prefix}_{name}.{function}')

                    if prefix == 'bash' and (function == 'error' or function == 'usage'):
                        fname = None

                    array[idx, j, k, 0] = (idx, step, function, fname)

        return array




