import os

class DataBasic:

    def __init__(self, ):
        self.nbeam = 36
        self.steps = [
            'GETDATA', 'UNTAR', 'FIXDATA', 'MODELING', 'IMGFAST', 'SELCAND', 'CLNDATA', 
        ]



class DataDir:
    '''
    Create VASTER folder saving system
    '''

    def __init__(self, sbid, parent_dir):
        self.sbid = sbid
        self.parent_dir = parent_dir 
        self.parent_dir_abs = os.path.abspath(parent_dir)


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
    

    
    # @property
    # def files(self, ):
    #     f = {}
    #     f_beam = {}
        
    #     for step in self.steps:
    #         for idx in range(self.nbeam):
    #             name = f'{step}_SB{self.sbid}_beam{idx:02d}'
    #             f_beam['log_usage'] = os.path.join(self.paths['path_logs'], f'slurm_{name}.usage')
    #             f_beam['log_output'] = os.path.join(self.paths['path_logs'], f'slurm_{name}.output')
    #             f_beam['log_error'] = os.path.join(self.paths['path_logs'], f'slurm_{name}.error')

    #             f[idx] = f_beam

    #     return f




