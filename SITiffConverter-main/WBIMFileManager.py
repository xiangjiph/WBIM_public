from lzma import CHECK_NONE
import os.path
import platform
import h5py
import tifffile
import time
import numpy as np

class WBIMFileManager:

    def __init__(self, hostname=None):
        if hostname is None:
            self.HOSTNAME = platform.node()
        else:
            self.HOSTNAME = hostname.lower()

        if self.HOSTNAME == 'pia':
            self.SCRIPT_PATH = 'C:\\Users\\dklab\\Documents\\GitHub\\WBIM'
            self.SCRATCH_ROOT_PATH = 'H:\\WBIM'
            self.DATA_ROOT_PATH = 'D:\\WBIM'
            self.SERVER_ROOT_PATH = 'Z:\\WBIM'
            self.PYTHON_SCRIPT_PATH = 'C:\\Users\\dklab\\Documents\\GitHub\\SITiffConverter'
        elif self.HOSTNAME == 'DKLab-Precision7530':
            self.SCRIPT_PATH = 'C:\\Users\\xij072\\Documents\\GitHub\\WBIM'
            self.SCRATCH_ROOT_PATH = 'C:\\data\\WBIM'
            self.DATA_ROOT_PATH = self.SCRATCH_ROOT_PATH
            self.PYTHON_SCRIPT_PATH = 'c:\\Users\\xij072\\Documents\\GitHub\\SITiffConverter'
        elif self.HOSTNAME == 'XiangJi-PC':
            self.DATA_ROOT_PATH = 'D:\\data\\Vessel\\WBIM'
            self.SCRATCH_ROOT_PATH = self.DATA_ROOT_PATH
            self.SCRIPT_PATH = 'D:\\Github\\WBIM'
        elif self.HOSTNAME in ['bird', 'bird.dk.ucsd.edu']:
            self.DATA_ROOT_PATH = '/nfs/birdstore/Vessel/WBIM'
            self.SCRIPT_PATH = '/home/xij072/Documents/Github/WBIM'
            self.SCRATCH_ROOT_PATH = '/scratch/Vessel/WBIM'
            self.SERVER_ROOT_PATH = '/nfs/birdstore/WBIM'

    @property
    def fp_acquisition_scratch(self):
        return os.path.join(self.SCRATCH_ROOT_PATH, 'Acquisition')

    @property
    def fp_acquisition_disk(self):
        return os.path.join(self.DATA_ROOT_PATH, 'Acquisition')
    
    @property
    def fp_acquisition_server(self):
        return os.path.join(self.SERVER_ROOT_PATH, 'Acquisition')


def load_h5(h5_fp):
    # Internal function
    def _load_h5_node(f):
        result = {}
        for key in f.keys():
            if isinstance(f[key], h5py.Dataset):
                result[key] = {}
                result[key]['data'] = f[key][:]
                # Deal with attributes

            else:
                result[key] = _load_h5_node(f[key])
        return result

    with h5py.File(h5_fp, 'r') as f:
        return _load_h5_node(f)

def get_si_metadata(fp):
    with tifffile.TiffFile(fp) as tif:
        metadata = tif.pages[0].tags
        tag_dic = {}
        for tag in metadata.values():
            tag_dic[tag.name] = tag.value
        info = tag_dic['Software'] + '\n\n' + tag_dic['Artist']
        dtype = tif.pages[0].dtype
        return (info, tif.pages[0].dtype)
    
def load_si_tiff_using_tifffile(fp):
    # This line can be replaced by a lower-performance one. 
    im_vol = tifffile.imread(fp)
    # Reshape to match scanimage output format
    if im_vol.ndim == 4 and im_vol.shape[1] <= 4: 
        im_vol = im_vol.reshape(-1, *im_vol.shape[2:])
    info, dtype = get_si_metadata(fp)
    return im_vol, info, dtype

def load_si_tiff_slow(fp, pause_step=50, pause_time_s=0.1):

    with tifffile.TiffFile(fp) as tif:
        im_vol = []
        for i, page in enumerate(tif.pages):
            im_vol.append(page.asarray())
            if i % pause_step == 0:
                time.sleep(pause_time_s)
        im_vol = np.stack((im_vol))

    if im_vol.ndim == 4 and im_vol.shape[1] <= 4: 
        im_vol = im_vol.reshape(-1, *im_vol.shape[2:])
    info, dtype = get_si_metadata(fp)
    return im_vol, info, dtype



def check_folder_existence(fp, createQ=True):
    parent_dir, dir_name = os.path.split(fp)
    parent_dir_exist_Q = os.path.isdir(parent_dir)
    if createQ and (not parent_dir_exist_Q):
        os.makedirs(parent_dir, exist_ok=True)
    return parent_dir_exist_Q