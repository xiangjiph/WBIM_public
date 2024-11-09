from calendar import day_name
from keyword import iskeyword
from tkinter import image_names
import ScanImageTiffReader as SITR
import numpy as np
import os
import logging
import traceback
import os.path
import json
import platform
import h5py
import tifffile
import shutil
import argparse
from time import sleep
import WBIMFileManager
import TileAnalysis as TA
from RowShiftDetector import RowShiftDetector


class SITiffParser:
    __PIA_ACQ_DISK = "D:\\WBIM\\Acquisition"
    __PIA_SERVER_DISK = "Z:\\WBIM\\Acquisition"
    __num_read_attempt = 2
    __read_attempt_delay_time_s = 3
    __try_save_h5_to_server_Q = True

    def __init__(self, info_json_fp):
        self.DM = WBIMFileManager.WBIMFileManager()
        # Identify system
        self.hostname = platform.node().lower()
        # Processing parameters
        self.SI_CONVERT_IM_TYPE = np.uint16
        self.downsample_pixel_size_zxy_um = np.array([5, 5, 5])
        # Disable lineshift correction at the moment
        self.apply_lineshift_correction_Q = False
        self.est_lineshift_num_col = 15
        self.doneQ = False
        self.rm_last_frame_in_scan_mode_Q = True
        self.row_shift_detectedQ = False
        # Row shift detection 
        self.rs_min_fg_val = 3e3
        self.rs_sm_wd = 9
        # Inter-percentile for determine outlier for each channel
        self.rs_ipr = [3.0, 3.0, 3.0, 3.0]  

        # Read json file
        self.info_file_path = info_json_fp
        self.info_folder_path = os.path.dirname(info_json_fp)
        # Load metadata and set logger
        with open(info_json_fp) as f:
            self.info = json.load(f)

        self.parse_tile_info()
        self.configure_logger()
        try:
            # Save to H5 if the data is loaded from SI Tiff
            self.load_raw_data()
            self.compute_mip_and_global_lineshift()
            # 3 for 1 um resolution, 5 for 0.5 um resolution, assume x and y have
            # the same pixel size
            # medfilt2_size = int(2 * (1 / self.info['pixel_size_um'][1]) + 1)
            # self.remove_shot_noise_2d_in_place(medfilt2_size)
            # Recompute MIP after median filtering?

            # Use the same function to process both the exploration and scanning data
            self.process_scan_data()
            # if self.info['acq_mode'] == 'Scan':
            #     self.process_scan_data()
            # elif self.info['acq_mode'] == 'Explore':
            #     self.process_exploration_data()
            TA.save_dict_as_mat_file(self.info['stat_fp'], self.stat)
            # Need to verify all the processings have been completed properly
            self.delete_si_tiff()
            self.logger.info("Successfully finish processing SI Tiff file")
            self.doneQ = True
        except:
            self.logger.exception('Got exception in the constructor')

    def parse_tile_info(self):
        # For local testing
        if self.hostname != 'pia':
            self.info['save_root_folder'] = self.info['save_root_folder'].replace(
                self.__PIA_ACQ_DISK, self.DM.fp_acquisition_disk)
        else:
            # On PIA, try to save to the network drive first?  
            pass

        if type(self.info['channel']) == int:
            self.info['channel'] = [self.info['channel']]
        self.info['pixel_size_um'] = np.array(self.info['pixel_size_um'])
        self.info['pixel_size_zxy_um'] = self.info['pixel_size_um'][[2, 1, 0]]
        self.downsample_pixel_size_zxy_um = np.maximum(self.downsample_pixel_size_zxy_um,
                                                       self.info['pixel_size_zxy_um'])

        self.info['h5_fp'] = os.path.join(self.info['save_root_folder'],
                                          self.info['fprr_tile'])
        self.info['h5_group'] = []
        for group_name in self.info['h5_dataset']:
            # The '/' is added for H5 IO in MATLAB
            if group_name.startswith('/'):
                group_name = group_name[1::]
                self.info['h5_group'].append(group_name)

        # Initialize stat
        if 'fprr_stat' not in self.info.keys():
            self.info['fprr_stat'] = self.info['fprr_info'].replace(
                'info.mat', 'stat.mat')
        self.info['stat_fp'] = os.path.join(self.info['save_root_folder'],
                                            self.info['fprr_stat'])
        self.stat = {}
        for ch_id in self.info['channel']:
            self.stat['CH{0}'.format(ch_id)] = {}

    def configure_logger(self):
        if self.hostname.lower() == 'pia':
            log_filename = os.path.join(self.info['save_root_folder'],
                                        self.info['fprr_pplog'])
        else:
            log_filename = os.path.join(self.info_folder_path,
                                        'log.txt')
        log_folder = os.path.dirname(log_filename)
        if os.path.isdir(log_folder):
            Warning("Log folder does not exist. Create folder")
            os.makedirs(log_folder, exist_ok=True)
        # err_filename = os.path.join(self.info['save_root_folder'], \
        #     self.info['fprr_pperr'])

        log_format = '%(asctime)s %(message)s'
        # logging.basicConfig(\
        #     level=logging.DEBUG, format=log_format, \
        #     handlers=[logging.FileHandler(log_filename),
        #             logging.StreamHandler()])
        logging.basicConfig(
            level=logging.DEBUG, format=log_format,
            handlers=[logging.FileHandler(log_filename)])

        self.logger = logging.getLogger("WBIMSIParser")
        self.logger.debug("Finish initializing logger")

    def load_raw_data(self, si_fp=None):
        if si_fp is None:
            si_fp = self.info['SI_filepath']

        si_fp_2 = os.path.join(self.info['save_root_folder'],
                               self.info['fprr_raw'])
        if os.path.isfile(si_fp):
            self.convert_SI_tiff(si_fp)
        elif os.path.isfile(si_fp_2):
            self.convert_SI_tiff(si_fp_2)
        else:
            h5_data = WBIMFileManager.load_h5(self.info['h5_fp'])
            self.im_stack = []
            for chid in self.info['channel']:
                tmp_ch = h5_data['CH{}'.format(chid)]
                if 'raw' in tmp_ch.keys():
                    self.im_stack.append(
                        tmp_ch['raw']['data'])
                else:
                    self.im_stack.append(tmp_ch['data'])

    def convert_SI_tiff(self, si_fp=None):
        # The loaded image stack is in the correct axis for visualization, for
        # writing tiff mip. However, for images, the y and x axis need to be
        # transposed after loading into MATLAB. When ImageJ load the image from
        # the h5 file, the image is also transposed.

        if si_fp is None:
            si_fp = self.info['SI_filepath']
            if not os.path.isfile(si_fp):
                si_fp = os.path.join(self.info['save_root_folder'],
                                     self.info['fprr_raw'])
                assert os.path.isfile(
                    si_fp), "SI raw tiff {} does not exist".format(si_fp)
        si_metadata_fp = os.path.join(self.info['save_root_folder'],
                                      self.info['fprr_si_info'])

        attempt_count = 0
        while attempt_count < self.__num_read_attempt:
            try:
                # Fast version
                # [raw_data, self.si_metadata,
                #    self.im_dtype] = self.load_si_tiff(si_fp)
                # Slow version - to reduce competition with SI writer 
                [raw_data, self.si_metadata,
                    self.im_dtype] = WBIMFileManager.load_si_tiff_slow(si_fp)
                # The raw data from SI reader is in [z, y, x] order -> no need
                # to transpose the first two dimension
                self.logger.info("Finish parsing SI tiff stack")
                self.write_txt_to_file(si_metadata_fp, self.si_metadata)
                self.logger.info("Finish writing SI metadata")
                self.im_stack = []
                # Transpose the array - not necessary after changing the axis?
                # self.im_stack = np.transpose(self.im_stack, axes=[0,2,1])
                num_ch = np.size(self.info['channel'])
                if raw_data.dtype is not self.SI_CONVERT_IM_TYPE:
                    raw_data = TA.im_int16_to_uint16(raw_data)
                    self.logger.info("Finish converting to uint16")

                if num_ch == 1:
                    self.im_stack = [raw_data]
                else:
                    num_channel = np.size(self.info['channel'])
                    assert np.mod(
                        raw_data.shape[0], num_channel) == 0, "Mismatch image stack shape"
                    for i in range(num_ch):
                        self.im_stack.append(raw_data[i::num_ch, :, :])
                self.save_SI_data()
                break
            except:
                attempt_count += 1
                if attempt_count < self.__num_read_attempt:
                    sleep_time = self.__read_attempt_delay_time_s
                    self.logger.warning(
                        "Cannot open the existing SI TIFF file for {count:d} times.".format(count=attempt_count))
                    sleep(sleep_time)
                else:
                    raise

    def save_SI_data(self, filepath=None):
        if filepath is None:
            filepath = self.info['h5_fp']

        if self.__try_save_h5_to_server_Q and os.path.isdir(self.__PIA_SERVER_DISK):
            new_filepath = filepath.replace(self.__PIA_ACQ_DISK, self.__PIA_SERVER_DISK)
            self.logger.info("Server mounted. Try saving h5 data to the server")
            try:
                self._write_raw_h5_data(new_filepath)
            except Exception as e:
                self.logger.error("Error message: %s", e)
                self.logger.error(traceback.format_exc())
                self.logger.info("Fail to write h5 data to the server. Try to save to local drive.")
                self._write_raw_h5_data(filepath)
        else:
            self._write_raw_h5_data(filepath)

    def _write_raw_h5_data(self, filepath):
        WBIMFileManager.check_folder_existence(filepath, createQ=True)
        with h5py.File(filepath, 'a') as f:
            for i, group_name in enumerate(self.info['h5_group']):
                if group_name not in f.keys():
                    h5_group = f.create_group(group_name)
                elif isinstance(f[group_name], h5py.Dataset):
                    del f[group_name]
                    h5_group = f.create_group(group_name)
                else:
                    h5_group = f[group_name]

                if 'raw' in h5_group.keys():
                    self.logger.warning(
                        "Target dataset already exist. Delete current data")
                    del h5_group['raw']
                dset = h5_group.create_dataset(
                    'raw', data=self.im_stack[i], chunks=True)
        self.logger.info("Finish writing SI raw data")

    def rm_last_section(self):
        # To be tested on PIA
        # For removing the piezo flyback artefact
        if self.rm_last_frame_in_scan_mode_Q and self.info['acq_mode'] == 'Scan':
            for i, im in enumerate(self.im_stack):
                self.im_stack[i] = im[:-1]

    def compute_mip_and_global_lineshift(self):
        self.mip = []
        self.lineshift = []
        for i, im in enumerate(self.im_stack):
            if im.ndim == 3:
                ch_mip = np.amax(im, axis=0)
            else:
                ch_mip = im
            # Estimate global lineshift
            try:
                assert ch_mip.shape[0] > 2 * self.est_lineshift_num_col,\
                    f"The image has less than {2 * self.est_lineshift_num_col} rows!"
                ch_shift = TA.compute_global_column_shift(
                    ch_mip[:, self.est_lineshift_num_col:-self.est_lineshift_num_col])
                self.logger.info(
                    "Global lineshift is {0} in {1}".format(ch_shift, self.info['h5_group'][i]))
                if ch_shift != 0 and self.apply_lineshift_correction_Q:
                    # im and self.im_stack[i] point to the same value
                    # im = TA.correct_global_lineshift(
                    #     im, ch_shift)
                    ch_mip = TA.correct_global_lineshift(ch_mip, ch_shift)
            except:
                ch_shift = np.nan
                self.logger.exception("Failed to estimate the line shift")
            # Save MIP tiff
            tmp_mip_fp = os.path.join(self.info['save_root_folder'],
                                      self.info['fprr_mip'][i])
            tifffile.imwrite(tmp_mip_fp, ch_mip)
            self.mip.append(ch_mip)
            self.stat['CH{0}'.format(
                self.info['channel'][i])]['line_shift'] = ch_shift

    def process_exploration_data(self):
        for i, im in enumerate(self.im_stack):
            ch_name = 'CH{0}'.format(self.info['channel'][i])
            self.stat[ch_name]['z_stat_raw'] = TA.comptue_z_statistics(im)
            medfilt2_size = np.uint8(
                2 * (1 / self.info['pixel_size_um'][1]) + 1)
            im_sm = TA.remove_shot_noise_2d_in_plane(im, medfilt2_size)
            self.stat[ch_name]['z_stat_sm'] = TA.comptue_z_statistics(im_sm)
            stat_bin_edge = np.linspace(0, np.iinfo(self.mip[i].dtype).max, 10)
            self.stat[ch_name]['mip_stat_raw'] = TA.compute_basic_statistics(
                self.mip[i], stat_bin_edge)
            self.stat[ch_name]['mip_stat_sm'] = TA.compute_basic_statistics(
                np.amax(im_sm, axis=0), stat_bin_edge)

    def process_scan_data(self):
        for i, im in enumerate(self.im_stack):
            ch_name = 'CH{0}'.format(self.info['channel'][i])
            medfilt2_size = np.uint8(np.minimum(
                5, 2 * (1 / self.info['pixel_size_um'][1]) + 1))
            im = TA.remove_shot_noise_2d_in_plane(im, medfilt2_size)
            step_mip_step = np.round(
                self.downsample_pixel_size_zxy_um[0] / self.info['pixel_size_zxy_um'][0])
            self.stat[ch_name]['step_mip'] = TA.compute_step_mip(
                im, step_mip_step)
            # Resize 2d image
            target_size_2d = np.round(im.shape[1:3] * self.info['pixel_size_zxy_um'][1:3] /
                                      self.downsample_pixel_size_zxy_um[1:3])

            self.stat[ch_name]['step_mip'] = TA.resize_image_stack_2d(self.stat[ch_name]['step_mip'],
                                                                      target_size_2d.astype(np.int32))
            self.stat[ch_name]['step_mip_pixel_yxz_um'] = self.downsample_pixel_size_zxy_um[[
                2, 1, 0]]
            if self.info['acq_mode'] == 'Scan':
                self.stat[ch_name]['mip'] = TA.compute_three_view_mip(im)
            self.stat[ch_name]['z_stat'] = TA.comptue_z_statistics(im)

            # Detect row shift
            self.stat[ch_name]['flyback_row_ind'] = np.repeat(np.nan, (2, ))
            try:
                mip_im = self.mip[i][:, self.est_lineshift_num_col:-self.est_lineshift_num_col]
                rsd = RowShiftDetector(mip_im, sm_wd=self.rs_sm_wd, min_fg_val=self.rs_min_fg_val, 
                                       ipr=self.rs_ipr[i])
                if rsd.max_oi_ep_ind is not None:
                    self.stat[ch_name]['flyback_row_ind'] = rsd.max_oi_ep_ind
            except:
                self.logger.exception("Fail to estimate the row shift artefact")
  
    @property
    def num_row_shift(self):
        # Estimation from the vascular channel is more reliable 
        num_row = self.info['stack_size'][0]
        for i, chid in enumerate(self.info['channel']):
            ch_name = 'CH{0}'.format(chid)
            ep_ind = self.stat[ch_name]['flyback_row_ind']
            ep_med = (ep_ind[0] + ep_ind[1]) / 2
            if ~np.any(np.isnan(ep_ind)):
                if ep_med < (num_row - ep_med):
                    # Closer to the top of the image
                    return ep_ind[1]
                else: 
                    return ep_ind[0] - num_row
        return 0
    
    # About files
    def move_si_tiff(self, source_fp=None, dest_fp=None):
        if source_fp is None:
            source_fp = self.info['SI_filepath']
        if dest_fp is None:
            dest_fp = os.path.join(self.info['save_root_folder'],
                                   self.info['fprr_raw'])

        if os.path.isfile(source_fp) and not os.path.isfile(dest_fp):
            attempt_count = 0
            while attempt_count < self.__num_read_attempt:
                try:
                    shutil.move(source_fp, dest_fp)
                    self.logger.info("Finish moving scanimage raw tiff file")
                    break
                except:
                    attempt_count += 1
                    if attempt_count <= self.__num_read_attempt:
                        self.logger.warning(
                            "Cannot move file for {count:d} times.".format(count=attempt_count))
                        sleep(self.__read_attempt_delay_time_s)
                    else:
                        raise

    def delete_si_tiff(self):
        if os.path.isfile(self.info['SI_filepath']):
            os.remove(self.info['SI_filepath'])
            self.logger.info("Finish deleting the SI Tiff file")

    # Static methods
    @staticmethod
    def write_txt_to_file(filepath, txt):
        f = open(filepath, "w")
        f.write(txt)
        f.close()

    @staticmethod
    def load_si_tiff(si_file_path):
        if not os.path.isfile(si_file_path):
            print("File does not exist")
            return
        else:
            with SITR.ScanImageTiffReader(si_file_path) as reader:
                vol = reader.data()
                # vol is in order [z, y, x]
                # Need to check whether it is consistent with the one viewed in imageJ
                info = reader.metadata()
                d_type = reader.dtype()
                reader.close()
                return [vol, info, d_type]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert scanimage raw tiff stack to HDF5")
    parser.add_argument('-f', '--filepath', metavar='info_fp', required=True,
                        help='Absolute file path to the json file of the tile metadata')
    args = parser.parse_args()
    converter = SITiffParser(args.filepath)
    # Testing on Dell
    # test_folder = "C:\\data\\WBIM\\Acquisition\\TestSample\\WBIM20220616\\00001\\Scan\\00049_00072\\00000"
    # info_json_fp = os.path.join(test_folder, "info.json")
    # wbim_parser = SITiffParser(info_json_fp)

    # Testing on PIA
