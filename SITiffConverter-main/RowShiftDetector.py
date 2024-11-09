import numpy as np
import scipy as sp
import cv2
from scipy import ndimage 
import TileAnalysis as TA

# To do: 
# 1. Replace percentile outlier estimation with moving percentile 
# 2. 

class RowShiftDetector():
    def __init__(self, im, min_fg_val=5e3, sm_wd=9, max_width_s=0.003, ipr=3.0):
        im = im.astype(np.single)
        self.im = im
        self.min_fg_val = min_fg_val
        self.sm_wd = sm_wd
        self.ipr = ipr
        # 2 ms for the maximum flyback artefact
        self.max_artefact_width = max_width_s *  7920 * 2

        self.compute()

    def compute(self):
        self.suppress_background()
        self.row_stat = TA.compute_row_stat(self.im)
        self.detect_stat_outlier()
        self.select_outlier_intervals()
    
    def binarize_image(self, th=1.5e4):
        self.im = self.im > th

    def suppress_background(self):
        self.im[self.im <self.min_fg_val] = 0

    def detect_stat_outlier(self, min_corr_th=0.01, min_ig_th=0.0, min_ig2_th=0.0):
        # Merge small nearby connected components

        self.corr_th = max(min_corr_th, TA.compute_percentile_outlier_threshold(self.row_stat['adj_corr_da'], ipr=self.ipr)[1])
        corr_outlier_Q = TA.binary_close_1d(self.row_stat['adj_corr_da'] > self.corr_th, self.sm_wd)

        self.int_grad_th = max(min_ig_th, TA.compute_percentile_outlier_threshold(self.row_stat['row_abs_diff'], ipr=self.ipr)[1])
        int_grad_outlier_Q = TA.binary_close_1d(self.row_stat['row_abs_diff'] > self.int_grad_th, self.sm_wd)

        # self.int_grad2_th = max(min_ig2_th, TA.compute_percentile_outlier_threshold(self.row_stat['row_abs_diff2'], ipr=self.ipr)[1])
        # int_grad2_outlier_Q = TA.binary_close_1d(self.row_stat['row_abs_diff2'] > self.int_grad2_th, self.sm_wd)
        # int_grad_outlier_Q = np.logical_or(int_grad2_outlier_Q, int_grad_outlier_Q)

        self.is_outlier = np.logical_and(corr_outlier_Q, int_grad_outlier_Q)
        self.outlier_ep_ind = TA.find_cc_1d(self.is_outlier)
        
    def select_outlier_intervals(self, max_dist_2_edge=None):
        if self.outlier_ep_ind is not None:
            int_length = self.outlier_ep_ind[:, 1] - self.outlier_ep_ind[:, 0]
            selected_Q = int_length <= self.max_artefact_width
            # Select outliers close to the edge of the image 
            if max_dist_2_edge is not None:
                num_row = self.im.shape[0]
                in_range_Q = np.logical_or(self.outlier_ep_ind < max_dist_2_edge, \
                                           self.outlier_ep_ind > (num_row - max_dist_2_edge))
                in_range_Q = in_range_Q.all(axis=1)
                selected_Q = np.logical_and(selected_Q, in_range_Q)

            # maybe need more here 
            if np.any(selected_Q):
                self.outlier_ep_ind = self.outlier_ep_ind[selected_Q]
                self.outlier_int_length = int_length[selected_Q]
            else:
                print(f"Outlier interval removed: {self.outlier_ep_ind}")
                self.outlier_ep_ind = None

        if self.outlier_ep_ind is not None:
            # Reject outliers
            self.num_outlier_int = self.outlier_ep_ind.shape[0]
            self.outlier_adj_cad_sum = np.array([np.sum(self.row_stat['adj_corr_da'][ep[0]:ep[1]]) for ep in self.outlier_ep_ind])
            self.max_oi_idx =  np.argmax(self.outlier_adj_cad_sum)
            self.max_oi_length = self.outlier_int_length[self.max_oi_idx]
            self.max_oi_ep_ind = self.outlier_ep_ind[self.max_oi_idx]
            self.num_row_shifted = self.im.shape[0] - self.max_oi_ep_ind[0]
        else:
            self.max_oi_ep_ind = None
            self.num_row_shifted = None

    def correct_row_shift(self):
        if self.max_oi_ep_ind is not None:
            row_shift = self.im.shape[0] - self.max_oi_ep_ind[0]
            self.im_corrected = np.roll(self.im, row_shift, axis=0)
        else:
            self.im_corrected = self.im

    def compute_edge_correlation(self):
        self.edge_corr = np.corrcoef(self.im[0], self.im[-1])[0,1]


    def vis_stat(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        # ax.plot(row_adj_corr_ad, label='AdjDiff', color='k')
        ax.plot(self.row_stat['adj_corr_da'], color='k')
        ax.axhline(y=self.corr_th, color='k')
        ax.set_ylabel("AdjCorrD")
        ax2 = ax.twinx()
        ax2.plot(self.row_stat['row_abs_diff'], color='g')
        ax2.axhline(y=self.int_grad_th, color='g')
        ax2.set_ylabel("RowIntAbsGrad", color='g')
        ax2.tick_params(axis="y", labelcolor='g')
        # ax2.set_ylim([0,1])
        ax.set_xlabel("Row index")
        if self.max_oi_ep_ind is not None:
            ax.axvline(x=self.max_oi_ep_ind[0], ymin=0, ymax=0.25)
            ax.axvline(x=self.max_oi_ep_ind[1], ymin=0, ymax=0.25)
        # ax.xlabel("Row index - 2")
        # ax.ylabel("Abs Diff Adj Corr")
        
        plt.show()
    
    def vis_im_with_detection(self, title=None):
        import matplotlib.pyplot as plt
        fig, [ax, ax2] = plt.subplots(ncols=2)
        ax.imshow(self.im)
        if self.max_oi_ep_ind is not None:
            ax.axhline(y=self.max_oi_ep_ind[0], color='r')
            ax.axhline(y=self.max_oi_ep_ind[1], color='r')
        ax.set_ylim((self.im.shape[0], 0))
        if title is not None:
            ax.set_title(title)
        self.correct_row_shift()
        ax2.imshow(self.im_corrected)
        ax2.set_title("Corrected image")
        plt.show()


    

