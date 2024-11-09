import numpy as np
import scipy as sp
import cv2
from scipy.io import savemat
from scipy import ndimage
# import pandas as pd
from scipy import interpolate

def compute_basic_statistics(data, bins=None):
    # bins:
    #   If scalar, should be an integer that specifies the number of bins used for the histogram calculation.
    #   If a list or numpy array, the PDF, CDF calculation will be based on the bin edge. Data outside the edge will be ignored.
    # Takes 5.7 seconds for data.shape = (150, 1280, 1024)
    if not isinstance(data, np.ndarray):
        data = np.array(data)
    stat = {}
    stat['num_data'] = data.size
    stat['mean'] = np.mean(data)
    mean_x2 = np.mean(data ** 2)
    stat['std'] = np.sqrt(mean_x2 - stat['mean'] ** 2)
    stat['cv'] = stat['std'] / stat['mean']

    if bins is None:
        stat['num_bins'] = 10
    elif np.isscalar(bins):
        stat['num_bins'] = bins
    else:
        if not isinstance(bins, np.ndarray):
            bins = np.array(bins)

    if isinstance(bins, np.ndarray):
        stat['num_bins'] = bins.size - 1
        stat['hist_count'], stat['bin_edge'] = np.histogram(data, bins=bins)
    else:
        stat['hist_count'], stat['bin_edge'] = np.histogram(
            data, bins=stat['num_bins'])

    stat['bin_val'] = (stat['bin_edge'][:-1] + stat['bin_edge'][1:]) / 2
    stat['bin_width'] = stat['bin_edge'][1:] - stat['bin_edge'][:-1]
    stat['probability'] = stat['hist_count'] / np.sum(stat['hist_count'])
    stat['pdf'] = stat['probability'] / stat['bin_width']
    stat['cdf'] = np.cumsum(stat['probability'])

    stat['prctile_th'] = np.array(
        [0, 0.1, 1, 5, 10, 25, 50, 75, 90, 95, 99, 99.9, 100])
    stat['prctile_val'] = np.percentile(data, q=stat['prctile_th'])
    stat['max'] = stat['prctile_val'][-1]
    stat['min'] = stat['prctile_val'][0]
    return stat


def compute_three_view_mip(data):
    # data is a 3-dimension array in the order of zyx
    mip = {}
    mip['yx'] = np.amax(data, axis=0)
    mip['zx'] = np.amax(data, axis=1)
    mip['zy'] = np.amax(data, axis=2)
    return mip


def comptue_z_statistics(data):
    # Take 1.4 second for a (150, 1280, 1024) array
    z_stat = {}
    z_stat['mean'] = np.mean(data, axis=(1, 2))
    z_stat['std'] = np.std(data, axis=(1, 2))
    # z_stat['']
    return z_stat


def normalized_integer_data(data):
    # This function return a new array and does not alter the input outside the function
    assert isinstance(
        data, np.ndarray), 'The input data should be a numpy array'
    data_type = data.dtype
    if np.issubdtype(data_type, np.integer):
        type_info = np.iinfo(data_type)
        data = data.astype(np.float32) / type_info.max
    return data


def save_dict_as_mat_file(fp, data):
    savemat(fp, data, appendmat=True, format='5',
            long_field_names=True, do_compression=False, oned_as='column')


def im_int16_to_uint16(data, method='truncate'):
    # Set all negative values to 0
    assert data.dtype == np.int16, 'The input data is suppoosed to be an int16 array'
    if method == 'truncate':
        data = np.maximum(0, data).astype(np.uint16)
    elif method == 'shift':
        data = (data.astype(np.float32) +
                np.iinfo(np.int16).max).astype(np.uint16)

    return data


def compute_z_fraction_above_th(data, th):
    assert isinstance(
        data, np.ndarray), 'The input data is supposed to be an numpy array'
    assert np.isscalar(th), 'The threshold should be a numeric scalar'
    fraction = np.ones(data.shape[0])
    for i in range(data.shape[0]):
        fraction[i] = np.count_nonzero(data[i] >= th) / np.size(data[i])
    return fraction


def compute_step_mip(data, step):
    step = int(step)
    if data.ndim == 2:
        return data[None, :, :]
    elif data.ndim == 3:
        if step > 1:
            num_sec = data.shape[0]
            num_step = np.ceil(num_sec / step).astype(int)
            ch_mip = np.zeros((num_step, data.shape[1], data.shape[2]),
                              dtype=data.dtype)
            for i in range(num_step):
                start_idx = i * step
                end_idx = min((i+1) * step, num_sec)
                ch_mip[i, :, :] = np.amax(
                    data[start_idx:end_idx, :, :], axis=0)
            return ch_mip
        elif step == 1:
            return data
        else:
            # To be implemented. Interpolation?
            raise ValueError('step size should be >= 1')


def correct_global_lineshift(im, shift):
    if shift != 0:
        if im.ndim == 2:
            im[1::2, :] = np.roll(im[1::2, :], shift, axis=1)
        elif im.ndim == 3:
            for i in range(im.shape[0]):
                im[i, 1::2, :] = np.roll(im[i, 1::2, :], shift, axis=1)
    return im


def downsampling_by_stepping(im, step):
    assert im.ndim == np.size(step), 'Dimension mismatch'
    if im.ndim == 2:
        return im[::step[0], ::step[1]]
    elif im.ndim == 3:
        return im[::step[0], ::step[1], ::step[2]]


def compute_global_column_shift(im, search_radius=30):
    assert type(
        im).__module__ == np.__name__, 'The first input should be a numpy array'
    assert im.ndim == 2, 'The first input should be of dimension 2'
    # Rescale image
    if np.issubdtype(im.dtype, np.integer):
        im = (im / np.iinfo(im.dtype).max).astype(np.float16)
    else:
        im = ((im - im.mean()) / (im.max() - im.min())).astype(np.float16)
    im_1 = im[::2]
    im_2 = im[1::2]
    search_list = np.arange(-search_radius, search_radius + 1)
    corr_list = np.zeros(search_list.shape)
    for i, shift in enumerate(search_list):
        corr_list[i] = np.mean(im_1 * np.roll(im_2, shift, axis=1))
        # If tier, use the last one
    max_corr_idx = search_list[np.argmax(corr_list)]
    # Erhan's implementation uses interpolation
    # itp_x = np.linspace(search_list[0], search_list[-1], num=1000, endpoint=True)
    # itp_f = interp1d(search_list, corr_list, kind='cubic')
    # max_corr_idx_f = itp_x[np.argmax(itp_f(itp_x))]
    # max_corr_idx = int(np.round(max_corr_idx_f))
    # return max_corr_idx, [search_list, corr_list]
    return max_corr_idx


def compute_single_image_y_gradient_peak(im):
    im_grad = np.diff(im, axis=0)
    abs_avg_grad = np.abs(im_grad.mean(axis=1))
    grad_prctile = np.percentile(abs_avg_grad, [25, 50, 75])
    dev2int = (abs_avg_grad - grad_prctile[1]) / \
        (grad_prctile[2] - grad_prctile[0])
    return dev2int.max, dev2int.argmax


def remove_shot_noise_2d_in_plane(data, wd_sz):
    data_dim = data.ndim
    if data_dim == 2:
        data = data[None, :, :]
    for i in range(data.shape[0]):
        data[i, :, :] = cv2.medianBlur(data[i, :, :], wd_sz)
        # data[i, :, :] = ndimage.median_filter(data[i, :, :], size=wd_sz)

    if data_dim == 2:
        data = np.squeeze(data)
    return data


def imgaussfilt2(data, wd_sz):
    data_dim = data.ndim
    if data_dim == 2:
        data = data[None, :, :]
    for i in range(data.shape[0]):
        data[i, :, :] = cv2.GaussianBlur(data[i, :, :], wd_sz, 0)
    if data_dim == 2:
        data = np.squeeze(data)
    return data


def resize_image_stack_2d(im, target_size_2d):
    if im.ndim == 2:
        im = im[None, :, :]
    target_stack_size = (im.shape[0], target_size_2d[0], target_size_2d[1])
    im_rz = np.zeros(target_stack_size, dtype=im.dtype)
    for i in range(im.shape[0]):
        # The size in opencv is (width, height)
        im_rz[i, :, :] = cv2.resize(
            im[i, :, :], (target_size_2d[1], target_size_2d[0]), interpolation=cv2.INTER_AREA)
    return im_rz

def row_shift_detection(tmp_im, min_th=0.05, sm_wd_sz=3):
    
    row_mean = np.mean(tmp_im, axis=1, keepdims=True)
    # This need to be updated later. We might need to change the acquisition format 
    row_mean_n = row_mean / np.iinfo(np.int16).max 
    row2_mean = np.mean(tmp_im.astype(np.double) ** 2, axis=1, keepdims=True)
    row_std = np.sqrt(row2_mean - row_mean**2)
    row_demean = tmp_im.astype(np.double) - row_mean
    # Compute the row-adjacent correlation value 
    row_adj_corr = np.mean(row_demean[0:-1, :] * row_demean[1:, :],
                        axis=1, keepdims=True) / (row_std[0:-1, :] * row_std[1:, :])
    row_adj_corr_ad = np.abs(np.diff(row_adj_corr, axis=0)) 
    # Mean-weighted adjacent correlation 
    adj_corr_ad_mw = row_adj_corr_ad * row_mean_n[1:-1]
    # Outlier detection 
    p25, p50, p75 = np.percentile(adj_corr_ad_mw, [25, 50, 75])
    pth = max(p50 + 3 * (p75 - p25), min_th)
    is_outlier = adj_corr_ad_mw > pth
    is_outlier = ndimage.binary_closing(is_outlier.flatten(), np.ones((sm_wd_sz, )))
    # Todo: find the most prominant one 
    # Analyze the distance to the button of the frame 
    if np.any(is_outlier):
        outlier_ep_Q = np.concatenate(([False], is_outlier, [False]), axis=None)
        outlier_ep_Q = np.diff(outlier_ep_Q, )
        outlier_ep_idx = np.nonzero(outlier_ep_Q)
    else:
        outlier_ep_idx = None
    
    return outlier_ep_idx
    

def compute_row_stat(im, visQ=False):
    im = np.double(im)
    result = {}
    result['row_mean'] = np.mean(im, axis=1)
    row2_mean = np.mean(im ** 2, axis=1)
    std = np.sqrt(row2_mean - result['row_mean'] ** 2)
    std[std == 0] = np.nan
    std = fill_edge_nan_1d(std, method='nearest')
    std = fill_internal_nan(std, method='linear')
    result['std'] = std
    im_dm = im - result['row_mean'].reshape((-1, 1))
    
    # Adjacent row correlation 
    adj_corr = np.mean(im_dm[0:-1, :] * im_dm[1:, :], axis=1) / \
        (result['std'][0:-1] * result['std'][1:])
    adj_corr = np.concatenate(([adj_corr[0]], np.maximum(adj_corr[1:], adj_corr[:-1]), [adj_corr[-1]]))
    # adj_corr = np.concatenate(([adj_corr[0]], (adj_corr[1:] + adj_corr[:-1]) / 2, [adj_corr[-1]]))
    result['adj_corr'] = adj_corr
    result['adj_corr_d'] = np.gradient(adj_corr)
    result['adj_corr_da'] = np.abs(result['adj_corr_d'])
    
    # Average adjacent row intensity gradient absolute difference 
    im_grad = np.abs(im[1:, :] - im[0:-1, :])
    im_grad = np.concatenate((im_grad[0].reshape((1, -1)), (im_grad[0:-1, :] + im_grad[1:, ])/2, im_grad[-1].reshape((1, -1))), axis=0)
    im_grad2 = np.abs(im_grad[1:, :] - im_grad[0:-1, :])
    im_grad2 = np.concatenate((im_grad2[0].reshape((1, -1)), (im_grad2[0:-1, :] + im_grad2[1:, ])/2, im_grad2[-1].reshape((1, -1))), axis=0)

    # im_grad = np.concatenate((im_grad[0].reshape((1, -1)), np.maximum(im_grad[0:-1, :], im_grad[1:, ]), im_grad[-1].reshape((1, -1))), axis=0)
    result['row_abs_diff'] = np.mean(im_grad, axis=1)
    result['row_abs_diff2'] = np.mean(im_grad2, axis=1)
    # result['row_abs_diff'] = np.mean(np.abs(np.gradient(im, axis=0)), axis=1)

    return result

def find_cc_1d(mask_1d):
    # Find connected component in 1D logical array 
    if np.any(mask_1d):
        outlier_ep_Q = np.concatenate(([False], mask_1d, [False]), axis=None)
        outlier_ep_Q = np.diff(outlier_ep_Q)
        outlier_ep_idx = np.nonzero(outlier_ep_Q)[0]
        outlier_ep_idx = outlier_ep_idx.reshape((-1, 2))
        return outlier_ep_idx
    else:
        return None    


def compute_percentile_outlier_threshold(x, ipr = 3):
    p25, p50, p75 = np.nanpercentile(x, [25, 50, 75])
    width = ipr * (p75 - p25)
    return [p50 - width, p50 + width]

# def compute_moving_percentile(data, wd_sz, percentile):
#     percentile = float(percentile) / 100
#     data = pd.DataFrame(data)
#     wd_mp = round(wd_sz / 2)
#     data = data.rolling(window=wd_sz, center=True).quantile(percentile)
#     data = data.to_numpy().flatten()
#     # Replace the boundary nan to the nearest non nan value 
#     data = fill_edge_nan_1d(data, method='nearest')
#     return data


def fill_internal_nan(data, method='linear'):
    mask = np.isnan(data)
    if np.any(mask) and not np.all(mask):
        itp_hdl = interpolate.interp1d(np.flatnonzero(~mask), data[~mask], kind=method)
        data[mask] = itp_hdl(np.flatnonzero(mask))
    return data

def fill_edge_nan_1d(data, method='nearest'):
    if method == 'nearest':
        ind = np.flatnonzero(~np.isnan(data))
        if len(ind):
            data[:ind[0]] = data[ind[0]]
            data[ind[-1]:] = data[ind[-1]]

    return data 

def binary_close_1d(x, wd_sz):
    # Scipy implementation of binary close removes the connected elements near the edge (smaller than the window size)
    return np.logical_or(x, ndimage.binary_closing(x, np.ones(wd_sz, )))


if __name__ == '__main__':
    test_data = np.arange(101)
    data_stst = compute_basic_statistics(test_data, np.arange(0, 100, 20))
    print(data_stst)
