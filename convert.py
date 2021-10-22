import os
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from matplotlib.widgets import RangeSlider

from pydicom import dcmread
from pydicom.data import get_testdata_file

# fetch the path to the test data
path = get_testdata_file('DICOMDIR')
path = Path('CDVIEWER') / 'DICOMDIR'
ds = dcmread(path)
root_dir = Path(ds.filename).resolve().parent

print(f'Root directory: {root_dir}\n')

# Iterate through the PATIENT records
patients = ds.patient_records
print(f'Number of patients: {len(patients)}')
for patient in patients:
    print(f"PatientID {patient.PatientID}\nPatientName {patient.PatientName}")
    studies = [x for x in patient.children if x.DirectoryRecordType == "STUDY"]
    print(f'Number of studies: {len(studies)}')
    for study in studies:
        all_series = [ii for ii in study.children if ii.DirectoryRecordType == "SERIES"]
        print(f'Number of series: {len(all_series)}')
        for si, series in enumerate(all_series):
            print(f'Series {si + 1}/{len(all_series)}')
            images = [x for x in series.children if x.DirectoryRecordType == "IMAGE"]
            # plural = ('', 's')[len(images) > 1]
            # print(images)
            elems = [x["ReferencedFileID"] for x in images]
            paths = [[x.value] if x.VM == 1 else x.value for x in elems]
            paths = [Path(*p) for p in paths]
            print(f'Number of images: {len(paths)}')
            if len(paths) > 1:
                fig, axs = plt.subplots(1, 2)  # rows, cols
                # All axial, coronal, sagittal
                ar = np.array([dcmread(Path(root_dir) / x).pixel_array for x in paths])
                min_ar, max_ar = ar.min(), ar.max()
                print(min_ar, max_ar)
                ar = (ar - ar.min()) / (ar.max() - ar.min())  # normalize
                print(ar.min(), ar.max())
                v = 150
                # sub_ar, t, max_v = ar[v, :, :], 'axial', ar.shape[0]
                sub_ar, t, max_v = ar[:, v, :], 'coronal', ar.shape[1]
                # sub_ar, t, max_v = ar[:, :, v], 'sagittal', ar.shape[2]
                sub_ar_min, sub_ar_max = sub_ar.min(), sub_ar.max()
                axs[1].hist(sub_ar.flatten(), bins=100, color='k')
                low_color, high_color = (1, 0.618, 0), (0, 0.382, 1)
                mid_colors = [(0, 0, 0), (1, 1, 1)]
                min_flt, max_flt = 0.2, 0.8  # filters
                min_thr, max_thr = 0.3, 0.5  # thresholds
                img = axs[0].imshow(sub_ar, vmin=0, vmax=1)
                # Create the RangeSlider
                slider_ax = plt.axes([0.05, 0.05, 0.3, 0.03])
                thr_slider = RangeSlider(slider_ax, "thr",
                                         valinit=(min_thr, max_thr),
                                         valstep=0.001, valmin=0,
                                         valmax=1, color='k')
                flt_ax = plt.axes([0.05, 0.01, 0.3, 0.03])
                flt_slider = RangeSlider(flt_ax, "flt",
                                         valinit=(min_flt, max_flt),
                                         valstep=0.001, valmin=0,
                                         valmax=1, color='grey')
                # Create the Vertical lines on the histogram
                thr_lower_limit_line = axs[1].axvline(
                    thr_slider.val[0], color='k', linestyle='--')
                thr_upper_limit_line = axs[1].axvline(
                    thr_slider.val[1], color='k', linestyle='--')
                flt_lower_limit_line = axs[1].axvline(
                    flt_slider.val[0], color=low_color, linestyle='--')
                flt_upper_limit_line = axs[1].axvline(
                    flt_slider.val[1], color=high_color, linestyle='--')

                def update(*args, **kwargs):
                    # Update the image's colormap
                    min_thr, max_thr = thr_slider.val[0], thr_slider.val[1]
                    min_flt, max_flt = flt_slider.val[0], flt_slider.val[1]
                    if min_flt > min_thr:
                        min_thr = min_flt
                        thr_slider.set_val((min_flt, thr_slider.val[1]))
                    if max_flt < max_thr:
                        max_thr = max_flt
                        thr_slider.set_val((thr_slider.val[0], max_flt))
                    cm = LinearSegmentedColormap.from_list(
                        'rkwb', colors=[
                            (0., low_color),
                            (min_flt, low_color),
                            (min_thr, mid_colors[0]),
                            (max_thr, mid_colors[1]),
                            (max_flt, high_color),
                            (1., high_color)
                        ], N=400)
                    img.cmap = cm
                    # Update the position of the vertical lines
                    thr_lower_limit_line.set_xdata([min_thr, min_thr])
                    thr_upper_limit_line.set_xdata([max_thr, max_thr])
                    flt_lower_limit_line.set_xdata([min_flt, min_flt])
                    flt_upper_limit_line.set_xdata([max_flt, max_flt])
                    # Redraw the figure to ensure it updates
                    fig.suptitle(f'{t}, {v}/{max_v}, '
                                 f'global {min_ar} < v < {max_ar}, '
                                 f'local {sub_ar_min:.3f} < v < {sub_ar_max:.3f}, '
                                 f'{sub_ar_min * (max_ar - min_ar) + min_ar:.0f} < v < '
                                 f'{sub_ar_max * (max_ar - min_ar) + min_ar:.0f}, '
                                 f'threshold {min_thr} < v < {max_thr}, '
                                 f'{min_thr * (max_ar - min_ar) + min_ar:.0f} < v < '
                                 f'{max_thr * (max_ar - min_ar) + min_ar:.0f}')
                    fig.canvas.draw_idle()
                update()
                thr_slider.on_changed(update)
                flt_slider.on_changed(update)
                # fig.tight_layout()
                plt.show()
            # for p in paths:
            #     img = dcmread(Path(root_dir) / p)
            #     print()
