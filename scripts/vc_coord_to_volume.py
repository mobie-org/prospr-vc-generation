#! /g/kreshuk/zinchenk/miniconda3/envs/platyneris/bin/python

import csv
import os
import subprocess
import argparse
import numpy as np
from vigra.analysis import extractRegionFeatures, labelVolumeWithBackground
from vigra.filters import distanceTransform
from vigra.impex import writeVolume


def map_coordinates_to_volume(coord_file):
    # let this be hardcoded
    volume_shape = (251, 550, 320)
    prospr_volume = np.zeros(volume_shape, dtype='int16')
    with open(coord_file) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        label = 0
        for line in tsv_reader:
            label += 1
            is_muscle = line[0].startswith('MHCL4')
            line_split = line[1].split(';')
            for coord_string in line_split:
                if not coord_string:
                    continue
                coord = []
                for item in coord_string.replace(',', ' ').split(' '):
                    if item.isdigit():
                        coord.append(int(item))
                assert len(coord) == 3, coord_string
                # in tsv the coordinates are stored as xyz
                if prospr_volume[coord[2]][coord[0]][coord[1]] != 0:
                    if is_muscle:
                        continue
                prospr_volume[coord[2]][coord[0]][coord[1]] = label
    return prospr_volume


def upsample_the_volume(prospr_volume):
    # no, I couldn't just dilate, we have labels
    for axis in range(3):
        prospr_volume = np.moveaxis(prospr_volume, axis, 0)
        for x in range(2, prospr_volume.shape[0], 3):
            # copy the middle slice to both sides for all dims
            prospr_volume[x-1] = prospr_volume[x-1] + prospr_volume[x]
            prospr_volume[x+1] = prospr_volume[x+1] + prospr_volume[x]
        prospr_volume = np.moveaxis(prospr_volume, 0, axis)
    return prospr_volume


def get_features(data):
    # compute the relevant vigra region features
    features = extractRegionFeatures(data.astype('float32'), data.astype('uint32'),
                                     features=['Coord<Maximum >', 'Coord<Minimum >',
                                               'Count', 'RegionCenter'])
    # compute bounding boxes from features
    mins = features['Coord<Minimum >'].astype('uint32')
    maxs = features['Coord<Maximum >'].astype('uint32') + 1
    bbs = [tuple(slice(mi, ma) for mi, ma in zip(min_, max_))
           for min_, max_ in zip(mins, maxs)]
    sizes = features['Count'].squeeze().astype('uint64')
    mass_centers = np.rint(features['RegionCenter']).astype("uint32")
    return bbs, sizes, mass_centers


def remove_small_clusters(cluster_data, size_thres):
    sv_len = 3  # I guess this value will not change
    size_thres = size_thres * sv_len**3
    cluster_bbs, cluster_sizes, _ = get_features(cluster_data)
    labels_to_remove = np.where(cluster_sizes < size_thres)[0]
    for label in labels_to_remove:
        bb = cluster_bbs[label]
        mask = cluster_data[bb] == label
        cluster_data[bb][mask] = 0


def remove_vc_labels(labels_to_remove, cluster_data, vc_data, vc_bbs):
    assert cluster_data.shape == vc_data.shape
    for label in labels_to_remove:
        bb = vc_bbs[label]
        mask = vc_data[bb] == label
        cluster_data[bb][mask] = 0
        vc_data[bb][mask] = 0


def map_vc_to_cluster_labels(cluster_data, vc_data, vc_bbs):
    assert cluster_data.shape == vc_data.shape
    vc_label_dict = {}
    vc_labels = np.unique(vc_data)
    vc_labels = vc_labels[vc_labels != 0]
    for vc_label in vc_labels:
        vc_bb = vc_bbs[vc_label]
        vc_mask = vc_data[vc_bb] == vc_label
        cluster_label = np.unique(cluster_data[vc_bb][vc_mask])
        assert len(cluster_label) == 1, cluster_label
        vc_label_dict[vc_label] = cluster_label[0]
    return vc_label_dict


def filter_tiny_nonsymmetric(data, min_vc_to_keep=8, radius=4, min_vc=5):
    sv_len = 3  # I guess this value will not change
    radius = radius * sv_len
    min_vc = min_vc * sv_len**3
    min_vc_to_keep = min_vc_to_keep * sv_len**3
    pad_data = np.pad(data, (radius, radius), 'constant')
    vc_split = labelVolumeWithBackground(pad_data.astype('uint32'), neighborhood=26)
    vc_bbs, vc_sizes, vc_centers = get_features(vc_split)
    vc_to_clust_labels = map_vc_to_cluster_labels(pad_data, vc_split, vc_bbs)
    noise_to_remove = list(np.where(vc_sizes < min_vc)[0])
    remove_vc_labels(noise_to_remove, pad_data, vc_split, vc_bbs)
    tiny_vc_labels = np.where((vc_sizes >= min_vc) & (vc_sizes < min_vc_to_keep))[0]
    to_remove = []
    for vc_label in tiny_vc_labels:
        # get the cluster label of this vc
        cluster_label = vc_to_clust_labels[vc_label]

        # get a mask-sphere with a given radius
        center = vc_centers[vc_label]
        mask = np.zeros([radius*2+1] * 3, dtype='uint32')
        mask[(radius,)*3] = 1  # mark the center
        dist = distanceTransform(mask)
        sphere_mask = (dist <= radius)

        # get a bb coordinates for the symm side
        bb_sphere = [slice(i - radius, i + radius + 1) for i in center]
        bb_sphere[2] = slice(320 + 2*radius - bb_sphere[2].stop,
                             320 + 2*radius - bb_sphere[2].start)
        bb_sphere = tuple(bb_sphere)

        # find labels within a sphere mask on the symm side
        sym_labels = np.unique(pad_data[bb_sphere][sphere_mask])
        # if there's no symm vc with the same cluster label, remove
        if cluster_label not in sym_labels:
            to_remove.append(vc_label)
    remove_vc_labels(to_remove, pad_data, vc_split, vc_bbs)
    print("Removed", len(noise_to_remove) + len(to_remove), "pieces")
    print(len(np.unique(vc_split)), 'pieces left')
    pad_data = pad_data[(slice(radius, -radius),) * 3]
    vc_split = vc_split[(slice(radius, -radius),) * 3]
    return pad_data, vc_split, vc_to_clust_labels


def remap_to_dense_labels(data):
    labels = np.unique(data)
    new_map_dict = {label: i for i, label in enumerate(labels)}
    new_data = np.array([[[new_map_dict.get(label, 0) for label in y] for y in x] for x in data])
    return new_data


def save_coord_file(data, old_file, new_file, label_dict=None):
    # to save the vc labels provide the label_dict
    with open(old_file) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        vc_names = [line[0] for line in tsv_reader]
    # to keep sequential numbering
    keep_counts = np.ones(np.max(len(label_dict.values())), dtype='int') if label_dict else None
    with open(new_file, 'w') as tsv_file:
        csv_writer = csv.writer(tsv_file, delimiter='\t')
        for label in np.unique(data):
            if label == 0:
                continue
            coords = np.where(data == label)
            # take only the central pixel of the supervoxel, and change to xyz
            centr_coords = [str((x, y, z)) for z, x, y in zip(*coords)
                            if np.all(np.array((z, x, y)) % 3 == 2)]
            row = ' ; '.join(centr_coords).replace('(', '( ').replace(')', ' )')
            if label_dict:
                cluster_label = label_dict[label]
                name = vc_names[cluster_label-1] + '_' + str(keep_counts[cluster_label])
                keep_counts[cluster_label] += 1
            else:
                name = vc_names[label-1]
            _ = csv_writer.writerow([name, row])


def save_tif(data, out_path, resolution):
    # the courtesy of Constantin
    # write initial tif with vigra
    out_path = out_path + '.tif'
    writeVolume(data.T.astype('float'), out_path, '', dtype='FLOAT')
    # encode the arguments for the imagej macro:
    # imagej macros can only take a single string as argument, so we need
    # to comma seperate the individual arguments
    assert "," not in out_path, "Can't encode pathname containing a comma"
    arguments = "%s,%i,%f,%f,%f" % (out_path, data.shape[0],
                                    resolution[0], resolution[1], resolution[2])
    # call the imagej macro
    fiji_executable = "/g/kreshuk/zinchenk/Fiji.app/ImageJ-linux64"
    script = "/g/kreshuk/zinchenk/cell_match/scripts/other/set_voxel_size.ijm"
    cmd = [fiji_executable, '-batch', '--headless', script, arguments]
    subprocess.run(cmd)


#coord_file = '/g/kreshuk/zinchenk/cell_match/data/genes/CellModels_ALL_coordinates.tsv'
#out_tif = '/g/kreshuk/zinchenk/cell_match/data/genes/vc_volume_prospr_space_all_vc15_ms8_nt4'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert vc coordinates tsv into a volume with labels.')
    parser.add_argument('coord_file', type=str, help='path to the coordinates (tsv) file')
    parser.add_argument('output_file', type=str, help='path to the output tiff file (no extension required)')
    parser.add_argument('--min_clust_size', type=int, default=16,
                        help='all clusters below this size (in supervoxels) will be removed')
    parser.add_argument('--min_vc_size', type=int, default=5,
                        help='all the VCs below this size (in supervoxels) will be removed')
    parser.add_argument('--min_vc_size_to_keep', type=int, default=8,
                        help='all the VCs below this size (in supervoxels) will be checked for symmetry')
    parser.add_argument('--search_radius', type=int, default=4,
                        help='the search radius (in supervoxels) when looking for a symmetric VC')
    parser.add_argument('--save_vcs', type=bool, default=False,
                        help='if True saves the coordinates of VCs, if False - of clusters')
    args = parser.parse_args()

    # the new coordinate file after VC filtering
    tag = '_vc_curated' if args.save_vcs else '_clust_curated'
    new_coord_file = os.path.splitext(args.coord_file)[0] + tag + os.path.splitext(args.coord_file)[1]
    # making the pixel coordinates from tsv and labeling the pixel in the volume with the row number
    volume = map_coordinates_to_volume(args.coord_file)
    # we need to upsample, because the coordinates in tsv file are for the central pixel of 3*3 cube
    volume = upsample_the_volume(volume)
    remove_small_clusters(volume, args.min_clust_size)
    curated_clust_volume, curated_vc_volume, vc_to_clust = \
        filter_tiny_nonsymmetric(volume, min_vc=args.min_vc_size, radius=args.search_radius,
                                 min_vc_to_keep=args.min_vc_size_to_keep)
    if args.save_vcs:
        save_coord_file(curated_vc_volume, args.coord_file, new_coord_file, vc_to_clust)
        remapped_volume = remap_to_dense_labels(curated_vc_volume)
    else:
        save_coord_file(curated_clust_volume, args.coord_file, new_coord_file)
        remapped_volume = remap_to_dense_labels(curated_clust_volume)
    # flip the volume horizontally
    flipped_volume = remapped_volume[:, :, ::-1].copy()

    save_tif(flipped_volume, args.output_file, (0.55,)*3)
