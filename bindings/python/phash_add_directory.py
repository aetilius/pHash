"""
Hash all images in the specified directory
"""
import os
import pickle
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from phash import DCTImageHash, MHImageHash


def dct_hash_image(path):
    print(path, flush=True)
    return DCTImageHash.from_path(Path(path))


def mh_hash_image(path):
    print(path, flush=True)
    return MHImageHash.from_path(Path(path))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('directory')
    p.add_argument('cache_file', default='phash.cache', nargs='?')
    p.add_argument('-e', '--extension', default='.jpg')
    p.add_argument('-m', '--hash-method', default='dct', choices=('dct', 'mh'))
    return p.parse_args()


def main():
    opts = parse_args()
    if not os.path.isdir(opts.directory):
        raise RuntimeError('directory %s does not exist' % opts.directory)
    image_files = [os.path.join(opts.directory, p) for p in os.listdir(opts.directory) if p.endswith(opts.extension)]
    with Pool(processes=cpu_count()) as p:
        if opts.hash_method == 'dct':
            hashes = p.map(dct_hash_image, image_files)
        elif opts.hash_method == 'mh':
            hashes = p.map(mh_hash_image, image_files)
        else:
            raise RuntimeError('unknown hash method %s' % opts.hash_method)
    if os.path.exists(opts.cache_file):
        with open(opts.cache_file, 'rb') as fin:
            image_files_in, hashes_in = pickle.load(fin)
            image_files += image_files_in
            hashes += hashes_in
    with open(opts.cache_file, 'wb') as fout:
        pickle.dump((image_files, hashes), fout, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
