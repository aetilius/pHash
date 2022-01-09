"""
Print a list of similar pictures
optionally delete it (maybe in a second pass)
"""
import os
import sys
import argparse
import pickle


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('cache_file')
    p.add_argument('-t', '--threshold', type=int, default=10)
    p.add_argument('-d', '--delete-threshold', type=float)
    return p.parse_args()


def main():
    opts = parse_args()
    with open(opts.cache_file, 'rb') as fin:
        image_files, hashes = pickle.load(fin)
    print('hashed', len(hashes), 'files with', hashes[0].__class__.__name__, file=sys.stderr)

    print('<title>', opts.cache_file, '</title>')
    print('<table border=1>')
    print('<tr>')
    print('<td></td>')
    for i in range(opts.threshold):
        print('<th>%d</th>' % i)
    print('</tr>')
    for i1, h1 in enumerate(hashes):
        similar = {}
        for i2, h2 in enumerate(hashes[i1+1:], start=i1+1):
            distance = h1.hamming_distance(h2)
            if opts.delete_threshold is not None and distance < opts.delete_threshold:
                print('deleting file', image_files[i2], file=sys.stderr)
                os.unlink(image_files[i2])
            elif distance < opts.threshold:
                similar.setdefault(distance, []).append(i2)
        if len(similar):
            print('<tr>')
            print('<th><a href="%s">[%03d]</a></th>' % (image_files[i1], i1))
            for i in range(opts.threshold):
                if i in similar:
                    print('<td>')
                    for j in similar[i]:
                        print('<a href="%s">[%03d]</a>' % (image_files[j], j))
                    print('</td>')
                else:
                    print('<td></td>')
            print('</tr>')
    print('</table>')


if __name__ == '__main__':
    main()
