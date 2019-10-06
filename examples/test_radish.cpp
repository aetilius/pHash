/*

    pHash, the open source perceptual hash library
    Copyright (C) 2009 Aetilius, Inc.
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Evan Klinger - eklinger@phash.org
    David Starkweather - dstarkweather@phash.org

*/

#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "pHash.h"

using namespace std;

int main(int argc, char **argv) {
    const char *msg = ph_about();
    printf(" %s\n", msg);

    if (argc < 1) {
        printf("no input args\n");
        printf("expected: %s <image>\n", argv[0]);
        exit(1);
    }
    const char *img1 = argv[1];
    Digest digest;
    int ret = ph_image_digest(img1, 1.5, 3.5, digest, 180);
    printf("ret: %d\n", ret);

    for (int i = 0; i < digest.size; ++i) printf("%d ", digest.coeffs[i]);
    return 0;
}
