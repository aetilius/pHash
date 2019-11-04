#include <errno.h>
#include <stdio.h>
#include "dirent.h"
#include "pHash.h"

using namespace std;

#define TRUE 1
#define FALSE 0

// data structure for a hash and id
struct ph_imagepoint {
    ulong64 hash;
    char *id;
};

// aux function to create imagepoint data struct
struct ph_imagepoint *ph_malloc_imagepoint() {
    return (struct ph_imagepoint *)malloc(sizeof(struct ph_imagepoint));
}

// auxiliary function for sorting list of hashes
bool cmp_lt_imp(struct ph_imagepoint dpa, struct ph_imagepoint dpb) {
    int result = strcmp(dpa.id, dpb.id);
    if (result < 0) return TRUE;
    return FALSE;
}
/** TEST for image DCT hash function
 *  The program reads all images from the two directories given on the command
 *line. The program expects the same number of image files in each directory.
 *Each image should be perceptually similar to a corresponding file in the other
 *directory and have the same exact name.  For example, one directory could
 *contain the originals, and the other directory blurred counterparts.  The
 *program calculates the hashes. First, the hamming distances are calculated
 *between all similar image hashes (-i.e. the intra compares), and then hamming
 *distances for different images hashes (-i.e. the inter compares).
 **/
int main(int argc, char **argv) {
    const char *msg = ph_about();
    printf(" %s\n", msg);

    if (argc < 3) {
        printf("no input args\n");
        printf("expected: \"test_imagephash [dir name] [dir_name]\"\n");
        exit(1);
    }
    const char *dir_name = argv[1];
    const char *dir_name2 = argv[2];
    printf("dir1: %s\n", dir_name);
    printf("dir2: %s\n", dir_name2);
    struct dirent *dir_entry;
    ph_imagepoint *dp = NULL;

    // first directory
    DIR *dir = opendir(dir_name);
    if (!dir) {
        printf("unable to open directory\n");
        exit(1);
    }
    errno = 0;
    int i = 0;
    ulong64 tmphash;
    char path[100];
    path[0] = '\0';
    while ((dir_entry = readdir(dir)) != 0) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, "..")) {
            strcat(path, dir_name);
            strcat(path, "/");
            strcat(path, dir_entry->d_name);
            printf("file[%d]=%s\n", i, path);
            i++;
        }
        errno = 0;
        path[0] = '\0';
    }

    if (errno) {
        printf("error reading directory\n");
        exit(1);
    }

    // second directory
    dir_entry = NULL;
    DIR *dir2 = opendir(dir_name2);
    if (!dir) {
        printf("unable to open directory\n");
        exit(1);
    }
    errno = 0;
    path[0] = '\0';
    i = 0;
    while ((dir_entry = readdir(dir2)) != 0) {
        if (strcmp(dir_entry->d_name, ".") && strcmp(dir_entry->d_name, "..")) {
            strcat(path, dir_name2);
            strcat(path, "/");
            strcat(path, dir_entry->d_name);
            printf("file[%d]=%s\n", i, path);
            i++;
        }
        errno = 0;
        path[0] = '\0';
    }

    if (errno) {
        printf("error reading directory\n");
        exit(1);
    }

    printf("done\n");
    return 0;
}
