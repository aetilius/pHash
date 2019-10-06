#include "mman.h"

/* getpagesize for windows */
__declspec(dllexport) long getpagesize(void) {
    static long g_pagesize = 0;
    if (!g_pagesize) {
        SYSTEM_INFO system_info;
        GetSystemInfo(&system_info);
        g_pagesize = system_info.dwPageSize;
    }
    return g_pagesize;
}
__declspec(dllexport) long getregionsize(void) {
    static long g_regionsize = 0;
    if (!g_regionsize) {
        SYSTEM_INFO system_info;
        GetSystemInfo(&system_info);
        g_regionsize = system_info.dwAllocationGranularity;
    }
    return g_regionsize;
}

/* mmap for windows */
__declspec(dllexport) void *mmap(void *ptr, size_t size, int prot, int flags,
                                 int fd, char *fm_name, off_t offset) {
    long alloc_size = getregionsize();
    DWORD alloc_mask = ~(alloc_size - 1);
    DWORD alloc_offset_mask = alloc_size - 1;
    LPCTSTR fm_objname = (LPCTSTR)fm_name;
    HANDLE fh = (HANDLE)_get_osfhandle(fd);
    HANDLE fmaph = CreateFileMapping(fh, NULL, prot, 0, 0, fm_objname);

    if (fmaph == NULL) {
        DWORD err = GetLastError();
        fprintf(stderr,
                "mmap:unable to get file mapping object: errorcode %d\n", err);
        return NULL;
    }

    DWORD offs = (DWORD)offset;
    long alloc_unit = offs & alloc_mask;
    long alloc_offset = offs & alloc_offset_mask;

    DWORD offsetlow = alloc_unit;
    DWORD offsethigh = 0;

    char *buf = (char *)MapViewOfFile(fmaph, flags, offsethigh, offsetlow,
                                      alloc_offset + size);

    if (buf == NULL) {
        DWORD err = GetLastError();
        fprintf(stderr, "mmap:unable to create view of file: errorcode %d\n",
                err);
        return buf;
    }
    buf += alloc_offset;

    return buf;
}
/* munmap for windows */
__declspec(dllexport) long munmap(void *ptr, long size) {
    long res = 0;

    if (!UnmapViewOfFile(ptr)) res = MUNMAP_FAILURE;

    return res;
}

__declspec(dllexport) int ftruncate(int fd, off_t fp) {
    return _chsize(fd, fp);
}

__declspec(dllexport) int msync(void *start, size_t length, int mflags) {
    int res = 0;
    if (!FlushViewOfFile(start, length)) return -1;
    return res;
}
