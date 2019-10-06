#ifndef _MMAN_H_
#define _MMAN_H_

#include <Windows.h>
#include <fcntl.h>
#include <io.h>
#include <share.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#define MMAP_FAILURE -1
#define MUNMAP_FAILURE -1
#define MS_SYNC 0

/* determine system page size */
__declspec(dllexport) long getpagesize(void);

/* determine system memory allocation size */
__declspec(dllexport) long getregionsize(void);

/* !memory map a file
 *  @param ptr - void ptr - NULL
 *  @param size - size_t size to map
 *  @param prot - int protection flags (flags passed to CreateFileMapping in
 * windows) available flags: PAGE_READONLY, PAGE_READWRITE, PAGE_WRITECOPY
 *  @param flags - int extra flags (flags passed to MapView of File in windows)
 *                available flags: FILE_MAP_ALL_ACCESS, FILE_MAP_COPY,
 * FILE_MAP_READ, FILE_MAP_FILE_MAP_WRITE
 *  @param fd - int file descriptor
 *  @return void ptr to mapped memory - null for error
 */
__declspec(dllexport) void *mmap(void *ptr, size_t size, int prot, int flags,
                                 int fd, char *fm_name, off_t offset);

/* !unmap memory
 * @param ptr - ptr to memory to release
 * @param size - long size of memory block to free (ignored)
 * @return long value 0 for success, -1 (MUNMAP_FAILURE) for error
 */
__declspec(dllexport) long munmap(void *ptr, long size);

/* !truncate file to set length
 *  @param fd - int file descriptor
 *  @param fp - off_t file pos at which to truncate the file
 *  @return int  -  0 for success, -1 for error (
 *  (errno is set to EACCESS if file is locked against access, EBADF if file is
 * readonly or handle is invalid, ENOSPC if not enough space left on the device)
 */
__declspec(dllexport) int ftruncate(int fd, off_t fp);

/* !synchronize the file with persistent memory
 * @param start - ptr to start of memory to sync
 * @param length - size_t number of bytes to flush
 * @param flags - int msync flags (ignored)
 * @return int value - 0 for success, -1 for error (use GetLastError() to get
 * extended information)
 */
__declspec(dllexport) int msync(void *start, size_t length, int flags);

#endif
