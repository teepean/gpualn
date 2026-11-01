/* LD_PRELOAD shim for galn_pagecache.sh: turns posix_fadvise into a no-op, to measure G-ALN 0.2.6+
 * without the POSIX_FADV_DONTNEED call at the end of its load_bwt (bilalbioinfo/galn#2). */
#include <sys/types.h>
int posix_fadvise(int fd, off_t off, off_t len, int adv) { (void)fd; (void)off; (void)len; (void)adv; return 0; }
int posix_fadvise64(int fd, off_t off, off_t len, int adv) { (void)fd; (void)off; (void)len; (void)adv; return 0; }
