#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H

extern "C" {
    void  bl_mempool_init ();
    void  bl_mempool_finalize ();
    void* bl_mempool_alloc (size_t n);
    void  bl_mempool_free (void* p);
    void  bl_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  // min, max & tot in MB
    void  bl_real_array_init (double* p, size_t nelems); 
    void  bl_array_init_snan (double* p, size_t nelems);
}

#endif
