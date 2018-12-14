#ifndef   _ParSHUM_INTERNAL_MEM_H 
#define   _ParSHUM_INTERNAL_MEM_H


typedef struct _ParSHUM_internal_mem *ParSHUM_internal_mem;


typedef enum _ParSHUM_alloc {
  ParSHUM_malloc, 
  ParSHUM_aligned, 
  ParSHUM_calloc, 
} ParSHUM_alloc;


ParSHUM_internal_mem ParSHUM_internal_mem_create(size_t size, size_t alignment, ParSHUM_alloc type);
void                 ParSHUM_internal_mem_destroy(ParSHUM_internal_mem self);
size_t  ParSHUM_internal_mem_alloc(ParSHUM_internal_mem self, void **mem, size_t size);

#endif // _ParSHUM_INTERNAL_MEM_H
