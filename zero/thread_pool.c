#include <gkyl_thread_pool.h>
#include <gkyl_alloc.h>

#include <thpool.h>

struct jp_thread_pool {
  struct gkyl_job_pool jp; // base job-pool object
  threadpool thpool; // thread-pool object
};

static void
thread_pool_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_job_pool *base = container_of(ref, struct gkyl_job_pool, ref_count);
  struct jp_thread_pool *th = container_of(base, struct jp_thread_pool, jp);
  thpool_destroy(th->thpool);
  gkyl_free(th);
}

struct gkyl_job_pool*
gkyl_thread_pool_new(int nthreads)
{
  struct jp_thread_pool *th = gkyl_malloc(sizeof(struct jp_thread_pool));

  th->jp.pool_size = nthreads;
  th->thpool = thpool_init(nthreads);

  // set reference counter
  th->jp.ref_count = (struct gkyl_ref_count) { thread_pool_free, 1 };
    
  return &th->jp;
}
