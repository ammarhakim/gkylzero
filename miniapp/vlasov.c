#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>

struct gkyl_vlasov_app {
};

gkyl_vlasov_app*
vlasov_app_new(struct gkyl_vm vm)
{
  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));

  return app;
}

void
vlasov_app_release(gkyl_vlasov_app* app)
{
  gkyl_free(app);
}
