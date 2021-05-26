
# You can set compiler to use as:
#
# make CC=mpicc 
#

CFLAGS = -O3 -g -Iminus -Izero -Iapps -Ikernels/basis -Ikernels/maxwell -Ikernels/vlasov
PREFIX = ${HOME}/gkylsoft

# Header dependencies
headers = $(wildcard minus/*.h) $(wildcard zero/*.h) $(wildcard apps/*.h) $(wildcard kernels/*/*.h)

# Object files to compile in library
libobjs = $(patsubst %.c,%.o,$(wildcard minus/*.c)) \
	$(patsubst %.c,%.o,$(wildcard zero/*.c)) \
	$(patsubst %.c,%.o,$(wildcard apps/*.c)) \
	$(patsubst %.c,%.o,$(wildcard kernels/*/*.c))

# Make targets: libraties, regression tests and unit tests
all: makeout/libgkylzero.a \
	$(patsubst %.c,makeout/%,$(wildcard regression/app_*.c)) \
	$(patsubst %.c,makeout/%,$(wildcard unit/ctest_*.c))

# Library archive
makeout/libgkylzero.a: ${libobjs}
	ar -crs makeout/libgkylzero.a ${libobjs}

# Regression tests
makeout/regression/twostream.ini: regression/twostream.ini
	cp regression/twostream.ini makeout/regression/twostream.ini

makeout/regression/app_5m_gem: regression/app_5m_gem.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_5m_gem.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_5m_gem

makeout/regression/app_5m_riem: regression/app_5m_riem.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_5m_riem.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_5m_riem

makeout/regression/app_esshock: regression/app_esshock.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_esshock.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_esshock

makeout/regression/app_euler_multiblock: regression/app_euler_multiblock.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_euler_multiblock.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_euler_multiblock

makeout/regression/app_euler_reflect_2d: regression/app_euler_reflect_2d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_euler_reflect_2d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_euler_reflect_2d

makeout/regression/app_euler_riem_2d: regression/app_euler_riem_2d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_euler_riem_2d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_euler_riem_2d

makeout/regression/app_euler_sodshock: regression/app_euler_sodshock.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_euler_sodshock.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_euler_sodshock

makeout/regression/app_iso_euler_sodshock: regression/app_iso_euler_sodshock.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_iso_euler_sodshock.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_iso_euler_sodshock

makeout/regression/app_maxwell_plane_wave_1d: regression/app_maxwell_plane_wave_1d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_maxwell_plane_wave_1d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_maxwell_plane_wave_1d

makeout/regression/app_maxwell_plane_wave_2d: regression/app_maxwell_plane_wave_2d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_maxwell_plane_wave_2d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_maxwell_plane_wave_2d

makeout/regression/app_maxwell_reflect_2d: regression/app_maxwell_reflect_2d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_maxwell_reflect_2d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_maxwell_reflect_2d

makeout/regression/app_pthread_proj: regression/app_pthread_proj.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_pthread_proj.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_pthread_proj

makeout/regression/app_twostream: regression/app_twostream.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_twostream.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_twostream

makeout/regression/app_vlasov_kerntm: regression/app_vlasov_kerntm.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_vlasov_kerntm.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_vlasov_kerntm

makeout/regression/app_weibel_3d: regression/app_weibel_3d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_weibel_3d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_weibel_3d

makeout/regression/app_weibel_4d: regression/app_weibel_4d.c makeout/libgkylzero.a
	${CC} ${CFLAGS} regression/app_weibel_4d.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/regression/app_weibel_4d

# Unit tests

makeout/unit/ctest_alloc: unit/ctest_alloc.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_alloc.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_alloc

makeout/unit/ctest_array: unit/ctest_array.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_array.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_array

makeout/unit/ctest_basis: unit/ctest_basis.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_basis.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_basis

makeout/unit/ctest_block_topo: unit/ctest_block_topo.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_block_topo.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_block_topo

makeout/unit/ctest_fv_proj: unit/ctest_fv_proj.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_fv_proj.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_fv_proj

makeout/unit/ctest_gauss_quad: unit/ctest_gauss_quad.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_gauss_quad.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_gauss_quad

makeout/unit/ctest_proj_on_basis: unit/ctest_proj_on_basis.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_proj_on_basis.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_proj_on_basis

makeout/unit/ctest_range: unit/ctest_range.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_range.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_range

makeout/unit/ctest_rect_apply_bc: unit/ctest_rect_apply_bc.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_rect_apply_bc.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_rect_apply_bc

makeout/unit/ctest_rect_decomp: unit/ctest_rect_decomp.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_rect_decomp.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_rect_decomp

makeout/unit/ctest_rect_grid: unit/ctest_rect_grid.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_rect_grid.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_rect_grid

makeout/unit/ctest_ref_count: unit/ctest_ref_count.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_ref_count.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_ref_count

makeout/unit/ctest_update_fsm: unit/ctest_update_fsm.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_update_fsm.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_update_fsm

makeout/unit/ctest_wv_euler: unit/ctest_wv_euler.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_wv_euler.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_wv_euler

makeout/unit/ctest_wv_iso_euler: unit/ctest_wv_iso_euler.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_wv_iso_euler.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_wv_iso_euler

makeout/unit/ctest_wv_maxwell: unit/ctest_wv_maxwell.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_wv_maxwell.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_wv_maxwell

makeout/unit/ctest_wv_ten_moment: unit/ctest_wv_ten_moment.c makeout/libgkylzero.a
	${CC} ${CFLAGS} unit/ctest_wv_ten_moment.c -I. -Lmakeout -lgkylzero -lm -lpthread -o makeout/unit/ctest_wv_ten_moment

# Run unit tests
check: $(patsubst %.c,makeout/%,$(wildcard unit/ctest_*.c))
	./makeout/unit/ctest_alloc
	./makeout/unit/ctest_array
	./makeout/unit/ctest_basis
	./makeout/unit/ctest_block_topo
	./makeout/unit/ctest_fv_proj
	./makeout/unit/ctest_gauss_quad
	./makeout/unit/ctest_proj_on_basis
	./makeout/unit/ctest_range
	./makeout/unit/ctest_rect_apply_bc
	./makeout/unit/ctest_rect_decomp
	./makeout/unit/ctest_rect_grid
	./makeout/unit/ctest_ref_count
	./makeout/unit/ctest_update_fsm
	./makeout/unit/ctest_wv_euler
	./makeout/unit/ctest_wv_iso_euler
	./makeout/unit/ctest_wv_maxwell
	./makeout/unit/ctest_wv_ten_moment

clean:
	rm -rf makeout/libgkylzero.a ${libobjs} makeout/regression/app_* makeout/unit/ctest_*

