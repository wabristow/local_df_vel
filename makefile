CFLAGS= -D _GNU_SOURCE -D __USE_XOPEN  -I$(RSTPATH)/include -L /usr/local/lib -L/usr/intel -L$(RSTPATH)/lib -O3 -ipo -mcmodel medium -shared_intel -fPIC -Wall -D_GNU_SOURCE -D_LINUX -D_Float32=float -D_Float64=float -D_Float32x=float -D_Float64x=float

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn -I$(IPATH)/analysis

LIBS=-lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -laacgm_v2.1 -ligrf_v2.1 -lastalg.1 -lrpos.1  -lcnvmap.1 -loldcnvmap.1 -lshf.1 -lgrd.1 -loldgrd.1 -laacgm.1 -ldmap.1 -lrfile.1 -lopt.1 -lcfit.1 -ligrf.1

.c.o:
	icc $(CFLAGS) $(INCLUDE) -c -o $@ $<

OBJS=local_df_vel_medfilt.o sub_sphazm.o sub_sphcal.o sub_gcp.o find_fit_files.o make_local_pgrid.o pinv.o

local_df_vel:	$(OBJS)
	icc -o ~/bin/local_df_vel $(CFLAGS) $(INCLUDE) $(OBJS)  -lgsl -lgslcblas -lm  -mkl=parallel -qopenmp $(LIBS)

