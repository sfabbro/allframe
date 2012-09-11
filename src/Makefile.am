MAINTAINERCLEANFILES = Makefile.in
CLEANFILES = *~

bin_PROGRAMS = daophot allstar daomaster daomatch montage2 allframe \
	daogrow daopsf

noinst_HEADERS = arrays.inc

if WITH_FITS
if WITH_IRAF
io_source = bothsubs.f
else !WITH_IRAF
io_source = fitsubs.f
endif
endif

if WITH_IRAF
if !WITH_FITS
io_source = irasubs.f
endif
endif

daophot_SOURCES = \
	daophot.f pckpsf.f find.f fotometry.f \
	psf.f peak.f nstar.f fudge.f addstar.f substar.f \
	group.f sort.f lnxsubs.f iosubs.f mathsubs.f $(io_source)
daophot_LDADD = @CFITSIO_LIBS@

daopsf_SOURCES = daopsf.f mathsubs.f lnxsubs.f iosubs.f $(io_source)
daopsf_LDADD = @CFITSIO_LIBS@ @IRAF_LIBS@

allstar_SOURCES = allstar.f allstsubs.f lnxsubs.f iosubs.f \
	mathsubs.f $(io_source)
allstar_LDADD = @CFITSIO_LIBS@ @IRAF_LIBS@

allframe_SOURCES = allframe.f lnxsubs.f iosubs.f mathsubs.f $(io_source)
allframe_LDADD = @CFITSIO_LIBS@ @IRAF_LIBS@

montage2_SOURCES = montage2.f mathsubs.f iosubs.f lnxsubs.f $(io_source)
montage2_LDADD = @CFITSIO_LIBS@ @IRAF_LIBS@

daomaster_SOURCES = daomaster.f iosubs.f mathsubs.f lnxsubs.f dummysm.f
daomatch_SOURCES = daomatch.f iosubs.f mathsubs.f lnxsubs.f
daogrow_SOURCES = daogrow.f iosubs.f mathsubs.f lnxsubs.f

doc_DATA = ccdpck.man cook3.1.ps daophotii.ps daophotii.tex README.txt

EXTRA_DIST = airmass.f bothsubs.f ccdave.f ccdstd.f irasubs.f newtrial.f \
	redsubs.f unxsubs.f $(doc_DATA)
