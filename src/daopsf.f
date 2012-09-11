      PROGRAM DAOPSF

      IMPLICIT NONE

      INTEGER MAXPSF, MAXPAR, MAXEXP, MAXROW, NOPT
      PARAMETER (MAXPSF=207, MAXPAR=6, MAXEXP=10, MAXROW=1024,NOPT=5)

      INTEGER IER, UNIT, RDPSF
      INTEGER NPAR, IPSTYP, NPSF, NEXP, NFRAC
      INTEGER IARG, ARGC, I, J
      INTEGER PSFSIZE, IMSIZE, NAXES(2)
      REAL PSFMAG, BRIGHT, XPSF, YPSF, DVDX, DVDY, SCALE
      REAL PSFX, PSFY, MAG, DX, DY, DXFROMPSF, DYFROMPSF
      REAL PSFRAD, PSFRADSQ, START, DELTA, DYSQ, DRSQ
      REAL PSF(MAXPSF,MAXPSF,MAXEXP), PAR(MAXPAR)
      REAL PSFROW(MAXROW), USEPSF
      CHARACTER PSFFILE*256, FITSFILE*260, SWITCH*260
      CHARACTER ARGV*32, C
      CHARACTER USAGE(NOPT+3)*80
      CHARACTER ERRMSG*80, ERRTXT*30

      DATA USAGE/ 'Usage: daopsf [OPTIONS] <file>',
     .     'Write a FITS image of a DAOPHOT PSF <file>',
     .     'OPTIONS:',
     .     '  -o=<file>: output (default: <file>.psf.fits)',
     .     '  -s=<n>   : size of image (default: psf radius)',
     .     '  -x=<n>   : x position of PSF (default: from file)',
     .     '  -y=<n>   : y position of PSF (default: from file)',
     .     '  -m=<n>   : magnitude of PSF (default:  1)'/
C
C     Parse command line arguments
C
      IMSIZE=-1
      XPSF=-1.
      YPSF=-1.
      MAG=99.
      FITSFILE=''
      ARGC = IARGC()
      IER=0

      IF (ARGC .LT. 1) THEN
         WRITE (6,101) (USAGE(IARG), IARG=1,NOPT+3)
 101     FORMAT(A50)
         GO TO 1000
      END IF

      DO IARG=1,ARGC
         CALL GETARG(IARG, ARGV)
         IF (ARGV(1:1) .EQ. '-') THEN
            C=ARGV(2:2)
            IF (C .EQ. 'o') THEN
               FITSFILE = ARGV(4:)
            ELSE IF (C .EQ. 's') THEN
               READ(ARGV(4:),*,IOSTAT=IER) IMSIZE
               IF (IMSIZE .GT. MAXROW) THEN
                  WRITE (6,102) IMSIZE, MAXROW
 102              FORMAT(' image size too big: ',I6, 'max:',I5)
                  GO TO 1000
               END IF
            ELSE IF (C .EQ. 'x') THEN
               READ(ARGV(4:),*,IOSTAT=IER) XPSF
            ELSE IF (C .EQ. 'y') THEN
               READ(ARGV(4:),*,IOSTAT=IER) YPSF
            ELSE IF (C .EQ. 'm') THEN
               READ(ARGV(4:),*,IOSTAT=IER) MAG
            ELSE
               WRITE (*,*) ' bad argument: ', ARGV
               GO TO 1000
            END IF
         ELSE
            PSFFILE=ARGV
         END IF
      END DO
      IF (IER .NE. 0) THEN
         WRITE (*,*) ' error parsing arguments '
         GO TO 1000
      END IF
C
C     Read PSF file
C
      IER = RDPSF(PSFFILE, IPSTYP, PAR, MAXPAR, NPAR,
     .            PSF, MAXPSF, MAXEXP, NPSF, NEXP, NFRAC, 
     .            PSFMAG, BRIGHT, PSFX, PSFY)

      IF (IER .NE. 0) THEN
         WRITE (*,*) ' error reading psf file: ', PSFFILE
         GO TO 1000
      END IF


      PSFRAD = (REAL(NPSF - 1) / 2. - 1.) / 2.
      PSFRADSQ = PSFRAD ** 2
      PSFSIZE  = 2 * INT(PSFRAD) + 1 
C
C     Initialize parameters
C
      IF (IMSIZE .LT. 0) THEN
         IMSIZE = PSFSIZE
      END IF
      NAXES(1) = IMSIZE
      NAXES(2) = IMSIZE

      IF (XPSF .LT. 0) THEN
         DX = PSFX
         DXFROMPSF = (PSFX - 1.0) / PSFX - 1.0
      ELSE
         DX = XPSF
         DXFROMPSF = (XPSF - 1.0) / PSFX - 1.0
      END IF

      IF (YPSF .LT. 0) THEN
         DY = PSFY
         DYFROMPSF = (PSFY - 1.0) / PSFY - 1.0
      ELSE 
         DY = YPSF
         DYFROMPSF = (YPSF - 1.0) / PSFY - 1.0
      END IF

      IF (MAG .GT. 90) THEN
         SCALE = 1.0
         MAG = PSFMAG
      ELSE
         SCALE = 10.**(0.4*(PSFMAG - MAG))
      END IF

      IF (FITSFILE .EQ. '') THEN
         FITSFILE = SWITCH(PSFFILE,'.psf.fits')
      END IF

      WRITE (6,103) DX, DY, MAG, IMSIZE
 103  FORMAT(' PSF evaluated at (', 
     .     F7.2, ', ', F7.2, ') Mag: ', F7.3, ' Size:', I5)

C
C     Initialize FITS array
C      
      CALL FTGIOU(UNIT, IER)
      CALL FTINIT(UNIT, FITSFILE, 1, IER)
      CALL FTPHPS(UNIT, -32, 2, NAXES, IER)
C
C     Loop over the PSF array to write each pixel
C
      START = - (PSFSIZE - 1) / 2.
      DELTA = REAL(PSFSIZE - 1) / REAL(IMSIZE - 1)

      DO J=1,IMSIZE
         DY= START + DELTA * (J - 1)
         DYSQ=DY**2
         DO I=1,IMSIZE
            DX = START + DELTA * (I - 1)
            DRSQ=DX**2 + DYSQ
            IF (DRSQ >= PSFRADSQ) THEN
               PSFROW(I) = 0.0
            ELSE
               PSFROW(I) = SCALE * USEPSF(IPSTYP, DX, DY, BRIGHT, 
     .              PAR, PSF, NPSF, NPAR, NEXP, NFRAC,
     .              DXFROMPSF, DYFROMPSF, DVDX, DVDY)               
            END IF
         END DO
         CALL FTPPRE(UNIT, 1, (J-1)*IMSIZE+1, IMSIZE, PSFROW, IER)
      END DO
C
C     Close up FITS array and print possible errors
C
      CALL FTCLOS(UNIT, IER)
      CALL FTFIOU(UNIT, IER)

      IF (IER .GT. 0) THEN
         CALL FTGERR(IER, ERRTXT)
         WRITE (*,*) 'FITSIO Error Status =', IER, ': ', ERRTXT
         CALL  FTGMSG(ERRMSG)
         DO WHILE (ERRMSG .NE. ' ')
            WRITE(*,*) ERRMSG
            CALL FTGMSG(ERRMSG)
         END DO
      END IF
      

 1000 CONTINUE
      END PROGRAM
