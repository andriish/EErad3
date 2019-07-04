      program eerad3_combine
      implicit real*8(a-h,o-z)
      character*20 fname,frooty,frooti
      character*2 ianame
      character*40 ibname
      character*40 ifile
      character*2 filetag
      dimension ivoid(100)
      common/voidlist/iv(100),nv
      common/outfile/filetag

      n = iargc()
      ifile = 'eerad3_combine.input'
      if (n.eq.2) then 
         call getarg(1,ianame)
         call getarg(2,ibname)
         if (ianame.eq.'-i') ifile = ibname
      endif 
 
      open(9,file=ifile)
      read(9,*) iaver
      read(9,*) frooty
      read(9,*) frooti
      read(9,*) filetag
      read(9,*) minfile,maxfile
      read(9,*) nv
      do i=1,nv
         read(9,*) iv(i)
      enddo
      close(9)

      numfiles = maxfile-minfile-nv
      if (numfiles.le.0) then
         write(6,*) 
     .   'Error in eerad3_combine: only single set of input files'
         stop
      endif
      
      fname = '.'//frooty(1:4)//'.'//frooti(1:3)//'.'

      if (iaver.eq.8.or.iaver.eq.6) then
      call combinehistoL(2,fname(1:10)//'Y3a',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'Y4a',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'Y5a',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'S3a',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'S4a',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'S5a',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'Y3b',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'Y4b',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'Y5b',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'S3b',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'S4b',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'S5b',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'Y3c',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'Y4c',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'Y5c',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'S3c',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'S4c',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'S5c',minfile,maxfile)
      endif
      if (iaver.eq.8) then
      call combinehistoL(2,fname(1:10)//'TLa',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'CLa',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'BLa',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'WLa',minfile,maxfile)
      call combinehistoL(2,fname(1:10)//'MLa',minfile,maxfile)
CAV*********************************************************************
      call combinehistoL(12,fname(1:10)//'ELa',minfile,maxfile)
CAV*********************************************************************
      call combinehistoL(3,fname(1:10)//'TLb',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'CLb',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'BLb',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'WLb',minfile,maxfile)
      call combinehistoL(3,fname(1:10)//'MLb',minfile,maxfile)
CAV*********************************************************************
      call combinehistoL(13,fname(1:10)//'ELb',minfile,maxfile)
CAV*********************************************************************
      call combinehistoL(4,fname(1:10)//'TLc',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'CLc',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'BLc',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'WLc',minfile,maxfile)
      call combinehistoL(4,fname(1:10)//'MLc',minfile,maxfile)
CAV*********************************************************************
      call combinehistoL(14,fname(1:10)//'ELc',minfile,maxfile)
CAV*********************************************************************
      endif
      if (iaver.eq.0.or.iaver.eq.1) then
      call combinehisto(1,fname(1:10)//'W1a',minfile,maxfile)
      call combinehisto(1,fname(1:10)//'W2a',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'W1b',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'W2b',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'W1c',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'W2c',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'W1d',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'W2d',minfile,maxfile)
      endif
      if (iaver.eq.0.or.iaver.eq.2) then
      call combinehisto(0,fname(1:10)//'C1a',minfile,maxfile)
      call combinehisto(0,fname(1:10)//'C2a',minfile,maxfile)
      call combinehisto(1,fname(1:10)//'C1b',minfile,maxfile)
      call combinehisto(1,fname(1:10)//'C2b',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'C1c',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'C2c',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'C1d',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'C2d',minfile,maxfile)
      endif
      if (iaver.eq.0.or.iaver.eq.3) then
      call combinehisto(1,fname(1:10)//'M1a',minfile,maxfile)
      call combinehisto(1,fname(1:10)//'M2a',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'M1b',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'M2b',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'M1c',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'M2c',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'M1d',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'M2d',minfile,maxfile)
      endif
      if (iaver.eq.0.or.iaver.eq.4) then
      call combinehisto(1,fname(1:10)//'T1a',minfile,maxfile)
      call combinehisto(1,fname(1:10)//'T2a',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'T1b',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'T2b',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'T1c',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'T2c',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'T1d',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'T2d',minfile,maxfile)
      endif
      if (iaver.eq.0.or.iaver.eq.5) then
      call combinehisto(1,fname(1:10)//'B1a',minfile,maxfile)
      call combinehisto(1,fname(1:10)//'B2a',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'B1b',minfile,maxfile)
      call combinehisto(2,fname(1:10)//'B2b',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'B1c',minfile,maxfile)
      call combinehisto(3,fname(1:10)//'B2c',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'B1d',minfile,maxfile)
      call combinehisto(4,fname(1:10)//'B2d',minfile,maxfile)
      endif




      stop
      end


      subroutine combinehistoL(inum,fname,minfile,maxfile)
      implicit real*8(a-h,o-z)
      character*22 fname
      parameter (nbins0=400,nbins1=200,nbins2=100,nbins3=50,nbins4=25)
      dimension t0(1:nbins0),t1(1:nbins1),t2(1:nbins2),
     .          t3(1:nbins3),t4(1:nbins4)
      dimension yi0(minfile:maxfile,1:nbins0)
     .         ,ei0(minfile:maxfile,1:nbins0)
     .         ,wi0(minfile:maxfile,1:nbins0)
     .         ,di0(minfile:maxfile,1:nbins0)
      dimension yi1(minfile:maxfile,1:nbins1)
     .         ,ei1(minfile:maxfile,1:nbins1)
     .         ,wi1(minfile:maxfile,1:nbins1)
     .         ,di1(minfile:maxfile,1:nbins1)
      dimension yi2(minfile:maxfile,1:nbins2)
     .         ,ei2(minfile:maxfile,1:nbins2)
     .         ,wi2(minfile:maxfile,1:nbins2)
     .         ,di2(minfile:maxfile,1:nbins2)
      dimension yi3(minfile:maxfile,1:nbins3)
     .         ,ei3(minfile:maxfile,1:nbins3)
     .         ,wi3(minfile:maxfile,1:nbins3)
     .         ,di3(minfile:maxfile,1:nbins3)
      dimension yi4(minfile:maxfile,1:nbins4)
     .         ,ei4(minfile:maxfile,1:nbins4)
     .         ,wi4(minfile:maxfile,1:nbins4)
     .         ,di4(minfile:maxfile,1:nbins4)
CAV*********************************************************************
      parameter (nbins0EC=720,nbins1EC=360,nbins2EC=180,nbins3EC=90,
     .nbins4EC=45)
      dimension t0EC(1:nbins0EC),t1EC(1:nbins1EC),t2EC(1:nbins2EC),
     .          t3EC(1:nbins3EC),t4EC(1:nbins4EC)
      dimension yi0EC(minfile:maxfile,1:nbins0EC)
     .         ,ei0EC(minfile:maxfile,1:nbins0EC)
     .         ,wi0EC(minfile:maxfile,1:nbins0EC)
     .         ,di0EC(minfile:maxfile,1:nbins0EC)
      dimension yi1EC(minfile:maxfile,1:nbins1EC)
     .         ,ei1EC(minfile:maxfile,1:nbins1EC)
     .         ,wi1EC(minfile:maxfile,1:nbins1EC)
     .         ,di1EC(minfile:maxfile,1:nbins1EC)
      dimension yi2EC(minfile:maxfile,1:nbins2EC)
     .         ,ei2EC(minfile:maxfile,1:nbins2EC)
     .         ,wi2EC(minfile:maxfile,1:nbins2EC)
     .         ,di2EC(minfile:maxfile,1:nbins2EC)
      dimension yi3EC(minfile:maxfile,1:nbins3EC)
     .         ,ei3EC(minfile:maxfile,1:nbins3EC)
     .         ,wi3EC(minfile:maxfile,1:nbins3EC)
     .         ,di3EC(minfile:maxfile,1:nbins3EC)
      dimension yi4EC(minfile:maxfile,1:nbins4EC)
     .         ,ei4EC(minfile:maxfile,1:nbins4EC)
     .         ,wi4EC(minfile:maxfile,1:nbins4EC)
     .         ,di4EC(minfile:maxfile,1:nbins4EC)
      if (inum.eq.10) then
      call loadhisto(nbins0EC,yi0EC,ei0EC,wi0EC,di0EC,t0EC,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins0EC,yi0EC,ei0EC,wi0EC,di0EC,t0EC,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.11) then
      call loadhisto(nbins1EC,yi1EC,ei1EC,wi1EC,di1EC,t1EC,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins1EC,yi1EC,ei1EC,wi1EC,di1EC,t1EC,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.12) then
      call loadhisto(nbins2EC,yi2EC,ei2EC,wi2EC,di2EC,t2EC,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins2EC,yi2EC,ei2EC,wi2EC,di2EC,t2EC,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.13) then
      call loadhisto(nbins3EC,yi3EC,ei3EC,wi3EC,di3EC,t3EC,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins3EC,yi3EC,ei3EC,wi3EC,di3EC,t3EC,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.14) then
      call loadhisto(nbins4EC,yi4EC,ei4EC,wi4EC,di4EC,t4EC,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins4EC,yi4EC,ei4EC,wi4EC,di4EC,t4EC,fname(1:13)
     .               ,minfile,maxfile)
      endif

CAV*********************************************************************
      if (inum.eq.0) then
      call loadhisto(nbins0,yi0,ei0,wi0,di0,t0,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins0,yi0,ei0,wi0,di0,t0,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.1) then
      call loadhisto(nbins1,yi1,ei1,wi1,di1,t1,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins1,yi1,ei1,wi1,di1,t1,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.2) then
      call loadhisto(nbins2,yi2,ei2,wi2,di2,t2,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins2,yi2,ei2,wi2,di2,t2,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.3) then
      call loadhisto(nbins3,yi3,ei3,wi3,di3,t3,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins3,yi3,ei3,wi3,di3,t3,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.4) then
      call loadhisto(nbins4,yi4,ei4,wi4,di4,t4,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins4,yi4,ei4,wi4,di4,t4,fname(1:13)
     .               ,minfile,maxfile)
      endif
      return
      end

      subroutine combinehisto(inum,fname,minfile,maxfile)
      implicit real*8(a-h,o-z)
      character*22 fname
      parameter (nbins0=400,nbins1=200,nbins2=100,nbins3=50,nbins4=25)
      dimension t0(1:nbins0),t1(1:nbins1),t2(1:nbins2),
     .          t3(1:nbins3),t4(1:nbins4)
      dimension yi0(minfile:maxfile,1:nbins0)
     .         ,ei0(minfile:maxfile,1:nbins0)
     .         ,wi0(minfile:maxfile,1:nbins0)
     .         ,di0(minfile:maxfile,1:nbins0)
      dimension yi1(minfile:maxfile,1:nbins1)
     .         ,ei1(minfile:maxfile,1:nbins1)
     .         ,wi1(minfile:maxfile,1:nbins1)
     .         ,di1(minfile:maxfile,1:nbins1)
      dimension yi2(minfile:maxfile,1:nbins2)
     .         ,ei2(minfile:maxfile,1:nbins2)
     .         ,wi2(minfile:maxfile,1:nbins2)
     .         ,di2(minfile:maxfile,1:nbins2)
      dimension yi3(minfile:maxfile,1:nbins3)
     .         ,ei3(minfile:maxfile,1:nbins3)
     .         ,wi3(minfile:maxfile,1:nbins3)
     .         ,di3(minfile:maxfile,1:nbins3)
      dimension yi4(minfile:maxfile,1:nbins4)
     .         ,ei4(minfile:maxfile,1:nbins4)
     .         ,wi4(minfile:maxfile,1:nbins4)
     .         ,di4(minfile:maxfile,1:nbins4)


      if (inum.eq.0) then
      call loadhisto(nbins0,yi0,ei0,wi0,di0,t0,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins0,yi0,ei0,wi0,di0,t0,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.1) then
      call loadhisto(nbins1,yi1,ei1,wi1,di1,t1,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins1,yi1,ei1,wi1,di1,t1,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.2) then
      call loadhisto(nbins2,yi2,ei2,wi2,di2,t2,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins2,yi2,ei2,wi2,di2,t2,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.3) then
      call loadhisto(nbins3,yi3,ei3,wi3,di3,t3,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins3,yi3,ei3,wi3,di3,t3,fname(1:13)
     .               ,minfile,maxfile)
      elseif (inum.eq.4) then
      call loadhisto(nbins4,yi4,ei4,wi4,di4,t4,fname(1:13)
     .               ,minfile,maxfile)
      call writehisto(nbins4,yi4,ei4,wi4,di4,t4,fname(1:13)
     .               ,minfile,maxfile)
      endif
      return
      end



      subroutine loadhisto(nbins,yi,ei,wi,di,t,froot,minfile,maxfile)
      implicit real*8(a-h,o-z)
      character*22 froot,fname
      logical fi
      dimension yi(minfile:maxfile,1:nbins)
     .         ,ei(minfile:maxfile,1:nbins)
     .         ,wi(minfile:maxfile,1:nbins)
     .         ,di(minfile:maxfile,1:nbins)
     .         ,t(1:nbins)
      common/voidlist/iv(100),nv

      do i=minfile,maxfile    
         iflag = 0
         do j=1,nv
            if (i.eq.iv(j)) iflag = 1
         enddo
         if (iflag.eq.1) goto 60   
         write(fname,'(I4)') i
         if (i.lt.10) fname='E000'//fname(4:4)//froot(1:13)
         if (i.ge.10.and.i.lt.100) 
     .             fname='E00'//fname(3:4)//froot(1:13)
         if (i.ge.100.and.i.lt.1000) 
     .             fname='E0'//fname(2:4)//froot(1:13)
         if (i.ge.1000.and.i.lt.10000) 
     .             fname='E'//fname(1:4)//froot(1:13)
         inquire(file=fname,EXIST=fi)
         if (.not.fi) then 
            write(6,*) 'Missing file: ', fname
            stop
         endif
         open(11,file=fname)
         do j=1,nbins
            read(11,*) t(j),yi(i,j),ei(i,j)
            wi(i,j) = 0d0
            if (ei(i,j).ne.0d0) wi(i,j)=1d0/ei(i,j)**2
c            wi(i,j) = 1d0
         enddo
         close(11)
 60      continue
      enddo
      return
      end

      subroutine writehisto(nbins,yi,ei,wi,di,t,froot,minfile,maxfile)
      implicit real*8(a-h,o-z)
      character*20 froot,fname
      dimension yi(minfile:maxfile,1:nbins)
     .         ,ei(minfile:maxfile,1:nbins)
     .         ,wi(minfile:maxfile,1:nbins)
     .         ,di(minfile:maxfile,1:nbins)
      dimension t(1:nbins),ymean(1:nbins)
     .          ,wsum(1:nbins),errymean(1:nbins)
      character*2 filetag
      common/voidlist/iv(100),nv
      common/outfile/filetag
      
      sum = 0d0
      errsum = 0d0
      open(11,file='E'//filetag//froot(1:13))
      do j=1,nbins
         ymean(j) = 0d0
         wsum(j) = 0d0
         do i=minfile,maxfile
            iflag = 0
            do k=1,nv
               if (i.eq.iv(k)) iflag = 1
            enddo
            if (iflag.eq.1) goto 70   
            ymean(j) = ymean(j)+wi(i,j)*yi(i,j)
            wsum(j) = wsum(j)+wi(i,j)
 70         continue
         enddo
         if (wsum(j).ne.0d0) ymean(j) = ymean(j)/wsum(j)
         sum = sum+ymean(j)

         errymean(j) = 0d0
         do i=minfile,maxfile
            iflag = 0
            do k=1,nv
               if (i.eq.iv(k)) iflag = 1
            enddo
            if (iflag.eq.1) goto 80   
            di(i,j) = yi(i,j)-ymean(j)
            errymean(j) = errymean(j)+wi(i,j)*di(i,j)**2
 80         continue
         enddo
         if (wsum(j).ne.0d0) errymean(j) = errymean(j)/
     .               real(maxfile-minfile-nv)/wsum(j)
         errsum = errsum + errymean(j)

         write(11,101) t(j),ymean(j),dsqrt(errymean(j))
      enddo
      close(11)
c      write(6,*) sum/16d0,dsqrt(errsum)/16d0
 101  format(3x,f11.6,1pe12.4,1pe12.4)
      return
      end
