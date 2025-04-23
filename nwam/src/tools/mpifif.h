   !  include 'mpif-config.h'
   !   include 'mpif-constants.h'
   !   include 'mpif-handles.h'
   !   include 'mpif-io-constants.h'
   !   include 'mpif-io-handles.h'
   !   include 'mpif-externals.h'
   !   include 'mpif-sentinels.h'
   !   include 'mpif-sizeof.h'
   include 'mpif.h'
  interface
  subroutine mpi_bcast(buf,ni,datatype, root, comm,ierr)
  type(*), dimension(..),intent(inout) ::buf
  integer,intent(in):: ni,datatype,root,comm,ierr
  end subroutine 
  subroutine mpi_unpack(inbuf, ni, pos, outbuf, no, datatype, comm,ierr)
    type(*), dimension(..),intent(in) ::inbuf
    type(*), dimension(..) ::outbuf
    integer,intent(in):: ni,no,pos,datatype,comm,ierr
  end subroutine 
  subroutine mpi_pack(inbuf, ni,datatype, outbuf, no,pos,  comm,ierr)
    type(*), dimension(..),intent(in) ::inbuf
    type(*), dimension(..) ::outbuf
    integer,intent(in):: ni,no,pos,datatype,comm,ierr
  end subroutine 
  subroutine mpi_scatterv(sendbuf, sendcounts, displs,sendtype, recvbuf, recvcount,recvtype, root, comm,ierr)
    type(*), dimension(..),intent(in) ::sendbuf
    type(*), dimension(..) ::recvbuf
    integer,intent(in)::sendcounts(*), displs(*)
    integer,intent(in)::sendtype, recvcount,recvtype, root, comm,ierr
  end subroutine 
  subroutine mpi_gatherv(sendbuf, sendcount,sendtype, recvbuf, recvcounts, displs,recvtype, root, comm,ierr)
    type(*), dimension(..),intent(in) ::sendbuf
    type(*), dimension(..) ::recvbuf
    integer,intent(in)::recvcounts(*), displs(*)
    integer,intent(in)::sendtype, sendcount,recvtype, root, comm,ierr
  end subroutine 
  subroutine mpi_scatter(sendbuf, sendcount,sendtype, recvbuf, recvcount, recvtype, root, comm,ierr)
    type(*), dimension(..),intent(in) ::sendbuf
    type(*), dimension(..) ::recvbuf
    integer,intent(in)::sendtype, sendcount,recvtype,recvcount,root, comm,ierr
  end subroutine 
  subroutine mpi_gather(sendbuf, sendcount,sendtype, recvbuf, recvcount, recvtype, root, comm,ierr)
    type(*), dimension(..),intent(in) ::sendbuf
    type(*), dimension(..) ::recvbuf
    integer,intent(in)::sendtype, sendcount,recvtype,recvcount,root, comm,ierr
  end subroutine 
  subroutine mpi_allreduce(sendbuf, recvbuf, ni, dtype,op, comm,ierr)
    type(*), dimension(..),intent(in) ::sendbuf
    type(*), dimension(..) ::recvbuf
    integer,intent(in)::dtype, ni,op, comm,ierr
  end subroutine 
  subroutine mpi_reduce(sendbuf, recvbuf, ni, dtype,op,root, comm,ierr)
    type(*), dimension(..),intent(in) ::sendbuf
    type(*), dimension(..) ::recvbuf
    integer,intent(in)::dtype, ni,op, root,comm,ierr
  end subroutine 
  subroutine mpi_recv(buf, ni, dtype,src,tag, comm,stats,ierr)
    type(*), dimension(..)::buf
    integer,intent(in)::dtype, ni,src,tag,comm,ierr
    integer,intent(out)::stats(*)
  end subroutine 
  subroutine mpi_send(buf, ni, dtype,dest,tag, comm,ierr)
    type(*), dimension(..)::buf
    integer,intent(in)::dtype, ni,dest,tag,comm,ierr
  end subroutine 
  end interface
