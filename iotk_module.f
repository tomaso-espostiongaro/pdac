! INPUT/OUTPUT TOOLKIT
!
! This kit is intended to provide an easy access
! to binary and testual files formatted using
! some specific rule.
! The module defined below should be included
! in the calling routine.
! All public names exported frm this module
! has the "iotk_" prefix.

! This module contributed by Giovanni Bussi (www.unimore.it)

module iotk_module
implicit none

! All names are private ...
private

! ... except the names listed below
 public :: iotk_write_begin, &
           iotk_write_end,   &
           iotk_write_empty, &
           iotk_write_dat,   &
           iotk_write_attr,  &
           iotk_scan_begin,  &
           iotk_scan_end,    &
           iotk_scan_empty,  &
           iotk_scan_dat,    &
           iotk_scan_attr,   &
           iotk_taglenx,     &
           iotk_attlenx,     &
           iotk_namlenx,     &
           iotk_i2c

! Max number of controls
integer, parameter :: iotk_ncontrol = 255 ! (2**8)

! Max length of info tag (more exactly, the max length is iotk_taglenx-1)
integer, parameter :: iotk_taglenx =  65536 ! (2**16)
integer, parameter :: iotk_namlenx =  1000
integer, parameter :: iotk_attlenx =  iotk_taglenx - iotk_namlenx


! control = 1 <    >
! control = 2 </   >
! control = 3 <   />

integer, save :: iostat

interface iotk_scan_attr
module procedure iotk_scan_attr_char
module procedure iotk_scan_attr_int2
module procedure iotk_scan_attr_int4
module procedure iotk_scan_attr_real4
module procedure iotk_scan_attr_real8
module procedure iotk_scan_attr_complex4
module procedure iotk_scan_attr_complex8
module procedure iotk_scan_attr_logical
end interface iotk_scan_attr

interface iotk_write_attr
module procedure iotk_write_attr_char
module procedure iotk_write_attr_int2
module procedure iotk_write_attr_int4
module procedure iotk_write_attr_real4
module procedure iotk_write_attr_real8
module procedure iotk_write_attr_complex4
module procedure iotk_write_attr_complex8
module procedure iotk_write_attr_logical
end interface iotk_write_attr

interface iotk_write_dat
module procedure iotk_write_dat_logical_1
module procedure iotk_write_dat_integer2_1
module procedure iotk_write_dat_integer4_1
module procedure iotk_write_dat_complex4_1
module procedure iotk_write_dat_complex8_1
module procedure iotk_write_dat_real4_1
module procedure iotk_write_dat_real8_1
module procedure iotk_write_dat_real8
module procedure iotk_write_dat_integer
module procedure iotk_write_dat_logical
end interface iotk_write_dat

interface iotk_scan_dat
module procedure iotk_scan_dat_logical_1
module procedure iotk_scan_dat_integer2_1
module procedure iotk_scan_dat_integer4_1
module procedure iotk_scan_dat_complex4_1
module procedure iotk_scan_dat_complex8_1
module procedure iotk_scan_dat_real4_1
module procedure iotk_scan_dat_real8_1
end interface iotk_scan_dat

integer, parameter :: r4 = selected_real_kind(6)
integer, parameter :: r8 = selected_real_kind(13)
integer, parameter :: i2 = selected_int_kind(3)
integer, parameter :: i4 = selected_int_kind(6)


! Errors
! 101 :: iotk_unformatted | iostat/=0
! 102 :: write_tag        | iostat/=0
! 103-120 :: scan_tag     | iostat/=0
! 501 :: iotk_scan

contains

function iotk_i2c(i)
  integer, intent(in) :: i
  character(len=15)   :: iotk_i2c
  write(iotk_i2c,"(i15)") i
  iotk_i2c = adjustl(iotk_i2c)
end function iotk_i2c

function killpar(str)
  character(*), intent(in) :: str
  character(len(str))      :: killpar
  integer :: i
  killpar = " "
  do i = 1,len(str)
    if(str(i:i) /= "(" .and. str(i:i) /= ")") then
      killpar(i:i)=str(i:i)
    else
      killpar(i:i) = " "
    end if
  end do
end function killpar

subroutine iotk_error(unit,ierr)
  integer, intent(in) :: unit
  integer, intent(in) :: ierr
  character(1000) :: filename
  logical :: named
  inquire(unit=unit,named=named,name=filename)
  if(.not.named) filename="unnamed file"
  if(ierr==0) return
  write(0,*) "#################################"
  write(0,*) "Error in Input/Output Tool Kit"
  write(0,*) "Code        : ",ierr
  write(0,*) "Unit        : ",unit
  write(0,*) "File        : "//trim(filename)
  write(0,*) "Last iostat : ",iostat
  write(0,*) "#################################"
  stop
end subroutine iotk_error

function iotk_unformatted(unit,ierr)
  logical              :: iotk_unformatted
  integer, intent(in)  :: unit
  integer, intent(out) :: ierr
  character(30) :: form
  ierr = 0
  inquire(unit=unit,form=form,iostat=iostat)
  if(iostat/=0) then
    ierr = 101
    return
  end if
  if(form=="unformatted" .or. form=="UNFORMATTED") then
    iotk_unformatted = .true.
  else
    iotk_unformatted = .false.
  end if
end function iotk_unformatted

subroutine write_tag(unit,control,tag,unformatted,ierr)
  integer,                   intent(in)  :: unit
  integer,                   intent(in)  :: control
  character(iotk_taglenx-1), intent(in)  :: tag
  logical,                   intent(in)  :: unformatted
  integer,                   intent(out) :: ierr

  integer :: header,taglen
  character(2) :: begin,end

  ierr = 0
  taglen = len_trim(tag)
  if(unformatted) then
    header = control + taglen*iotk_ncontrol
    write(unit,iostat=iostat) header,tag(1:taglen)
  else
    select case(control)
    case(1)
      begin = "<"
      end   = ">"
    case(2)
      begin = "</"
      end   = ">"
    case(3)
      begin = "<"
      end   = "/>"
    end select
    write(unit,"(a)",iostat=iostat) &
            trim(begin)//tag(1:taglen)//trim(end)
  end if
  if(iostat/=0) ierr = 102
end subroutine write_tag

subroutine scan_tag(unit,direction,control,tag,unformatted,lall,ierr)
  integer,                   intent(in)  :: unit
  integer,                   intent(in)  :: direction
  integer,                   intent(out) :: control
  character(iotk_taglenx-1), intent(out) :: tag
  logical,                   intent(in)  :: lall
  logical,                   intent(in)  :: unformatted
  integer,                   intent(out) :: ierr

  integer :: header,taglen,pos,pos1,res
  character(2) :: begin,end
  character(4096) :: line
  logical :: found
  ierr = 0
  tag = " "
  if(unformatted) then
    found = .false.
    do
      if(direction<0) then
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 103
          return
        end if
      end if
      read(unit,iostat=iostat) header
      if(iostat/=0) then
        ierr = 104
        return
      end if
      control = modulo(header,iotk_ncontrol)
      if(control/=0) then
        found = .true.
        taglen  = modulo(header/iotk_ncontrol,iotk_taglenx)
        if(lall) then
          backspace(unit,iostat=iostat)
          if(iostat/=0) then
            ierr = 105
            return
          end if
          read(unit,iostat=iostat) header,tag(1:taglen)
          if(iostat/=0) then
            ierr = 106
            return
          end if
        end if
      end if
      if(direction<0) then
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 107
          return
        end if
      end if
      if(found) exit
    end do
  else
    tag = " "
    if(direction>=0) then
      do
        read(unit,"(a)",iostat=iostat) line
        if(iostat/=0) then
          ierr = 108
          return
        end if
        pos = scan(line,"<")
        if(pos/=0) exit
      end do
      do
        pos1 = scan(line(pos+1:),">") + pos
        if(pos1/=pos) exit
        tag(len_trim(tag)+2:) = trim(adjustl(line(pos+1:)))
        pos = 0
        read(unit,"(a)",iostat=iostat) line
        if(iostat/=0) then
          ierr = 109
          return
        end if
      end do
      tag(len_trim(tag)+2:) = trim(adjustl(line(pos+1:pos1-1)))
      res = len_trim(line(pos1+1:))
      if(res>0) then
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 110
          return
        end if
        read(unit,"(a)",iostat=iostat) line
        if(iostat/=0) then
          ierr = 111
          return
        end if
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 112
          return
        end if
        res = len_trim(line)-res
        read(unit,"(a)",iostat=iostat,advance='no') line(1:res)
        if(iostat/=0) then
          ierr = 113
          return
        end if
      end if
    else
      read(unit,"(a)",iostat=iostat) line
      if(iostat/=0) then
        ierr = 114
        return
      end if
      res = len_trim(line)
      do
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 115
          return
        end if
        read(unit,"(a)",iostat=iostat) line
        if(iostat/=0) then
          ierr = 116
          return
        end if
        pos = len_trim(line) - res
        pos = scan(line(1:pos),">",back=.true.)
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 117
          return
        end if
        if(pos/=0) exit
        res = 0
      end do
      do
        pos1 = scan(line(1:pos-1),"<",back=.true.)
        res = verify(tag," ")
        if(res==0) res=len(tag)
        if(pos1>0) exit
        tag(res-pos:res-2) = line(1:pos-1)
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 118
          return
        end if
        read(unit,"(a)",iostat=iostat) line
        if(iostat/=0) then
          ierr = 119
          return
        end if
        backspace(unit,iostat=iostat)
        if(iostat/=0) then
          ierr = 120
          return
        end if
        pos = len_trim(line)+1
      end do
      read(unit,"(a)",iostat=iostat,advance="no") line(1:pos1-1)
      tag(res-pos+pos1:res-2) = line(pos1+1:pos-1)
      tag = adjustl(tag)
    end if
    pos = verify(tag," ")
    pos1 = len_trim(tag)
    if(tag(pos:pos)=="/" .and. tag(pos1:pos1)/="/") then
      control = 2
      tag = tag(pos+1:pos1)
    else if(tag(pos:pos)/="/" .and. tag(pos1:pos1)=="/") then
      control = 3
      tag = tag(pos:pos1-1)
    else
      control = 1
      tag = tag(pos:pos1)
    end if
  end if
end subroutine scan_tag

subroutine iotk_scan(unit,direction,control,name,attr,unformatted,ierr)
  integer,                   intent(in)  :: unit
  integer,                   intent(in)  :: direction
  integer,                   intent(in)  :: control
  character(iotk_namlenx-1), intent(out) :: name
  character(iotk_attlenx-1), intent(out) :: attr
  logical,                   intent(in)  :: unformatted
  integer,                   intent(out) :: ierr

  character(iotk_taglenx-1) :: tag
  character(iotk_namlenx-1) :: r_name
  integer :: level,r_control,pos,pos1
  logical :: lall,match

  if(control==2 .and. direction<0) then
    ierr=501
    return
  end if
  level = 0
  ierr = 0
  do
    lall=.false. ! per ora leggo sempre
    if(direction>=0 .and. level==0) lall=.true.
    if(direction<0  .and. level==0 .and. control/=1) lall=.true.
    if(direction<0  .and. level==1 .and. control==1) lall=.true.
    call scan_tag(unit,direction,r_control,tag,unformatted,lall,ierr)
    if(ierr/=0) return
    pos = verify(tag," ")
    pos1 = scan(tag(pos:)," ") + pos - 2
    r_name = tag(pos:pos1)
    match = lall .and. r_control==control .and. r_name==name
    select case(direction)
    case(0:)
      select case(r_control)
      case(1)
        if(level==0 .and. match) exit
        level = level + 1
      case(2)
        if(level==0 .and. match) exit
        if(level==0) then
          call scan_tag(unit,-1,r_control,tag,unformatted,.false.,ierr)
          if(ierr==0) ierr=-1
          return
        end if
        level = level - 1
      case(3)
        if(level==0 .and. match) exit
      end select
    case(:-1)
      select case(r_control)
      case(2)
        level = level + 1
      case(1)
        if(level==1 .and. match) exit
        if(level==0) then
          call scan_tag(unit,+1,r_control,tag,unformatted,.false.,ierr)
          if(ierr==0) ierr=-1
          return
        end if
        level = level - 1
      case(3)
        if(level==0 .and. match) exit
      end select
    end select
  end do
  pos = verify(tag(pos1+1:)," ") + pos1
  if(pos/=pos1) then
    attr = trim(tag(pos:))
  else
    attr=" "
  end if
  if(direction<0) then
    call scan_tag(unit,+1,r_control,tag,unformatted,.false.,ierr)
    return
  end if
end subroutine iotk_scan

subroutine iotk_write_begin(unit,name,attr,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(in)  :: attr
  integer,      optional, intent(out) :: ierr
  character(iotk_taglenx-1) :: tag
  logical :: unformatted
  integer :: ierrl
  ierrl = 0
  if(present(attr)) then
    tag = trim(name)//" "//trim(attr)
  else
    tag = trim(name)
  end if
  unformatted = iotk_unformatted(unit,ierrl)
  if(ierrl==0) call write_tag(unit,1,tag,unformatted,ierrl)
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_begin
    
subroutine iotk_write_end(unit,name,attr,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(in)  :: attr
  integer,      optional, intent(out) :: ierr
  character(iotk_taglenx-1) :: tag
  logical :: unformatted
  integer :: ierrl
  ierrl = 0
  if(present(attr)) then
    tag = trim(name)//" "//trim(attr)
  else
    tag = trim(name)
  end if
  unformatted = iotk_unformatted(unit,ierrl)
  if(ierrl==0) call write_tag(unit,2,tag,unformatted,ierrl)
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_end
    
subroutine iotk_write_empty(unit,name,attr,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(in)  :: attr
  integer,      optional, intent(out) :: ierr
  character(iotk_taglenx-1) :: tag
  logical :: unformatted
  integer :: ierrl
  ierrl = 0
  if(present(attr)) then
    tag = trim(name)//" "//trim(attr)
  else
    tag = trim(name)
  end if
  unformatted = iotk_unformatted(unit,ierrl)
  if(ierrl==0) call write_tag(unit,3,tag,unformatted,ierrl)
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_empty

subroutine iotk_scan_begin(unit,name,attr,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(out) :: attr
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_namlenx-1) :: namel
  character(iotk_attlenx-1) :: attrl
  logical :: unformatted
  integer :: ierrl
  namel = name
  ierrl = 0
  unformatted = iotk_unformatted(unit,ierrl)
  if(ierrl==0) call iotk_scan(unit,1,1,namel,attrl,unformatted,ierrl)
  if(ierrl<0) write(0,*) "GIRO"
  if(ierrl<0) call iotk_scan(unit,-1,1,namel,attrl,unformatted,ierrl)
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(ierrl==0 .and. present(attr)) attr=attrl
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_begin

subroutine iotk_scan_end(unit,name,attr,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(out) :: attr
  integer,      optional, intent(out) :: ierr
  character(iotk_namlenx-1) :: namel
  character(iotk_attlenx-1) :: attrl
  logical :: unformatted
  integer :: ierrl
  namel = name
  ierrl = 0
  unformatted = iotk_unformatted(unit,ierrl)
  if(ierrl==0) call iotk_scan(unit,1,2,namel,attrl,unformatted,ierrl)
  if(ierrl==0 .and. present(attr)) attr=attrl
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_end

subroutine iotk_scan_empty(unit,name,attr,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  character(*), optional, intent(out) :: attr
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_namlenx-1) :: namel
  character(iotk_attlenx-1) :: attrl
  logical :: unformatted
  integer :: ierrl
  namel = name
  ierrl = 0
  unformatted = iotk_unformatted(unit,ierrl)
  if(ierrl==0) call iotk_scan(unit,1,3,namel,attrl,unformatted,ierrl)
  if(ierrl<0) call iotk_scan(unit,-1,3,namel,attrl,unformatted,ierrl)
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(ierrl==0 .and. present(attr)) attr=attrl
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_empty

function type_name(type)
  character(*), intent(in) :: type
  character(len(type))     :: type_name
  integer :: pos
  pos = verify(type,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_")
  if(pos==0) pos = len(type)+1
  type_name = type(1:pos-1)
end function type_name

function type_byte(type)
  character(*), intent(in) :: type
  character(len(type))     :: type_byte
  integer :: pos
  pos = verify(type,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_")
  if(pos==0) then
    type_byte =" "
  else
    type_byte = type(pos:pos)
  end if
end function type_byte

subroutine iotk_write_dat_logical_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  logical,                intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","logical")
    call iotk_write_attr(attr,"size",size(dat))
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(24l3)",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_logical_1


subroutine iotk_scan_dat_logical_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  logical,                intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=size(dat) .or. type_name(rtype)/="logical") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case(" ")
          read(unit,iostat=iostat) idummy,dat
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_logical_1

subroutine iotk_write_dat_integer4_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  integer(i4),            intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","integer4")
    call iotk_write_attr(attr,"size",size(dat))
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(4(i16))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_integer4_1


subroutine iotk_scan_dat_integer4_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  integer(i4),            intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  integer(i2), allocatable :: dat_i2(:)
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=size(dat) .or. type_name(rtype)/="integer") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case("4")
          read(unit,iostat=iostat) idummy,dat
          if(idummy/=0) ierrl = 1
        case("2")
          allocate(dat_i2(ubound(dat,1)))
          read(unit,iostat=iostat) idummy,dat_i2
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
          dat = dat_i2
          deallocate(dat_i2)
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_integer4_1

subroutine iotk_write_dat_integer2_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  integer(i2),            intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","integer2")
    call iotk_write_attr(attr,"size",size(dat))
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(6(i9))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_integer2_1

subroutine iotk_scan_dat_integer2_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  integer(i2),            intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  integer(i4), allocatable :: dat_i4(:)
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=size(dat) .or. type_name(rtype)/="integer") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case("2")
          read(unit,iostat=iostat) idummy,dat
          if(idummy/=0) ierrl = 1
        case("4")
          allocate(dat_i4(ubound(dat,1)))
          read(unit,iostat=iostat) idummy,dat_i4
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
          dat = dat_i4
          deallocate(dat_i4)
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_integer2_1

subroutine iotk_write_dat_complex4_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  complex(r4),               intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","complex4")
    call iotk_write_attr(attr,"size",2*size(dat))
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(4(D16.7))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_complex4_1

subroutine iotk_scan_dat_complex4_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  complex(r4),               intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  complex(r8), allocatable :: dat_r8(:)
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=2*size(dat) .or. type_name(rtype)/="complex") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case("4")
          read(unit,iostat=iostat) idummy,dat
          if(idummy/=0) ierrl = 1
        case("8")
          allocate(dat_r8(ubound(dat,1)))
          read(unit,iostat=iostat) idummy,dat_r8
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
          dat = dat_r8
          deallocate(dat_r8)
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_complex4_1

subroutine iotk_write_dat_complex8_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  complex(r8),               intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","complex8")
    call iotk_write_attr(attr,"size",2*size(dat))
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(4(D16.7))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_complex8_1

subroutine iotk_scan_dat_complex8_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  complex(r8),               intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  complex(r4), allocatable :: dat_r4(:)
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=size(dat) .or. type_name(rtype)/="complex") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case("8")
          read(unit,iostat=iostat) idummy,dat
          if(idummy/=0) ierrl = 1
        case("4")
          allocate(dat_r4(ubound(dat,1)))
          read(unit,iostat=iostat) idummy,dat_r4
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
          dat = dat_r4
          deallocate(dat_r4)
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_complex8_1

subroutine iotk_write_dat_real4_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  real(r4),               intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","real4")
    call iotk_write_attr(attr,"size",size(dat))
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(4(D16.7))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_real4_1

subroutine iotk_scan_dat_real4_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  real(r4),               intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  real(r8), allocatable :: dat_r8(:)
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=size(dat) .or. type_name(rtype)/="real") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case("4")
          read(unit,iostat=iostat) idummy,dat
          if(idummy/=0) ierrl = 1
        case("8")
          allocate(dat_r8(ubound(dat,1)))
          read(unit,iostat=iostat) idummy,dat_r8
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
          dat = dat_r8
          deallocate(dat_r8)
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_real4_1

subroutine iotk_write_dat_real8_1(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  real(r8),               intent(in)  :: dat(:)
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl, dsiz
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    dsiz = size(dat)
    call iotk_write_attr(attr,"type","real8")
    call iotk_write_attr(attr,"size",dsiz)
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(2(D23.14))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 1
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_real8_1


subroutine iotk_scan_dat_real8_1(unit,name,dat,found,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  real(r8),               intent(out) :: dat(:)
  logical,      optional, intent(out) :: found
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl,rsize,idummy,ierrltmp
  real(r4), allocatable :: dat_r4(:)
  character(50) :: rtype
  logical :: unformatted
  attr = " "
  ierrl = 0
  do 
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_scan_begin(unit,name,attr,ierr=ierrl)
    if(ierrl/=0) exit
    do
      call iotk_scan_attr(attr,"type",rtype,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      call iotk_scan_attr(attr,"size",rsize,ierr=ierrl)
      if(ierrl<0) ierrl=1
      if(ierrl/=0) exit
      if(rsize/=size(dat) .or. type_name(rtype)/="real") ierrl = 1
      if(ierrl/=0) exit
      if(unformatted) then
        select case(type_byte(rtype))
        case("8")
          read(unit,iostat=iostat) idummy,dat
          if(idummy/=0) ierrl = 1
        case("4")
          allocate(dat_r4(ubound(dat,1)))
          read(unit,iostat=iostat) idummy,dat_r4
          if(iostat/=0) ierrl = 1
          if(idummy/=0) ierrl = 1
          dat = dat_r4
          deallocate(dat_r4)
        case default
          ierrl = 1
        end select
      else
        read(unit,fmt=*,iostat=iostat) dat
        if(iostat/=0) ierrl = 1
      end if
    exit
    end do
    ierrltmp=ierrl
    call iotk_scan_end(unit,name,ierr=ierrl)
    if(ierrltmp/=0) ierrl=ierrltmp
    exit
  end do
  if(present(found)) then
    found = .false.
    if(ierrl==0) found = .true.
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(unit,ierrl)
  end if
end subroutine iotk_scan_dat_real8_1

subroutine iotk_write_attr_char(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  character(*), intent(in)    :: val
  integer :: pos
  pos = len_trim(attr)+2
  if(pos==2) pos=1
  attr(pos:) = name//'="'//val//'"'
end subroutine iotk_write_attr_char

subroutine iotk_scan_attr_char(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  character(*),      intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: equal,pos,pos1,ierrl
  logical :: foundl
  equal = 0
  foundl = .false.
  ierrl = 0
  do
    pos = scan(attr(equal+1:),"=")
    if(pos<=0) exit
    equal = equal + pos
    pos  = verify(attr(1:equal-1)," ",back=.true.)
    pos1 = scan(attr(1:pos)," ",back=.true.)
    if(attr(pos1+1:pos)==name) then
      foundl = .true.
      exit
    end if
  end do
  if(.not.foundl) ierrl = -1
  if(foundl) then
    pos  = scan(attr(equal+1:),'"') + equal
    pos1 = scan(attr(pos+1:),'"') + pos
    if(pos==pos1) ierrl = 1
    val = attr(pos+1:pos1-1)
  end if
  if(present(found)) found = foundl
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_char

subroutine iotk_write_attr_int4(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  integer(i4),  intent(in)    :: val
  character(20) :: valc
  write(valc,"(i20)") val
  call iotk_write_attr(attr,name,trim(adjustl(valc)))
end subroutine iotk_write_attr_int4

subroutine iotk_scan_attr_int4(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  integer(i4),       intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(100) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  if(ierrl==0) then
    valc = adjustl(valc)
    read(valc,"(i20)",iostat=iostat) val
    if(iostat/=0) ierrl=1
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_int4

subroutine iotk_write_attr_int2(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  integer(i2),  intent(in)    :: val
  character(20) :: valc
  write(valc,"(i20)") val
  call iotk_write_attr(attr,name,trim(adjustl(valc)))
end subroutine iotk_write_attr_int2

subroutine iotk_scan_attr_int2(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  integer(i2),       intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(100) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  if(ierrl==0) then
    valc = adjustl(valc)
    read(valc,"(i20)",iostat=iostat) val
    if(iostat/=0) ierrl=1
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_int2

subroutine iotk_write_attr_real4(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  real(r4),     intent(in)    :: val
  character(23) :: valc
  write(valc,"(D16.7)") val
  call iotk_write_attr(attr,name,trim(adjustl(valc)))
end subroutine iotk_write_attr_real4
 
subroutine iotk_scan_attr_real4(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  real(r4),          intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(100) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  if(ierrl==0) then
    valc = adjustl(valc)
    read(valc,"(G100.100)",iostat=iostat) val
    if(iostat/=0) ierrl=1
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_real4

subroutine iotk_write_attr_real8(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  real(r8),     intent(in)    :: val
  character(23) :: valc
  write(valc,"(D23.14)") val
  call iotk_write_attr(attr,name,trim(adjustl(valc)))
end subroutine iotk_write_attr_real8
 
subroutine iotk_scan_attr_real8(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  real(r8),          intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(100) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  if(ierrl==0) then
    valc = adjustl(valc)
    read(valc,"(G100.100)",iostat=iostat) val
    if(iostat/=0) ierrl=1
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_real8

subroutine iotk_write_attr_complex4(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  complex(r4),  intent(in)    :: val
  character(23) :: valc1
  character(23) :: valc2
  write(valc1,"(D16.7)") real(val)
  write(valc2,"(D16.7)") aimag(val)
  valc1 = adjustl(valc1)
  valc2 = adjustl(valc2)
  call iotk_write_attr(attr,name,"("//trim(valc1)//","//trim(valc2)//")")
end subroutine iotk_write_attr_complex4
 
subroutine iotk_scan_attr_complex4(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  complex(r4),       intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(100) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  valc = adjustl(killpar(valc))
  if(ierrl==0) read(valc,"(2G100.100)") val
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_complex4

subroutine iotk_write_attr_complex8(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  complex(r8),  intent(in)    :: val
  character(23) :: valc1
  character(23) :: valc2
  write(valc1,"(D23.14)") real(val)
  write(valc2,"(D23.14)") aimag(val)
  valc1 = adjustl(valc1)
  valc2 = adjustl(valc2)
  call iotk_write_attr(attr,name,"("//trim(valc1)//","//trim(valc2)//")")
end subroutine iotk_write_attr_complex8
 
subroutine iotk_scan_attr_complex8(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  complex(r8),       intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(100) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  valc = adjustl(killpar(valc))
  if(ierrl==0) read(valc,"(2G100.100)") val
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_complex8

subroutine iotk_write_attr_logical(attr,name,val)
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  logical,      intent(in)    :: val
  character(20) :: valc
  write(valc,"(l20)") val
  call iotk_write_attr(attr,name,trim(adjustl(valc)))
end subroutine iotk_write_attr_logical

subroutine iotk_scan_attr_logical(attr,name,val,found,ierr)
  character(*),      intent(in)  :: attr
  character(*),      intent(in)  :: name
  logical,           intent(out) :: val
  logical, optional, intent(out) :: found
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  logical :: foundl
  character(20) :: valc
  ierrl = 0
  call iotk_scan_attr(attr,name,valc,foundl,ierrl)
  valc = adjustl(valc)
  if(ierrl==0) read(valc,"(l20)") val
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl>0 .or. .not.present(found)) call iotk_error(0,ierrl)
  end if
end subroutine iotk_scan_attr_logical

! added by carlo

subroutine iotk_write_dat_real8(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  real(r8),               intent(in)  :: dat
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","real8")
    call iotk_write_attr(attr,"size",1)
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(2(D23.14))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 1
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_real8


subroutine iotk_write_dat_integer(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  integer,            intent(in)  :: dat
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","integer")
    call iotk_write_attr(attr,"size",1)
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(4(i16))",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_integer


subroutine iotk_write_dat_logical(unit,name,dat,ierr)
  integer,                intent(in)  :: unit
  character(*),           intent(in)  :: name
  logical,                intent(in)  :: dat
  integer,      optional, intent(out) :: ierr
  character(iotk_attlenx-1) :: attr
  integer :: ierrl
  logical :: unformatted
  attr = " "
  ierrl = 0
  do
    unformatted = iotk_unformatted(unit,ierrl)
    if(ierrl/=0) exit
    call iotk_write_attr(attr,"type","logical")
    call iotk_write_attr(attr,"size",1)
    call iotk_write_begin(unit,name,attr,ierrl)
    if(ierrl/=0) exit
    if(unformatted) then
      write(unit,iostat=iostat) 0,dat
    else
      write(unit,fmt="(24l3)",iostat=iostat) dat
    end if
    if(iostat/=0) ierrl = 150
    if(ierrl/=0) exit
    call iotk_write_end(unit,name,ierr=ierrl)
    exit
  end do
  if(present(ierr)) then
    ierr = ierrl
  else
    call iotk_error(unit,ierrl)
  end if
end subroutine iotk_write_dat_logical


end module iotk_module




