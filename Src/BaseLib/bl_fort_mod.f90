module bl_fort_module

  use iso_c_binding, only : c_size_t
  implicit none

  integer, parameter :: bl_real = kind(0.d0)
  integer (kind=c_size_t), parameter :: bl_real_size = 8_c_size_t

end module bl_fort_module
