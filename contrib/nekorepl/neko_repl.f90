module neko_repl
  use neko
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
contains

  subroutine neko_repl_init() bind(c, name='repl_init')

  
    call neko_init()
    
  end subroutine neko_repl_init
  
end module neko_repl
