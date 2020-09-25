module interfaces_m
  
  use type_m   , only: DP
  use param_m  , only: param_t
  use paralel_m, only: parll_t
  
  implicit none
  
  abstract interface
  subroutine coef(points, param)
    import :: DP, param_t
    real(DP)     , intent(inout) :: points(:,:)
    type(param_t), intent(in   ) :: param
  end subroutine
  end interface

  
end module interfaces_m

