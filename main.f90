program main
!---------------------------------------------------------------------
!   Modules
!---------------------------------------------------------------------

    use iso_fortran_env, only: input_unit
    use matrix_m       , only: matrix_t
    use real_kind      , only: DP   
    use types_m        , only: problem_t
    use sol_m          , only: solution_t, independent_t, act_Pe_q
    use vector_m       , only: omega_t, vector_t
    implicit none
 
!---------------------------------------------------------------------
!   Variables
!---------------------------------------------------------------------
     
    real(DP) :: t_start, t_finish
    integer, pointer :: n
    real(DP) :: error, q_p, Pe_p
    integer :: i, j,k, channel = 2 
    
    
    type(problem_t), target     :: problem
    type(matrix_t)      :: channels(2)
    type(omega_t)       :: omega(2)
    type(solution_t)    :: sol
    type(independent_t) :: b
    CALL CPU_TIME(t_start)

!---------------------------------------------------------------------
!   Initialization Of The Problem
!---------------------------------------------------------------------
     
    ! read input constants file
    call problem%initialize("datainput")

    n => problem%n 
    if(problem%fixt <= problem%xmax) then 
        call sol%initialize(problem%h, n, channel, &
        problem%xmin, problem%q, problem%fixt, problem%i_ft)
    else
        call sol%initialize(problem%h, n, channel, problem%xmin, &
             problem%q)
    endif
    
    if (problem%i_data == 1) then
        call sol%read(n)
    end if

    do j = 1, channel
        call omega(j)%initialize(n)
    end do
    
    call b%initialize(problem%h, n, channel, problem%xmin, problem%q)
    call b%act_indep(sol%T(:), sol%Z(:), problem, channel)

    do j = 1, channel
        call channels(j)%matrix_allocate(problem,j)
        call omega(j)%act_omega(sol%T(j), n, problem%beta, &
        problem%gamma)
        call channels(j)%matrix_act(problem,omega(j),      &
        sol%F(j), sol%Z(j))
    end do
        
    
    call sol%write_temp(problem, omega(:), 0.0_DP, 0)
        
    q_p = problem%q
    Pe_p = problem%Pe(1)

!---------------------------------------------------------------------
!   Solver
!---------------------------------------------------------------------
    do i = 1, problem%m_i
        
        call sol%act_sol(channels(:), b, problem, channel)
        call sol%cal_error(error, channel)
        
        if (problem%q_Pe /= 0) then
            error = max(error, abs(problem%q-q_p),    &
            abs(problem%Pe(1)-Pe_p))
        endif
        
        if (error < problem%tol .or. i > problem%m_i) then
            if (error == 0.0 .or. error /= error ) then
                print *, "SOLUTION NOT FOUND"
                exit
            end if
            call sol%write_temp(problem, omega(:), error, i)
            CALL CPU_TIME(t_finish)
            call sol%write_end(problem, omega, t_start, t_finish)
            exit
        endif
        
        call b%act_indep(sol%T(:),sol%Z(:), problem, channel)
        do j = 1,2
            call omega(j)%act_omega(sol%T(j), n, problem%beta, &
            problem%gamma)
            call channels(j)%matrix_act(problem, omega(j), &
            sol%F(j), sol%Z(j))
        end do
        
        if (mod(i,problem%show) == 0) then
            call sol%write_temp(problem, omega(:), error, i)
            
            Pe_p = problem%Pe(1)
            q_p = problem%q
        end if
     
    enddo


    
    do j = 1, channel
        call omega(j)%free_vector()
        call channels(j)%free_matrix()
    end do
    
    call b%free_indep()
    call sol%free_indep()
    call problem%free_problem()

    
end program main
