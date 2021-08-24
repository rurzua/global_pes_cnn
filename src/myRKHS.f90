module myRKHS
    implicit none
    private
    !declare all methods/variables available from the outside
    public :: q_ker, dq_ker, k_ker, dk_ker, ytrans, ytransinv
    
    !functions and subroutines start here
    contains
    !the reproducing kernel q_ker for distance like variables
    real*8 function q_ker(x,r)
        implicit none
        real*8 :: x,r,x_l,x_s
        x_l=r
           x_s=x
        if(r < x) then
            x_s=r
            x_l=x
        endif       
        !for n=2 and m=5
        !q_ker = (2d0/(21d0*x_l**6)) - (1d0*x_s/(14d0*x_l**7))

        !for n=2 and m=6
        q_ker = (1d0/(14d0*x_l**7))*(1d0-7d0*x_s/(9d0*x_l))

        !for n=2 and m=4
        !q_ker = (2d0/(15d0*x_l**5)) * (1d0-5*x_s/(7*x_l))

        !for n=2 and m=3
        !q_ker = (1d0/(5d0*x_l**4)) - (2.d0*x_s/(15.d0*x_l**5))		

    end function q_ker
    
    !derivative of q_ker
    real*8 function dq_ker(x,r)
        implicit none
        real*8 :: x,r
        if(x < r) then
            dq_ker = -1d0/(18d0*r**8)
        else
            dq_ker = -(1d0/18d0)*(9*x-8*r)/(x**9)
        end if
    end function dq_ker


    !the reproducing kernel k_ker for angle like variables
    real*8 function k_ker(x,r)      !k=2 Anglelike
        implicit none
        real*8 :: x,r,x_l,x_s
        x_l=r !Large value of the coordinate
        x_s=x !Small value of the coordinate
        if(r < x) then !Here the smaller and larger value are chosen
            x_s=r
            x_l=x
        endif
        !Only two terms (up to n=1) were taken of the Gauss´ hypergeometric function expansion  
        k_ker=1d0+x_s*x_l+2d0*x_s**2d0*x_l*(1.d0-x_s/(3d0*x_l))
    end function k_ker

    !derivative of k_ker
    real*8 function dk_ker(x,r)
        implicit none
        real*8 :: x,r
        if(x < r) then
            dk_ker = r+4*x*r*(1.d0-x/(3*r))-(2*x**2)/3d0
        else
            dk_ker = r+2*r**2*(1.d0-r/(3*x))+(2*r**3)/(3*x)
        end if
    end function dk_ker

    real*8 function ytrans(alpha)
        implicit none
        real*8,intent(in) :: alpha
        ytrans = (1d0-DCOSD(alpha))/2d0
    end function ytrans

    real*8 function ytransinv(y)
        implicit none
        real*8,intent(in) :: y
        ytransinv = DACOSD(1d0-2d0*y)
    end function ytransinv
end module myRKHS

