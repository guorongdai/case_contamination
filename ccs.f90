!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve the estimating equation by Newton's method
! par is the output of the parameter
! abetainv is the estimated A(\beta)^{-1} in the asymptotic covariance 
! r is the missingness indicator. When s = 0, r = 1 always holds.
! s is the indicator of whether an individual is in the case pool. 
! d is the indicator of whether a case is a real case.
! u is the offset of the unlabelled set. v is the offset of the labelled set. 
! x is the covariate matrix of the association model. 
! w is the covariate matrix of the phenotyping model.
! p1 is the dimensionality of x. p2 is the dimensionality of w.
! p = p1 + p2.
! fi = 1 means full imputation; f1 = 0 means partial imputaiton.

subroutine eesolver(par, abetainv, intval, r, s, d, u, v, x, w, n, p1, p2, p, fi, eps, maxout)
implicit none

integer*4 :: n, p1, p2, maxout, i, j, p, fi
real*8 :: par(p), intval(p), r(n), s(n), d(n), u(n), v(n), &
          x(n, p1), w(n, p2), gra(p, p), invgra(p, p), ee(p, 1), eps, diff, &
          theta(p1), alpha(p2), hf, bh, dhf, dbh, xi(p1, 1), wi(p2, 1), &
          z1, z2, z3, z4, ss1, ss2, dir(p, 1), abetainv(p, p)



diff = 10.d0 ** 5
j = 1


do while((j <= maxout) .AND. (diff >= eps))

  ee = 0.d0
  gra = 0.d0
  
  theta = intval(1 : p1)
  alpha = intval( (p1 + 1) : p )


  do i = 1, n
  
  
    ss1 = sum(theta * x(i, :))
    ss2 = sum(alpha * w(i, :))
    
    if (fi == 0) then  
    
      z1 = s(i) * ( r(i) * d(i) + (1 - r(i)) * hf( u(i) + ss2 ) )
      
      z3 = s(i) * (1 - r(i)) * dhf( u(i) + ss2 ) * bh(ss1)
      
    else
      
      z1 = s(i) * hf( u(i) + ss2 ) 
      
      z3 = s(i) * dhf( u(i) + ss2 ) * bh(ss1)
      
    end if
    
      
    z2 = s(i) * r(i) * ( d(i) - hf( v(i) + ss2 ) )
    
    
    ee(1 : p1, 1) = ee(1 : p1, 1) +  z1 * bh( ss1 ) * x(i, :) - &
                    ( 1 - s(i) ) * hf( ss1 ) * x(i, :)
                
    
    
    ee((p1 + 1) : p, 1) =  ee((p1 + 1) : p, 1) + z2 * &
                           bh( ss1 ) * w(i, :)
    
                           
                                
    xi(:, 1) = x(i, :)
    wi(:, 1) = w(i, :)
                                
    
    gra(1 : p1, 1 : p1) = gra(1 : p1, 1 : p1) + z1 * dbh( ss1 ) * &
                          matmul(xi, transpose(xi)) - ( 1 - s(i) ) * & 
                          dhf( ss1 ) * matmul(xi, transpose(xi))
                          
    gra((p1 + 1) : p, 1 : p1) = gra((p1 + 1) : p, 1 : p1) + &
                                        z2 * dbh( ss1 ) * &
                                        matmul(wi, transpose(xi))
                                        
    
    
    
    
    gra(1 : p1, (p1 + 1) : p) = gra(1 : p1, (p1 + 1) : p) + &
                                z3 * matmul(xi, transpose(wi))
                                        
                                        
    
                          
    z4 = - s(i) * r(i) * dhf( v(i) + ss2 ) * bh( ss1 )
                 
      
    gra((p1 + 1) : p, (p1 + 1) : p ) = gra((p1 + 1) : p, (p1 + 1) : p ) +  &
                                       z4 * matmul(wi, transpose(wi))
    
  
  end do
  
  call inverse(gra, invgra, p)

  dir = matmul(invgra, ee)
  par = intval - dir(:, 1)

  j = j + 1
  diff = sum(abs(par - intval))
  intval = par

end do

abetainv = invgra 

return
end subroutine eesolver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nadaraya-Watson estimator with the bandwidth chosen by the leave-one-out cross 
! validation

subroutine nwecv(fit, x1, y1, x2, bset, m, n1, n2, p)
implicit none

integer*4 :: m, n1, n2, p
real*8 :: fit(n2), x1(n1 * p), y1(n1), x2(n2 * p), bset(m), b

call cv(b, bset, m, x1, y1, n1, p)

call nwe(fit, x1, y1, x2, b, n1, n2, p)


return
end subroutine nwecv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! leave-one-out cross validation for bandwidth
! The criterion is maximum likelihood.

subroutine cv(choice, bset, m, x, y, n, p)
implicit none

integer*4 :: m, n, p, i, j, ind(n - 1), z, location
real*8 :: choice, bset(m), x(n,p), y(n), x1(n-1, p), y1(n-1), x2(1, p), y2, b, &
          fit, ll(m), out(1)

ll = 0.d0

do i = 1, m
  
  b = bset(i)

  x1 = x(2:n, :)
  y1 = y(2:n)
  x2(1, :) = x(1, :)
  y2 = y(1)

  call nwe(out, x1, y1, x2, b, n - 1, 1, p)
  fit = out(1)
  

  ll(i) = ll(i) + y2 * log(fit) + (1.d0 - y2) * log(1.d0 - fit)
  
  
  do j = 2, (n-1)

    ind = (/ (z, z = 1, j-1), (z, z = j+1, n)/)
    
    x1 = x(ind, :)
    y1 = y(ind)
    x2(1, :) = x(j,:)
    y2 = y(j)

    call nwe(out, x1, y1, x2, b, n - 1, 1, p)
    fit = out(1)
  

    ll(i) = ll(i) + y2 * log(fit) + (1.d0 - y2) * log(1.d0 - fit)
  
  end do
  
  
  x1 = x(1:(n-1), :)
  y1 = y(1:(n-1))
  x2(1, :) = x(n, :)
  y2 = y(n)

  call nwe(out, x1, y1, x2, b, n - 1, 1, p)
  fit = out(1)
  

  ll(i) = ll(i) + y2 * log(fit) + (1.d0 - y2) * log(1.d0 - fit)


end do

location = maxloc(ll, dim = 1)
choice = bset(location)

return
end subroutine cv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nadaraya watson estimator

subroutine nwe(fit, x1, y1, x2, b, n1, n2, p)
implicit none

integer*4 :: n1, n2, p, i, j
real*8 :: x1(n1,p), y1(n1), x2(n2,p), b, dist, weight(n1), fit(n2), wa

do i = 1, n2

  do j = 1, n1

    dist = sum((x1(j,:) - x2(i,:)) ** 2)
    weight(j) = exp(-dist / (2.d0 * (b ** 2)))/((2.d0 * 3.1415926) ** (p / 2.d0))

  end do
  
  wa = sum(weight*y1)/sum(weight)
  
  if (wa == 0.d0) then
     
    wa = wa + 10.d0 ** (-4)
  
  end if
  
  
  if (wa == 1.d0) then
     
    wa = wa - 10.d0 ** (-4)
  
  end if

  fit(i) = wa

end do

return
end subroutine nwe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! matrix inverse
subroutine inverse(a,c,n)
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n), z(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
    z(i) = x(i)
    if(n==24) then
      z(i) = x(i)/sqrt(2.d0)
    end if
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = z(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! logistic function

function hf(xx) result(yy)

real*8, intent(in) :: xx
real*8 :: yy

yy = 1 / (1 + exp(-xx))

end function hf




! 1 - logistic function

function bh(xx) result(yy)

real*8, intent(in) :: xx
real*8 :: yy

yy = 1 / (1 + exp(xx))

end function bh




! derivative of logistic function

function dhf(xx) result(yy)

real*8, intent(in) :: xx
real*8 :: yy

yy = exp(xx) / ((1 + exp(xx)) ** 2)

end function dhf




! derivative of ( 1 -  logistic function )

function dbh(xx) result(yy)

real*8, intent(in) :: xx
real*8 :: yy

yy = - exp(xx) / ((1 + exp(xx)) ** 2)

end function dbh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
