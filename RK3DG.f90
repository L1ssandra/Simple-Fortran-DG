    subroutine RK3DG(Pol,PolValue,PolDx,PolDy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
    
    integer i,j,Nx,Ny,dim,dimpol,quad
    parameter(quad = 4)
    real Pol(Nx,Ny,dimpol,dim),PolValue(quad,quad,dimpol,dim),PolDx(quad,quad,dimpol,dim),PolDy(quad,quad,dimpol,dim)
    real Poledge(4,quad,dimpol,dim),Q(Nx,Ny,dim),diagpol(Nx,Ny,dimpol,dim),Polone(Nx,Ny,dimpol,dim),Poltwo(Nx,Ny,dimpol,dim)
    real LPol(Nx,Ny,dimpol,dim),lambdai(quad),lambdaj(quad),weight(quad)
    
    call Lh(Pol,LPol,PolValue,PolDx,PolDy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
    
    Polone = Pol + dt*LPol
    
    call Lh(Polone,LPol,PolValue,PolDx,PolDy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
    
    Poltwo = 0.75d0*Pol + 0.25d0*Polone + 0.25d0*dt*LPol
    
    call Lh(Poltwo,LPol,PolValue,PolDx,PolDy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
    
    Pol = (1d0/3d0)*Pol + (2d0/3d0)*Poltwo + (2d0/3d0)*dt*LPol
    
    call Pol_to_Q(Pol,Q,Nx,Ny,dimpol,dim)
    
    end subroutine RK3DG
    
    
    
    subroutine Lh(Pol,DPol,PolValue,PolDx,PolDy,Poledge,Q,hx,hy,dt,Nx,Ny,dim,dimpol,diagpol,lambdai,lambdaj,weight)
    
    integer i,j,Nx,Ny,dim,dimpol,quad,i1,j1,d1
    parameter(quad = 4)
    real Pol(Nx,Ny,dimpol,dim),PolValue(quad,quad,dimpol,dim),PolDx(quad,quad,dimpol,dim),PolDy(quad,quad,dimpol,dim)
    real Poledge(4,quad,dimpol,dim),Q(Nx,Ny,dim),diagpol(Nx,Ny,dimpol,dim),DPol(Nx,Ny,dimpol,dim),alphax,alphay,alphatest
    real PolR(Nx,Ny,dimpol,dim),PolL(Nx,Ny,dimpol,dim),PolU(Nx,Ny,dimpol,dim),PolD(Nx,Ny,dimpol,dim)
    real u,v,p,rho,gamma
    real IV(dim),IR(dim),IL(dim),IU(dim),ID(dim),lambdai(quad),lambdaj(quad),weight(quad),fu(dim,2)
    real uint(dim),uext(dim),hLF(dim),fuint(dim),fuext(dim),hx1,hy1
    integer :: boundconditionX = 1
    integer :: boundconditionY = 1
    
    call Pol_to_Q(Pol,Q,Nx,Ny,dimpol,dim)
    
    gamma = 1.4d0
    
    alphax = 1e-3
    alphay = 1e-3
    
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    
    ! 计算alphax,alphay，使用全局LF通量
    do i = 1,Nx
        do j = 1,Ny
            
            rho = Q(i,j,1)
            u = Q(i,j,2)/Q(i,j,1)
            v = Q(i,j,3)/Q(i,j,1)
            p = (gamma - 1)*(Q(i,j,4) - 0.5d0*rho*(u**2 + v**2))
            
            call max_wavespeed(alphatest,rho,u,v,p,gamma,1)
            
            if (alphatest > alphax) then
                alphax = alphatest
            end if
            
            call max_wavespeed(alphatest,rho,u,v,p,gamma,2)
            
            if (alphatest > alphay) then
                alphay = alphatest
            end if
            
        end do
    end do
    
    ! 赋边界条件
    if (boundconditionX == 1) then ! 周期边界
        PolR(1:Nx - 1,:,:,:) = Pol(2:Nx,:,:,:)
        PolR(Nx,:,:,:) = Pol(1,:,:,:)
        
        PolL(2:Nx,:,:,:) = Pol(1:Nx - 1,:,:,:)
        PolL(1,:,:,:) = Pol(Nx,:,:,:)
    end if
    
    if (boundconditionY == 1) then ! 周期边界    
        PolU(:,1:Ny - 1,:,:) = Pol(:,2:Ny,:,:)
        PolU(:,Ny,:,:) = Pol(:,1,:,:)
        
        PolD(:,2:Ny,:,:) = Pol(:,1:Ny - 1,:,:)
        PolD(:,1,:,:) = Pol(:,Ny,:,:)
    end if
        
    
    ! 迭代
    do i = 1,Nx
        do j = 1,Ny
            do d = 1,dimpol
                
                ! 体积分 \int_K (f(u)·▽v) dx
                IV = 0
                do i1 = 1,quad
                    do j1 = 1,quad
                        
                        ! 计算当前积分点上的解
                        uint = 0
                        do d1 = 1,dimpol
                            uint = uint + Pol(i,j,d1,:)*PolValue(i1,j1,d1,:)
                        end do
                        
                        call f(uint,fu)
                        
                        IV = IV + weight(i1)*weight(j1)*(fu(:,1)*PolDx(i1,j1,d,:) + fu(:,2)*PolDy(i1,j1,d,:))
                        
                    end do
                end do
                IV = hx1*hy1*IV
                
                ! 面积分 \int_(partial K) (h(uint,uext,n)·v) dSx
                
                ! 右
                IR = 0
                do j1 = 1,quad
                    
                    ! 计算当前积分点上内部和外部的解
                    uint = 0
                    uext = 0
                    do d1 = 1,dimpol
                        uint = uint + Pol(i,j,d1,:)*Poledge(1,j1,d1,:)
                        uext = uext + PolR(i,j,d1,:)*Poledge(2,j1,d1,:)
                    end do
                    
                    call f1(uint,fuint)
                    call f1(uext,fuext)
                    
                    hLF = 0.5d0*(fuint + fuext - alphax*(uext - uint))
                    
                    IR = IR + weight(j1)*hLF*Poledge(1,j1,d,:)
                    
                end do
                IR = hx1*IR
                
                ! 左
                IL = 0
                do j1 = 1,quad
                    
                    uint = 0
                    uext = 0
                    do d1 = 1,dimpol
                        uint = uint + Pol(i,j,d1,:)*Poledge(2,j1,d1,:)
                        uext = uext + PolL(i,j,d1,:)*Poledge(1,j1,d1,:)
                    end do
                    
                    call f1(uint,fuint)
                    call f1(uext,fuext)
                    
                    hLF = 0.5d0*(-(fuint + fuext) - alphax*(uext - uint))
                    
                    IL = IL + weight(j1)*hLF*Poledge(2,j1,d,:)
                    
                end do
                IL = hx1*IL
                
                ! 上
                IU = 0
                do i1 = 1,quad
                    
                    uint = 0
                    uext = 0
                    do d1 = 1,dimpol
                        uint = uint + Pol(i,j,d1,:)*Poledge(3,i1,d1,:)
                        uext = uext + PolU(i,j,d1,:)*Poledge(4,i1,d1,:)
                    end do
                    
                    call f2(uint,fuint)
                    call f2(uext,fuext)
                    
                    hLF = 0.5d0*(fuint + fuext - alphay*(uext - uint))
                    
                    IU = IU + weight(i1)*hLF*Poledge(3,i1,d,:)
                    
                end do
                IU = hy1*IU
                
                ! 下
                ID = 0
                do i1 = 1,quad
                    
                    uint = 0
                    uext = 0
                    do d1 = 1,dimpol
                        uint = uint + Pol(i,j,d1,:)*Poledge(4,i1,d1,:)
                        uext = uext + PolD(i,j,d1,:)*Poledge(3,i1,d1,:)
                    end do
                    
                    call f2(uint,fuint)
                    call f2(uext,fuext)
                    
                    hLF = 0.5d0*(-(fuint + fuext) - alphay*(uext - uint))
                    
                    ID = ID + weight(i1)*hLF*Poledge(4,i1,d,:)
                    
                end do
                ID = hy1*ID
                
                DPol(i,j,d,:) = IV - IR - IL - IU - ID
                
            end do
        end do
    end do
    DPol = DPol/diagpol
            
    
    end subroutine Lh