! CALCULATE BLOCK Central Force (only spring force and body force)
subroutine Get_Spring_Force( UU, tn, Fp, VVs )
	! use MSIMSLMD
	! use NUMERICAL_LIBRARIES
	use global_data
	implicit none
	
	integer :: i 
	real*8 :: UU( 3*TBN );		! Input 
	real*8 :: tn				! Input : time
	
	real*8 :: Fp( 3*(TBN-1) );	! Output : Spring Block Central Force
	real*8 :: VVs				! Output : spring potential energy
		
	real*8 :: xg, xgd, xg2d			
	integer :: bn1, bn2, ebn, vn1, vn2
	real*8 :: c, s, Rt1(2,2), Rt2(2,2), Xt1(2), Xt2(2);
	real*8 :: rr(2), r1_(2), r2_(2), dx, dy, theta
	real*8 :: M1, M2, Fx, Fy
	real*8 :: Tv(2), Nv(2), Rt(2,2), KS(2,2), FF(2);

	call getGround(tn, xg, xgd, xg2d);
	
	UU(1) = xg;
	
	Fp = 0;
	VVs = 0;
		
	! Calculate Block Central Force
	do i=1,cpN
		
		!! BLOCK I
		bn1 = cpdata(i)%bn(1);
		theta = UU(3*bn1);		
		
		c = dcos( theta );		s = dsin( theta );
		Rt1(1,1) = c;			Rt1(1,2) = -s;
		Rt1(2,1) = s;			Rt1(2,2) = c;  
		
		rr = cpdata(i)%p(1,:);
		rr = matmul( Rt1, rr );
		r1_ = (/-rr(2), rr(1)/);

		!! =========  Disp. & Vel. of contact point of block 1
		Xt1 = UU(3*bn1-2:3*bn1-1) + rr;

		!! BLOCK J
		bn2 = cpdata(i)%bn(2);
		theta = UU(3*bn2);		
		
		c = dcos( theta );		s = dsin( theta );
		Rt2(1,1) = c;			Rt2(1,2) = -s;
		Rt2(2,1) = s;			Rt2(2,2) = c;  

		rr = cpdata(i)%p(2,:);
		rr = matmul( Rt2, rr );
		r2_ = (/-rr(2), rr(1)/);

		!! =========  Disp. of contact point of block 2
		Xt2 = UU(3*bn2-2:3*bn2-1) + rr;	  ! Xc + r
		
		!! ###########################################################################
		! Slip 일 경우
		if( cpdata(i)%state ) then		
		!! Get the Tangential & Normal vector of contact point				
			ebn = cpdata(i)%en(1);	! edge block number
			vn1 = cpdata(i)%en(2);	! begin vertex number of edge block
			vn2 = vn1 + 1;			! end vertex number of edge block

			if( vn2 > Block(ebn)%vN ) then
				vn2 = 1;
			end if

			! Tangential Vector
			Tv = Block(ebn)%vtx(vn1,:) - Block(ebn)%vtx(vn2,:);
			Tv = Tv / dsqrt( Tv(1)**2 + Tv(2)**2 );
				
			! Block i (bn1) 에 대한 tangential vector 로 바꿔 줌. 
			if( ebn == cpdata(i)%bn(1) ) then				
				Tv = matmul(Rt1,Tv);
			else				
				Tv = -matmul(Rt2,Tv);
			end if 
							
			! Block i(1)에 대한 Normal Vector
			Nv = (/-Tv(2), Tv(1)/);

			! Transformation Matrix
			Rt(1,:) = Tv;	
			Rt(2,:) = Nv;
			
			KS = 0;
			KS(1,2) = cpdata(i)%sign*mu*kp;
			KS(2,2) = kp;

			KS = matmul( transpose(Rt), matmul(KS,Rt) );
						
		end if
				
				
		!! Calculate Spring Contact Force =======
		dx = Xt1(1) - Xt2(1);
		dy = Xt1(2) - Xt2(2);

		if( cpdata(i)%state ) then	! slip 이면						
			FF = -matmul(KS,(/dx, dy/));			
			Fx = FF(1);    	Fy = FF(2);	
		else
			Fx = -kp*dx;		
			Fy = -kp*dy;
		end if

		if(0 .and. cpdata(i)%sep_bit ) then
			Fx = 0;	  Fy = 0;
		end if
		
		! Spring potential energy
		VVs = VVs - 0.5*Fx*dx - 0.5*Fy*dy;	

		M1 = Fx * r1_(1) + Fy * r1_(2);
		M2 = Fx * r2_(1) + Fy * r2_(2);

		! Block i
		bn1=bn1-1; ! 자유도를 맞추기 위해 1을 빼줌
		if( bn1 /= 0 ) then
			Fp(3*bn1-2) = Fp(3*bn1-2) + Fx;
			Fp(3*bn1-1) = Fp(3*bn1-1) + Fy;
			Fp(3*bn1) = Fp(3*bn1) + M1;
		end if																	  

		! Block j
		bn2=bn2-1;
		if( bn2 /= 0 ) then			
			Fp(3*bn2-2) = Fp(3*bn2-2) - Fx;
			Fp(3*bn2-1) = Fp(3*bn2-1) - Fy;
			Fp(3*bn2) = Fp(3*bn2) - M2;
		end if
		!! ===========================================================

	end do

	!! Calculate Block Central Force
	do i=1, TBN-1
		Fp(3*i-1) = Fp(3*i-1) - MM(3*i-1,3*i-1)*g;
	end do
	!! ====================
	
end subroutine

! Newmark Method
subroutine Newmark( KK, RP, CC, Ut, Udt, U2dt, UU, UUd, UU2d )
	use global_data
	implicit none

	real*8 :: KK( 3*(TBN-1), 3*(TBN-1) ), RP( 3*(TBN-1) ), CC( 3*(TBN-1), 3*(TBN-1) );
	real*8, dimension(3*(TBN-1)) :: Ut, Udt, U2dt;

	real*8, dimension(3*(TBN-1)) :: UU, UUd, UU2d;

	real*8 :: alpha, delta
	real*8 :: a0, aa1, a2, a3, a4, a5, a6, a7
	real*8 :: K_hat( 3*(TBN-1), 3*(TBN-1) );
	real*8 :: R_hat( 3*(TBN-1) )

	alpha = 0.25d0;
	delta = 0.5d0;

	a0 = 1.0/alpha/dt**2;
	aa1 = delta/alpha/dt;
	a2 = 1.0/alpha/dt;
	a3 = 1.0/2.0/alpha - 1;
	a4 = delta/alpha - 1;
	a5 = dt/2.0 * (delta/alpha - 2 );
	a6 = dt*(1-delta);
	a7 = delta*dt;

	K_hat = KK + a0*MM + aa1*CC;
	R_hat = RP + matmul(MM, a0*Ut + a2*Udt + a3*U2dt ) + matmul(CC, aa1*Ut + a4*Udt + a5*U2dt );

	call DLSARG(3*(TBN-1),K_hat,3*(TBN-1),R_hat,1,UU)
	UU2d = a0*(UU - Ut) - a2*Udt - a3*U2dt;
	UUd =  Udt + a6*U2dt + a7*UU2d;
	
end subroutine


! CONSTRUCT MASS MATRIX
subroutine Construct_Mass
	use global_data
	implicit none

	integer :: i

	allocate( MM( 3*(TBN-1), 3*(TBN-1) ), FC( 3*(TBN-1) ) );
	
	MM = 0;
	do i=1,TBN-1
		MM(3*i-2, 3*i-2) = Block(i+1)%M;
		MM(3*i-1, 3*i-1) = Block(i+1)%M;
		MM(3*i, 3*i) = Block(i+1)%I;
	end do

end subroutine

! CONSTRUCT DAMPING AND GET THE a1 VALUE
subroutine Construct_Stiffness_Damping(UU, KK, Kfg, KKc)	
	use global_data
	implicit none
	
	real*8 :: UU( 3*TBN )		! Input

	real*8 :: KK(3*(TBN-1), 3*(TBN-1))	! Output	
	real*8 :: KKc(3*(TBN-1), 3*(TBN-1))	! Output
	real*8 :: Kfg(3*(TBN-1), 3)		! Output, K_fs matrix ...

	real*8 :: KG(3*TBN, 3*TBN)		! global stiffness matrix
	real*8 :: KCG(3*TBN, 3*TBN)		! global damping matrix
	integer :: n, i, j, bn1, bn2, ebn, vn1, vn2	
	real*8 :: theta, rct(2), r_(2);
	real*8 :: Tv(2), Nv(2), Rt(2,2), R(2,2), K(2,2), c, s, cl(2,3), Kb(3,3)	
								  
		
	cl = 0;		cl(1,1) = 1;	cl(2,2) = 1;
	
	KG = 0;		KK = 0;		KCG=0; KKc = 0;
	do n=1, cpN

		K = 0;		K(1,1) = kp;	K(2,2) = kp;

		! Slip 이면 K 를 따로 구성해 준다.
		if( cpdata(n)%state ) then

			!! ###########################################################################
			!! Get the Tangential & Normal vector of contact point				
			ebn = cpdata(n)%en(1);	! edge block number
			vn1 = cpdata(n)%en(2);	! begin vertex number of edge block
			vn2 = vn1 + 1;			! end vertex number of edge block

			if( vn2 > Block(ebn)%vN ) then
				vn2 = 1;
			end if

			! Tangential Vector
			Tv = Block(ebn)%vtx(vn1,:) - Block(ebn)%vtx(vn2,:);
			Tv = Tv / dsqrt( Tv(1)**2 + Tv(2)**2 );
				
			! Block i (bn1) 에 대한 tangential vector 로 바꿔 줌. 
			if( ebn == cpdata(n)%bn(1) ) then
				! block 1 에 관한것들
				theta = UU( 3*cpdata(n)%bn(1) );
				c = dcos(theta);	s = dsin(theta);
				Rt(1,1) =  c;		Rt(1,2) = -s;
				Rt(2,1) =  s;		Rt(2,2) =  c;
				Tv = matmul(Rt,Tv);
			else
				! block 2 에 관한것들
				theta = UU( 3*cpdata(n)%bn(2) );
				c = dcos(theta);	s = dsin(theta);
				Rt(1,1) =  c;		Rt(1,2) = -s;
				Rt(2,1) =  s;		Rt(2,2) =  c;
				Tv = -matmul(Rt,Tv);
			end if 
							
			! Block i(1)에 대한 Normal Vector
			Nv = (/-Tv(2), Tv(1)/);

			! Transformation Matrix
			Rt(1,:) = Tv;	
			Rt(2,:) = Nv;
			
			K = 0;
			K(1,2) = cpdata(n)%sign*mu*kp;
			K(2,2) = kp;

			K = matmul( transpose(Rt), matmul(K,Rt) );

		end if

		if( cpdata(n)%sep_bit ) then
			K = 0;
		end if
						
		!!!!!!!!!!!!!!!!!!!! Kl(6x6) matrix arrange		
		do i=1,2
			do j=1,2				
				
				! j block에 관한것들
				theta = UU( 3*cpdata(n)%bn(j) );
				c = dcos(theta);	s = dsin(theta);
				R(1,1) =  c;		R(1,2) = -s;
				R(2,1) =  s;		R(2,2) =  c;

				rct = matmul(R,cpdata(n)%p(j,:));
				r_(1) = -rct(2);	r_(2) = rct(1); 
				cl(:,3) = r_;

				! i block에 관한 것들
				theta = UU( 3*cpdata(n)%bn(i) );
				c = dcos(theta);	s = dsin(theta);
				R(1,1) =  c;		R(1,2) = -s;
				R(2,1) =  s;		R(2,2) =  c;

				rct = matmul(R,cpdata(n)%p(i,:));
				r_(1) = -rct(2);	r_(2) = rct(1); 
				
				! K_bar			   				
				Kb(1:2,1:3) = matmul(K,cl);
				Kb(3,:) = matmul(r_,matmul(K,cl));

				if( i /= j ) then
					Kb = -Kb;
				end if

				! Global Assemble
				bn1 = cpdata(n)%bn(i);		bn2 = cpdata(n)%bn(j);
				KG( 3*(bn1-1)+1:3*bn1, 3*(bn2-1)+1:3*bn2 ) = KG( 3*(bn1-1)+1:3*bn1, 3*(bn2-1)+1:3*bn2 ) + Kb;

				if( cpdata(n)%sep_bit /= 1 ) then
					KCG( 3*(bn1-1)+1:3*bn1, 3*(bn2-1)+1:3*bn2 ) = KCG( 3*(bn1-1)+1:3*bn1, 3*(bn2-1)+1:3*bn2 ) + Kb;
				end if	
			
			end do
		end do

	end do
	
	KK = KG(4:3*TBN, 4:3*TBN);
	KKc = KCG(4:3*TBN, 4:3*TBN);
	Kfg = KG(4:3*TBN, 1:3 );

end subroutine

! Get a1 value for the Stiffness Proportional Damping
subroutine Get_a1(KK, CC)	
	use global_data
	implicit none
	
	real*8 :: KK(3*(TBN-1), 3*(TBN-1))	! Input
	real*8 :: CC(3*(TBN-1), 3*(TBN-1))	! Output

	integer :: i
	real*8 :: evalt(3*(TBN-1)), eval( 3*(TBN-1) ), beta( 3*(TBN-1) ), xin(3*(TBN-1)), ev(3*(TBN-1), 3*(TBN-1));
	double complex :: alpha( 3*(TBN-1) ), EVEC(3*(TBN-1), 3*(TBN-1));

	real*8 :: Mn( 3*(TBN-1), 3*(TBN-1) ), Cn( 3*(TBN-1), 3*(TBN-1) ), Mn_i( 3*(TBN-1), 3*(TBN-1) );
	real*8 :: w1

	! for 2 dimensional analysis
	integer :: j
	real*8 :: eval2( 2*(TBN-1) ), K2( 2*(TBN-1), 2*(TBN-1) );
	real*8 :: M2( 2*(TBN-1), 2*(TBN-1) ), a2, ff(2*(TBN-1));


	
				
	!! GET EIGENVALUES  
    CALL DGVCRG (3*(TBN-1), KK, 3*(TBN-1), MM, 3*(TBN-1), ALPHA, BETA, EVEC, 3*(TBN-1))

	ev = dreal(EVEC);	
	do i=1,3*(TBN-1) 
		eval(i) = dabs( dreal(ALPHA(i)/BETA(i)) ) ;
		eval(i) = dsqrt(eval(i));				
	end do
	
	evalt = eval;

	!##########################################################3
	! 3 dof
	!! SORT THE EIGENVALUES
	call DSVRGN (3*(TBN-1), eval, evalt);
			
	do i=1, 3*(TBN-1)	!! 0이 아닌 가장 작은값 1개를 고른다.	
		if( evalt(i) > 1 ) then
			w1 = evalt(i);
			goto 100;
		end if
	end do

100	continue;
	
	a1 = xi*2/w1;		
	!##########################################################3

	!print *, "w1 =", w1
	
	
	if(0) then
	do i=1,3*(TBN-1) 
		
		if( eval(i) < 1 ) then
			xin(i) = 0;	
		else
			!xin(i) = a1 * eval(i) / 2;
			xin(i) = xi;
		end if
	
	end do
	
	!xin(1) = 2;

	!if( cpN == 2 ) then
	!	xin(3) = 9.46; xin(2) = 4.75; xin(3) = 1;		
	!end if

	!print *," w1 = ",w1		

	!##########################################################3
	! modal diaonal matrix of mass and damping	
	Mn = matmul(matmul( transpose(ev), MM ), ev);
	Cn = 0;
	Mn_i = 0;	
	do i=1,3*(TBN-1)
		Cn(i,i) = 2*xin(i)*eval(i)*Mn(i,i);
		Mn_i(i,i) = 1/Mn(i,i);
	end do
	
	! damping matrix
	CC = matmul( matmul(MM,ev), Mn_i );
	CC = matmul(CC, Cn);
	CC = matmul( CC, matmul( matmul(Mn_i,transpose(ev)), MM) ); 
	end if
	!##########################################################3


	!##############################################################
	!! 2 Dimensional
if(0) then
	j=1;
	do i=1,3*(TBN-1)
		if( mod(i,3) /= 0 ) then
			ff(j) = i;
			j=j+1;
		end if
	end do

	do i=1,2*(TBN-1)
		do j=1,2*(TBN-1)
			K2(i,j) = KK(ff(i),ff(j));
			M2(i,j) = MM(ff(i),ff(j));
		end do
	end do

	!! GET EIGENVALUES  
	call DGVLRG (2*(TBN-1), K2, 2*(TBN-1), M2, 2*(TBN-1), ALPHA, BETA)
    	
	do i=1,2*(TBN-1) 
		eval2(i) = dabs( dreal(ALPHA(i)/BETA(i)) ) ;
		eval2(i) = dsqrt(eval2(i));		
	end do
		

	!! SORT THE EIGENVALUES
	call DSVRGN (2*(TBN-1), eval2, eval2);


	do i=1, 2*(TBN-1)	!! 0이 아닌 가장 작은값 2개를 고른다.	
		if( eval2(i) > 1 ) then
			w1 = eval2(i);

			goto 160;
		end if
	end do

160	continue;

	a2 = 2*xi/w1;								  
	!a1 = a2;

end if
	!##############################################################




					
end subroutine

! GET THE GROUND MOTION DATA
subroutine getGround(tn, xg, xgd, xg2d)
	use global_data
	implicit none

	real*8 :: tn				! input
	real*8 :: xg, xgd, xg2d		! output

	real*8 :: tau	! 0 <= tau <= adt
	integer :: idnum
	real*8 :: u2d_i, u2d_i1

	
	! HARMORNIC MOTION INPUT TYPE
	if( input_type == 1 ) then
		xg = -Amp/((2*pi*frq)**2) * dsin( 2*pi*frq*tn );
		xgd = -Amp/(2*pi*frq) * dcos( 2*pi*frq*tn );
		xg2d = Amp * dsin( 2*pi*frq*tn );
	end if

	! EARTHQUAKE MOTION INPUT TYPE
	if( input_type == 2 ) then
	
		if( 0 .and. tn > 0.09998 ) then
			print *, "asdfasdfasdf"
		end if
		
		idnum = int(tn/adt);
		tau = tn - adt*(idnum);
		
		u2d_i = dacc(idnum+1);
		u2d_i1 = dacc(idnum+2);

		if( prev_idnum < idnum ) then
			prev_u_i = u_i;
			prev_ud_i = ud_i;
			
			u2d_i = dacc(idnum);
			u2d_i1 = dacc(idnum+1);

			u_i = u_i + ud_i * adt + adt**2 * ( u2d_i1/6 + u2d_i/3 );
			ud_i = ud_i + adt/2 * (u2d_i1 + u2d_i);
			
			prev_idnum = idnum;


		end if

		if( prev_idnum > idnum ) then
			u_i = prev_u_i;
			ud_i = prev_ud_i;
			prev_idnum = idnum;
		end if


		u2d_i = dacc(idnum+1);
		u2d_i1 = dacc(idnum+2);

		xg2d = u2d_i + tau/adt * (u2d_i1 - u2d_i);
		xgd = ud_i + u2d_i * tau + tau**2 / 2 / adt * (u2d_i1 - u2d_i);
		xg = u_i + ud_i * tau + u2d_i * tau**2 / 2 + tau**3 / 6 / adt * (u2d_i1 - u2d_i);
	
	end if



end subroutine


! Check the Stick or Slip or Separation &&&& Update Contact Data
subroutine Check_Stick_Slip_Separation( tn, UU, UUd, UU2d )
	use global_data
	implicit none

	real*8 :: UU( 3*TBN ), UUd( 3*TBN ), UU2d( 3*TBN ) ! Input 
	real*8 :: tn				! Input : time	

	byte :: cbit, cbit1, cbit2, stick_bit
	integer :: i, j, k, copy_cpN
	integer :: spN, spidx(cpN+1)  ! separation 개수 , separation index number ( program 참조 )
	integer :: vbn, v1, v2
	
	integer :: ebn, bn1, bn2, vn1, vn2, bbn1, bbn2
	real*8 :: c, s, Rt(2,2),Rt1(2,2), Rt2(2,2), Xt1(2), Xt2(2), Vt1(2), Vt2(2);
	real*8 :: At1(2), At2(2), dat(2);
	real*8 :: rr(2), r1_(2), r2_(2), dx, dy, theta
	real*8 :: Tv(2), Nv(2), lv(2), dxL, dyL, dvt(2)	
	real*8 :: FN, FT
	real*8 :: xg, xgd, xg2d
	!real*8 :: c_vel(cpN+1)
	!integer :: factor
		
	!c_vel = 1d-6; 
	call getGround(tn, xg, xgd, xg2d);
	
	UU(1) = xg;			
	UU(2) = Block(1)%Mc(2);	
	UUd(1) = xgd;
	UU2d(1) = xg2d;
	!Ug(4:3*TBN) = UU;	
	
	spN = 0;	spidx = 0;   ! variables for separation

	copy_cpN = cpN;		! 중간에 cpN이 변하는 경우가 있기 때문에..

	do i=1,copy_cpN
		
		!! BLOCK I
		bn1 = cpdata(i)%bn(1);
		theta = UU(3*bn1);		
		
		c = dcos( theta );		s = dsin( theta );
		Rt1(1,1) = c;			Rt1(1,2) = -s;
		Rt1(2,1) = s;			Rt1(2,2) = c;  

		rr = cpdata(i)%p(1,:);
		rr = matmul( Rt1, rr );
		r1_ = (/-rr(2), rr(1)/);

		!! =========  Disp. of contact point of block 1
		Xt1 = UU(3*bn1-2:3*bn1-1) + rr;		

		!! BLOCK J
		bn2 = cpdata(i)%bn(2);
		theta = UU(3*bn2);		
		
		c = dcos( theta );		s = dsin( theta );
		Rt2(1,1) = c;			Rt2(1,2) = -s;
		Rt2(2,1) = s;			Rt2(2,2) = c;  

		rr = cpdata(i)%p(2,:);
		rr = matmul( Rt2, rr );
		r2_ = (/-rr(2), rr(1)/);

		!! =========  Disp. of contact point of block 2
		Xt2 = UU(3*bn2-2:3*bn2-1) + rr;	  ! Xc + r		
				
		!! ###########################################################################
		!! Get the Tangential & Normal vector of contact point
				
		ebn = cpdata(i)%en(1);	! edge block number
		vn1 = cpdata(i)%en(2);	! begin vertex number of edge block
		vn2 = vn1 + 1;			! end vertex number of edge block

		if( vn2 > Block(ebn)%vN ) then
			vn2 = 1;
		end if

		! Tangential Vector
		Tv = Block(ebn)%vtx(vn1,:) - Block(ebn)%vtx(vn2,:);
		Tv = Tv / dsqrt( Tv(1)**2 + Tv(2)**2 );
		
		! Block i (bn1) 에 대한 tangential vector 로 바꿔 줌. 
		if( ebn == bn1 ) then
			Tv = matmul(Rt1,Tv);
		else
			Tv = -matmul(Rt2,Tv);
		end if 
							
		! Normal Vector
		Nv = (/-Tv(2), Tv(1)/);
				
		!! ###########################################################################
		!! Calculate Tangential & Normal Contact force of Block 1 (i)
		dx = Xt1(1) - Xt2(1);	! global
		dy = Xt1(2) - Xt2(2);	! global

		dxL = Tv(1)*dx + Tv(2)*dy;	! tangential disp.
		dyL = Nv(1)*dx + Nv(2)*dy;	! Normal disp.
		

		Vt1	= UUd(3*bn1-2:3*bn1-1) + UUd(3*bn1)*r1_; 
		Vt2	= UUd(3*bn2-2:3*bn2-1) + UUd(3*bn2)*r2_; 

		At1 = UU2d(3*bn1-2:3*bn1-1) + UU2d(3*bn1)*r1_; 
		At2	= UU2d(3*bn2-2:3*bn2-1) + UU2d(3*bn2)*r2_;

		Rt(1,:) = Tv;	Rt(2,:) = Nv;
		dvt = matmul(Rt,Vt1-Vt2);	! relative velocity
		dat = matmul(Rt,At1-At2);	! relative acc.
		
		if( 0 .and. tn > 0 ) then		
				
			if( 1 .and. cpdata(i)%vn(1) == 2 .and. cpdata(i)%vn(2) == 2 ) then
				open (unit=30, file="relv_out.dat", action="write", position="append" )
				write (30,15)  tn, dvt(2), dat(2), dyL, cpdata(i)%d_stable, cpdata(i)%stable
				15 FORMAT (5E14.7, I5)
				close(30);	
			end if

			!if( 1 .and. cpdata(i)%vn(1) == 4 .and. cpdata(i)%vn(2) == 2 ) then
			!	write (30,'(100E14.7)')  tn, dvt(2), dat(2), dyL, dyL, a1*kp*dvt(2)/(Block(2)%M*9.81);	
			!end	if
			
		
		end if
		
		! determine the contact stable state
		if( dabs( a1*kp*dvt(2)/cpdata(i)%light_mass/g ) < 0.01 .and. dabs( dat(2) ) < g ) then
			cpdata(i)%stable = 1;
			cpdata(i)%forced_contact = 0;
		end if

		! depth of the stable line
		if( cpdata(i)%stable ) then
			cpdata(i)%d_stable = dyL;
			cpdata(i)%p_stable = 1;
		end if
		
		! preventing the inordinate repeat of V.S and recontact
		if( cpdata(i)%sep_bit ) then
			if( dabs(cpdata(i)%prev_a) > g .and. dabs(dat(2)) > g ) then
				if( cpdata(i)%prev_a * dat(2) < 0 ) then
			
					cpdata(i)%sep_bit = 0;
					cpdata(i)%d_stable = dyL;
					cpdata(i)%forced_contact = 1;
					print *,"==> Forced Contact : Block=",cpdata(i)%vn(1), " Vertex=",cpdata(i)%vn(2);
					print *,"dyL=",dyL;
				end if
			end if
		end if
				
		! Checking the virtual separation
		if( 1 .and. cpdata(i)%prev_a /= -99999 .and. dabs(dat(2)) > g .and. cpdata(i)%prev_a /= 0 ) then			
			
			if( dabs( dat(2) / cpdata(i)%prev_a ) > 3 ) then								
				cpdata(i)%stable = 0;				
			end if
						
			if( cpdata(i)%stable == 0 .and. cpdata(i)%forced_contact == 0 ) then

				!	VIRTUAL SEPARATION OCCURS
				if( dat(2) < 0 .and. dvt(2) < 0 .and. 1*(dyL-cpdata(i)%d_stable) < 0 ) then
				!if( cpdata(i)%p_stable .and. dat(2) < 0 .and. dvt(2) < 0 .and. a1*kp*dabs(dvt(2))/cpdata(i)%light_mass/g  > 0.01  ) then
					cpdata(i)%p_stable = 0;
					if( cpdata(i)%sep_bit == 0 ) then
						print *,"Block=",cpdata(i)%vn(1), " Vertex=",cpdata(i)%vn(2);
						new_contact = 1;
						slip_bit = 1;  ! ==> 그냥 test 용임... 나중에 지우기 바람...
					end if

					cpdata(i)%sep_bit = 1;
					!back_bit = 1;										
					!write (30,*)  "asdfasdfasdf";

				end if

			end if

		end if			
		
		
		if( cpdata(i)%sep_bit .and. dvt(2) > 0 ) then
			
			cpdata(i)%sep_bit = 0;
			cpdata(i)%d_stable = dyL;	  
			print *," ====> V.S. is recontacted =>",cpdata(i)%vn;
			new_contact = 1;
			!back_bit = 1;
		end if

		cpdata(i)%prev_a = dat(2); 		

!		if( dvt(1) > 1d-3 .and. cpdata(i)%state == 0 ) then
!			 print *,"dxL=", dxL
!		end if

		!c_vel(i) = dabs( dvt(2) );	! normal 방향의 속도를 저장해서 가장 큰 걸로 새로운 dt를 결정.
  		
		FT = -kp*dxL !- a1*kp*dvt(1);	! block i 에 작용하는 힘
		FN = -kp*dyL !- a1*kp*dvt(2);
		
		!!!!! block 번호가 작은 걸 앞에 배치한다.
		if( bn1 < bn2 ) then
			bbn1 = bn1;		bbn2 = bn2;
		else
			bbn1 = bn2;		bbn2 = bn1;
		end if
		
		stick_bit = 1;
		! STICK CHECK  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		if( cpdata(i)%state ) then
			
			Vt1	= UUd(3*bn1-2:3*bn1-1) + UUd(3*bn1)*r1_; 
			Vt2	= UUd(3*bn2-2:3*bn2-1) + UUd(3*bn2)*r2_; 

			Rt(1,:) = Tv;	Rt(2,:) = Nv;
			dvt = matmul(Rt,Vt1-Vt2);	! relative tangential velocity

			! 상대속도 방향이 변하면 무조건 stick으로 해준다.
			if( int(dvt(1)/dabs(dvt(1))) /= cpdata(i)%sign ) then
				
				! dabs(dvt(1)) > 1d-7 .and.

				! contact point가 edge point들의 사이에 있는지 조사한다.
				call check_projection_inside(UU, (/cpdata(i)%vn,cpdata(i)%en/), lv, cbit );
				
				if( cbit ) then				
					! edge block의 contact point를 update한다.
					if( cpdata(i)%bn(1) == cpdata(i)%en(1) ) then
						cpdata(i)%p(1,:) = lv;
					else
						cpdata(i)%p(2,:) = lv;
					end if

				else
					! contact update를 한다.
					!call Contact_Search(UU);
					Contact_search_bit = 1;

				end if

				cpdata(i)%state = 0;
				stick_bit = 0;
				new_contact = 1;

				!print *, "slip -> stick, i= ", i, " time=",tn
				! print *, "dxL=", dxL
			
			! 그렇지 않으면( SLIP 인 상태 ) contact update를 해준다.
			else
				
				! contact point가 edge point들의 사이에 있는지 조사한다.
				call check_projection_inside(UU, (/cpdata(i)%vn,cpdata(i)%en/), lv, cbit );
				
				if( cbit ) then	
				! edge block의 contact point를 update한다.
					if( cpdata(i)%bn(1) == cpdata(i)%en(1) ) then
						cpdata(i)%p(1,:) = lv;
					else
						cpdata(i)%p(2,:) = lv;
					end if
						
				else
					! edge to corner contact 일 경우에는 separation이 발생한다.
					if( BCST(bbn1,bbn2,1) == 1 ) then
						spN = spN + 1;
						spidx(spN) = i;
						BCST(bbn1,bbn2,:) = 0;		! ==> separation
					
					! edge to edge contact 인 경우에는 contact update를 한다.										
					else if( BCST(bbn1,bbn2,1) == 2 ) then
						!call Contact_Search(UU);
						Contact_search_bit = 1;

					end if

				end if
				
			end if
			
		end if 

		! SLIP CHECK  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 		if( 1 .and. cpdata(i)%state == 0 .and. stick_bit ) then	! stick_bit //.and. cpdata(i)%en(1) /= 1 
			
			if( dabs(FT) > mu*dabs(FN) .and. tn > 0.1d-6 ) then
 				cpdata(i)%sign = - int(FT/dabs(FT));
				cpdata(i)%state = 1;
				new_contact = 1;
				slip_bit = 1;
				!print *, "stick -> slip, i=", i, " time =", tn					
			end if		
			
		end if
		

		! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		! Edge to Corner contact 일 경우 매 step마다 새로운 contact이 있는지 조사한다. (for Fast Contact Search)
		if( 0 .and. BCST(bbn1,bbn2,1) == 1 ) then

			vbn = cpdata(i)%vn(1);		! vertex block number
			v1 = cpdata(i)%vn(2)-1;		! vertex number 1
			v2 = cpdata(i)%vn(2)+1;		! vertex number 2

			if( v1 == 0 ) then
				v1 = Block(vbn)%vN;
			end if

			if( v2 > Block(vbn)%vN ) then
				v2 = 1;
			end if 

			call Check_Vertex_inside_Edge( UU, cpdata(i)%en(1), cpdata(i)%en(2), vbn, v1, cbit1 );
			call Check_Vertex_inside_Edge( UU, cpdata(i)%en(1), cpdata(i)%en(2), vbn, v2, cbit2 );

			if( cbit1 .or. cbit2 ) then
				
				!vertex_in = (/bbn1, bbn2/);
				call Contact_Search( UU );
				!vertex_in = 0;
				Contact_search_bit = 1;

				print *,"t=",tn,"block=",vbn,"v1=",v1,"v2=",v2

			end if

		end if
		

		!	SEPARATION  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		if( 1 .and. dyL < 0 .and. cpN /=0 ) then
			
			print *, "##SEP. time =", tn, "Block, vertex=", cpdata(i)%vn(1), cpdata(i)%vn(2);
			
			! Edge to Corner contact 인 경우 
			if( BCST(bbn1,bbn2,1) == 1 ) then
			
				BCST(bbn1,bbn2,:) = 0;		! ==> 완전한 separation
	
			! Edge to Edge contact 인 경우
			else

				BCST(bbn1,bbn2,1) = 1;		! ==> edge to corner contact으로 바꿔줌..

				if( BCST(bbn1,bbn2,2) == i ) then	! cpdata index 번호를 바꿔줌.
					BCST(bbn1,bbn2,2) = BCST(bbn1,bbn2,3);
				end if
					
				BCST(bbn1,bbn2,3) = 0;

			end if
				
			! separation 이 일어난 index를 저장하고 separation 개수를 하나 증가 시킨다.
			spN = spN + 1;
			spidx(spN) = i; 

		end if	! end of SEPARATION
		!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		

	end do

	
	!*********************************************************************************************
	! separation 이 일어난 개수만큼 cpN 을 빼주고 
	! cpdata 및 BCST를 update 한다.
	if( spN /=0 ) then
		
		! 편의상 spidx의 끝 값을 cpN+1로 해준다.
		spidx(spN+1) = cpN + 1;
		
		! cpdata update
		j=1;	k=1;
		do i=1, cpN
			if( i /= spidx(j) ) then
				cpdata(k) = cpdata(i);
				k = k + 1;
			else
				j = j + 1;
			end if
		end do
		cpN = cpN - spN;
		
		! BCST update
		do i=1, TBN
			do j=i+1, TBN
				
				cbit1 = 1;	  cbit2 = 1;
				do k=1, spN
					
					if( cbit1 .and. BCST(i,j,2) > spidx(k) .and. BCST(i,j,2) < spidx(k+1) ) then
						BCST(i,j,2) = BCST(i,j,2) - k;
						cbit1 = 0;
					end if   

					if( cbit2 .and. BCST(i,j,3) > spidx(k) .and. BCST(i,j,3) < spidx(k+1) ) then
						BCST(i,j,3) = BCST(i,j,3) - k;
						cbit2 =0;
					end if

				end do   
		
			end do
		end do
	
		!!print *, "separation cpN =", cpN 
		new_contact = 1;		
				
	end if
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! 다음 step의 시간간격을 상대속도 데이타로부터 구한다.
	!c_vel = jfix( dlog10( c_vel ) ); 
	!factor = 4 + maxval( c_vel ) ;
	!if( factor <= 0 ) then		
	!	dt = org_dt;
	!else
	!	dt = org_dt; ! / factor; 
	!end if
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

! 한 점이 edge의 안쪽이 있는지 조사한다.
subroutine Check_Vertex_inside_Edge( UU, bn1, eN, bn2, vN, cbit ) 
	use global_data
	implicit none
	
	real*8 :: UU( 3*TBN );	! Input : Block Displacement at the center
	integer :: bn1, bn2		! Input : bn1 : edge block , bn2 : vertex block	
	integer :: eN, vN		! Input : edge N, vertex N 
	byte :: cbit			! Output : 내부에 포함 여부를 나타냄
			
	integer :: bi, ei;		
	real*8 :: p(2), q(2)
	real*8 :: p1(2), p2(2), a
	real*8 :: c, s, Rt1(2,2), Rt2(2,2);

	! begin vertex 와 end vertex number를 구한다.				
	bi = eN;	ei = bi+1;
	if( bi == Block(bn1)%vN ) then
		ei = 1;
	end if

	if( bi == 0 ) then
		bi = Block(bn1)%vN;
		ei = 1;
	end if

	c = dcos( UU(3*bn1) );		s = dsin( UU(3*bn1) );
	Rt1(1,1) = c;		Rt1(1,2) = -s;
	Rt1(2,1) = s;		Rt1(2,2) = c;

	c = dcos( UU(3*bn2) );		s = dsin( UU(3*bn2) );
	Rt2(1,1) = c;		Rt2(1,2) = -s;
	Rt2(2,1) = s;		Rt2(2,2) = c;
			
	p1 = matmul( Rt1,Block(bn1)%vtx(bi,:) ) + UU(3*bn1-2:3*bn1-1); 
	p2 = matmul( Rt1,Block(bn1)%vtx(ei,:) ) + UU(3*bn1-2:3*bn1-1);
			
	! edge vector		
	p = p2 - p1;	
	
	! vertex vector from p1
	q = matmul( Rt2,Block(bn2)%vtx(vN,:) ) + UU(3*bn2-2:3*bn2-1) - p1;  
		
	call cross_product(p,q,a);
	
	if( a > err ) then
		cbit = 1;
	else
		cbit = 0;
	end if

end subroutine

subroutine check_projection_inside(UU, VE_data, ppt, cbit )
	use global_data
	implicit none

	real*8 :: UU( 3*TBN );		! Input : Block Displacement at the center
	integer :: VE_data(4);	! Input : (vertex block N, vertex N, edge block N, edge N)
	real*8 :: ppt(2)	! Output : projection point vector (local frame) 
	byte :: cbit		! Output : 안쪽에 있으면 1, 바깥쪽에 있으면 0

	integer :: vbn, ebn	! vertex block number, edge block number
	integer :: bv, ev	! begin & end vertex number
	real*8 :: ept1(2), ept2(2), cpt(2)	! Edge points & coner point	
	real*8 :: mevector		! edge vector의 크기
	real*8 :: evector(2)	! unit edge vector
	real*8 :: pvector(2)	! location vector of coner point
	real*8 :: mmv
	real*8 :: Rt1(2,2), Rt2(2,2), c, s
		
	vbn = VE_data(1);		! vertex block number
	ebn = VE_data(3);		! edge block number

	c = dcos( UU(3*vbn) );		s = dsin( UU(3*vbn) );
	Rt1(1,1) = c;		Rt1(1,2) = -s;
	Rt1(2,1) = s;		Rt1(2,2) = c;

	cpt = matmul(Rt1, Block(vbn)%vtx(VE_data(2),:) ) + UU(3*vbn-2:3*vbn-1);

	c = dcos( UU(3*ebn) );		s = dsin( UU(3*ebn) );
	Rt2(1,1) = c;		Rt2(1,2) = -s;
	Rt2(2,1) = s;		Rt2(2,2) = c;
		
	
	bv = VE_data(4);	ev = bv + 1;

	if( bv == Block(ebn)%vN ) then
		ev = 1;
	end if

	ept1 = matmul(Rt2, Block(ebn)%vtx(bv,:) ) + UU(3*ebn-2:3*ebn-1); 
	ept2 = matmul(Rt2, Block(ebn)%vtx(ev,:) ) + UU(3*ebn-2:3*ebn-1);
	
	evector = ept2 - ept1;
	mevector = dsqrt( evector(1)**2 + evector(2)**2 );	! evector 의 크기
	evector = evector / mevector;

	! edge vector의 시작점으로 부터 coner point의 위치 벡터를 구한다.
	pvector = cpt - ept1;
	
	! edge vector에 정사영된 vector의 크기
	mmv = pvector(1)*evector(1) + pvector(2)*evector(2);
	if( mmv >= 0 .and. mmv <= mevector ) then
		cbit = 1;
	else 
		cbit = 0;
	end if

	if( cbit ) then
		! projection point의 좌표를 구한다. (global)
		ppt = mmv*evector + ept1;
		
		! 블록 중심에 대한 위치 벡터로 바꾼다.
		ppt = ppt - UU(3*ebn-2:3*ebn-1);

		! edge block의 역변환 행렬을 구성한다.		
		Rt2 = transpose(Rt2);
		!Rt2(1,1) = c;		Rt2(1,2) = s;
		!Rt2(2,1) = -s;		Rt2(2,2) = c;

		! body fixed frame에서의 위치 벡터를 구한다. (local)
		ppt = matmul( Rt2, ppt );

	end if


end subroutine
