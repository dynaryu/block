!/////////////////////////////////////////////////////////////////
module global_data
	implicit none
	save
	
	real*8 :: pi = 3.14159265358979d0;
	real*8 ::  g = 9.81d0;

	integer, parameter :: STN = 1  ! do 문 안에 있는 또 다른 do 문의 step N
		
	integer :: TBN			! Total Block Number		
	integer :: cpN			! Total contact point number
	integer, parameter :: MVN = 20	! Maximum Vertex Number

	character :: in_f*12, out_f*12
	
	type STATIC_BLOCK_DATA
		
		integer :: vN		! Total vertex number of block
		real*8 :: M			! Mass
		real*8 :: I			! Rotation Inertia at Mass Center
		real*8 :: Mc(2)		! Mass Center : xc, yc
		real*8 :: vtx(20,2)	! Vertex Point Data
	
	end type STATIC_BLOCK_DATA

	type CONTACT_DATA
		
		integer :: bn(2)	! Block Number : i, j
		integer :: c_type	! Contact Type => 1: point to point, 2: edge to point, 3: edge to edge
		integer :: en(2)	! Edge Number  
		integer :: vn(2)	! Vertex Number
		byte :: state		! State => 0:stick, 1:slide
		byte :: sign		! slip sign
		real*8  :: p(2,2)	! Position Vector from body center
		real*8 :: prev_a	! previous relative normal accelleration
		real*8 :: d_stable	! stable line depth
		byte :: stable		! stable bit
		byte :: p_stable	! previous stable  ## for test..........
		byte :: sep_bit		! separation bit ( if the value is 1, then virtual separation occurs )
		real*8 :: light_mass	! light mass between two contact block
		byte :: forced_contact	! forced contact bit
		
	end type CONTACT_DATA 

	type(STATIC_BLOCK_DATA), allocatable :: Block(:), nBlock(:)
	type(CONTACT_DATA), allocatable :: cpdata(:)

	integer, allocatable :: BCST(:,:,:)	! Block contact state data : edge to edge(2), edge to corner(1), slip..

	real*8 :: kp, xi, a1, kkk	! penalty spring & damping ratio, C=a1*K 
	real*8 :: mu			! friction coefficient
	real*8 :: dt, org_dt	! time step
	integer :: tstepN		! time step number
	
	real*8, allocatable :: t(:)
	real*8, allocatable :: U(:,:), Ud(:,:), U2d(:,:)	! DOF => disp, vel, acc	
	real*8, allocatable :: MM(:,:), Fc(:)
							! Mass, Central Force(contact force + body force)
	
	real*8, allocatable :: PO(:)	! Total Potential energy

	integer :: input_type	! ground motion input type : 1 => harmonic, 2 => earthquake
	real*8 :: Amp, frq		! harmonic motion : amplitude & frequency

	byte :: new_contact		! new contact occured ==> 1, otherwise ==> 0
	integer :: vertex_in(2)	! edge to corner contact 에서 새로운 contact이 발생했을 때 블록번호를 저장함..

	byte :: slip_bit		! Slip 이 발생했을 경우 다시 뒤로 가서 계산하기 위해
	byte :: Contact_search_bit	! contact search bit

	real*8 :: err = 1d-10;

	! for earthquake input motion
	real*8, allocatable :: dacc(:);
	real*8 :: adt;
	real*8 :: u_i, ud_i, prev_u_i, prev_ud_i;
	integer :: prev_idnum;
	real*8 :: pga_factor;

	! for Fast contact search data
	integer, allocatable :: FCS_data(:,:); ! (bn1, v1, v2, bn2, v1, v2, cross_bit)
	integer :: fcsN;	! Fast contact search No. (조사하고자 하는 edge to edge contact 의 No.)

	integer :: now_state ! bug 잡이용


		
end module global_data

!! MAIN CODE
PROGRAM MAIN			
	! use MSIMSLMD
	! USE PORTLIB
	use global_data;
	implicit none

	integer :: n, k, i, j, ip_bit, time1, time2, ssN, modN, aabit;
	real*8, allocatable :: X(:,:), Xd(:,:), X2d(:,:); 
	real*8, allocatable :: Uc(:,:), Ucd(:,:), tstn(:); 	
	real*8, allocatable :: KK(:,:), KKc(:,:), Kfg(:,:), RP(:), Fp(:), CC(:,:), CC2(:,:);	
	real*8 :: xg, xgd, xg2d
	real*8 :: pxg;
		
	real*8 :: VVg, VVs, KKs;		! potential energy, spring potential & Kinetic energy
	!real*8 :: DampF(3);
	real*8 :: aa1
	integer :: nnn ! 임시용임... 나중에 찾아서 지울것.

	! time1 = time();
	call Read_Data;		! get structure data
	kkk = kp;
	
	allocate( X( 3*(TBN-1),2 ), Xd( 3*(TBN-1),2 ), X2d( 3*TBN,2 ) );	
	allocate( Uc( 3*TBN, 0:STN), Ucd( 3*TBN, 0:STN) ); 
	allocate( KK( 3*(TBN-1), 3*(TBN-1) ), Kfg(3*(TBN-1),3), RP(3*(TBN-1)), Fp(3*(TBN-1)) ); 
	allocate( KKc( 3*(TBN-1), 3*(TBN-1) ) );
	allocate( CC( 3*(TBN-1), 3*(TBN-1) ), CC2( 3*(TBN-1), 3*(TBN-1) ) );
	allocate( tstn(0:STN) );	
	
	call Construct_Mass;
	call Initial_Condition_Setting;

	!##################################
	! INITIAL CONDITION FOR VELOCITY

	! single block slip
	if(0) then
		Ud(6,0) =  -1.0*1
		Ud(4,0) =   0.3 !0.200247499*-1+0.5; ! rocking
		Ud(5,0) =  -0.047997533; 
	end if
	
	! 2 block (r2b1)
	if(0) then
		Ud(4,0) =   0.150592490;
		Ud(5,0) =  -0.058497025;  
		Ud(6,0) =  -1.0; 

		Ud(7,0) =   0.401209504;
		Ud(8,0) =  -0.026970390; 
		Ud(9,0) =  -0.5;
	end if

	! 2 block (r2b2)
	if(0) then
		Ud(7,0) =   0.2;
		Ud(8,0) =  -0.05; 
		Ud(9,0) =  -1.0;
	end if
	
	!##################################	


	Uc(:,0) = U(:,0);	Ucd(:,0) = Ud(:,0);

	call getGround(0d0, xg, xgd, xg2d);
	X(:,1) = 0;
	Xd(:,1) = Ucd(4:3*(TBN-1),0);
	X2d(:,1) = 0;

	pxg = U(1,0);

	! TIME MARCHING	
	! ==============================================================		
	tstn(0) = 0d0;
	new_contact = 0;
	slip_bit = 0;
	Contact_search_bit = 1;
	aabit = 1;	
		
	ssN = int(tstepN / 100);
	

	do n=0, tstepN-1 
		
		do k=0, STN-1

			if( n == 0  ) then !.and. k == 0 ) then
				
				! Initial Contact Search
				call Contact_search( Uc(:,k) );
				!fcsN = 1;
				!FCS_data(1,:) = (/1,3,4,2,1,2,1/);
				!cpdata(1)%state = 1; !cpdata(2)%state = 1;
				!cpdata(1)%sign =  -1; !cpdata(2)%sign = -1; 
				
			end if

			!print *,"t=",tstn(k),"k=",k
			
100			continue;
			
				!cpdata(1)%state=1; cpdata(1)%sign = -1;
				!cpdata(2)%state=1; cpdata(2)%sign = -1;
				!cpdata(3)%state=1; cpdata(3)%sign = 1;
				!cpdata(4)%state=1; cpdata(4)%sign = -1;


			! Stiffness
			call Construct_Stiffness_Damping(Uc(:,k), KK, Kfg, KKc);	

			if( 0 .and. tstn(k) > 0.02) then
			
				out_f = "debug.dat";
				open (unit=30, file=out_f, action="write", access="append" )

				write (30,'(F8.5,14I2)') tstn(k), (Fcs_data(1,nnn), nnn=1,7), (Fcs_data(2,nnn), nnn=1,7);	
	
				close(30);

			end if


			! 새로운 contact이 있으면 damping 계수(a1)를 다시 계산한다.
			if( 1 .and. new_contact ) then

				call Get_a1(KK, CC2);
				
				new_contact = 0;

				cpN = cpN;
				cpdata = cpdata;
				!print *,"cpN", cpN;

			    ip_bit = 1; 
				i = 0;

				!print *, " Time =", tstn(k), "  a1 =", a1, " cpN=",cpN;	
				print *, " Time =", tstn(k), "  a1 =", a1, " cpN=",cpN;	
								
			end if
				
			CC = a1*KK;
			!DampF = matmul(CC,Ucd(4:6,k));

	 		! Get the ground disp., vel, acc.
			call getGround(tstn(k)+dt, xg, xgd, xg2d);

			call Get_Spring_Force( Uc(:,k), tstn(k), Fp, VVs ); ! it contains the body force
																								
			RP = FP - matmul(Kfg, (/xg-pxg, 0d0, 0d0/) ) - a1*matmul(Kfg, (/xgd, 0d0, 0d0/) );
			
			! Newmark Algorithm			
			call Newmark( KK, RP, CC, X(:,1), Xd(:,1), X2d(4:3*TBN,1), X(:,2), Xd(:,2), X2d(4:3*TBN,2) );

			Uc(4:3*TBN,k+1) = Uc(4:3*TBN,k) + X(:,2);
			Ucd(4:3*TBN,k+1) = Xd(:,2);

			Uc(1,k+1) = xg;		Uc(2,k+1) = Block(1)%Mc(2);   Uc(3,k+1) = 0;
			Ucd(1,k+1) = xgd; Ucd(2,k+1) = 0; Ucd(3,k+1) = 0;
			
			X2d(1,2) = xg2d; X2d(2,2) = 0; X2d(3,2) = 0;
			
		
			!!!!!!												
			if( 1 .and. aabit ) then
				call Check_Stick_Slip_Separation( tstn(k)+dt, Uc(:,k+1), Ucd(:,k+1), X2d(:,2) );

				if(new_contact) then
					call Contact_search( Uc(:,k+1) );
				end if
							
			end if

			aabit = 1;
				
			call Contact_search( Uc(:,k+1) );

			
			!========== My test ==================
			now_state = 1;
			if( tstn(k) > 0.9833) then
				now_state=0;
				print *, "time=", tstn(k)," BCST(7,8,1)=",bcst(7,8,1)
				if(bcst(7,8,1)==2) then
					print *,"data1=",cpdata(BCST(7,8,2))%vn," data2=",cpdata(BCST(7,8,3))%vn
				end if
			end if


			!========== My test ==================


			!==============================================
			!call Fast_contact_search( Uc(:,k+1) );

			if( 0 .and. Contact_search_bit ) then
				call Contact_search( Uc(:,k+1) );
				Contact_search_bit = 0;
			end if
			!==============================================
			
			if( slip_bit ) then
				slip_bit = 0;
				aabit = 0;
				goto 100;
			end if


			! re-initializing
			X(:,1) = 0;
			Xd(:,1) = Xd(:,2);
			X2d(:,1) = X2d(:,2);

			pxg = xg;

			tstn(k+1) = tstn(k) + dt;
			
		end do
		
		modN = mod( n, ssN )
		if( 1 .and. modN == 0 ) then
			print *," Time = ", tstn(k)
		end if

		! Save to U, Ud and t
		U(:,n+1) = Uc(:,k);
		Ud(:,n+1) = Ucd(:,k);
		t(n+1) = tstn(k);
				
		!###### Total Potential Energy ###############3
		if(0) then
		! spring potential
		call Get_Spring_Force( Uc(:,k), tstn(k), Fp, VVs );

		! block potential about gravity
		VVg = 0;
		do j=2,TBN
			VVg = VVg + Block(j)%M*9.81*(Uc(3*j-1,k)-U(3*j-1,0));
		end do

		! kinetic energy
		KKs = 0;
		do j=2,TBN
			KKs = KKs + 0.5*Block(j)%M*Ucd(3*j-2,k)**2 + 0.5*Block(j)%M*Ucd(3*j-1,k)**2 + 0.5*Block(j)%I*Ucd(3*j,k)**2;
		end do

		PO(n+1) = VVs + VVg + KKs;
		end if
		! #############################################
				
		! Re-initializing for next STN loop
		Uc = 0;		Ucd = 0;	 tstn = 0;
		Uc(:,0) = U(:,n+1);
		Ucd(:,0) = Ud(:,n+1);
		tstn(0) = t(n+1);
		
			
	end do
	! ==============================================================
		
	call Print_Result;

	! time2 = time();

	! print *," Analysis Time = ", time2-time1;

	
	deallocate(X); 
	deallocate(Xd);
	deallocate( X2d );
	deallocate( Uc, Ucd, tstn );	
	deallocate( KK, Kfg, RP, Fp, CC, CC2, KKc );

	call Var_deallocate;	


end program main	 


subroutine INITIAL_SHAPE(endtime,UI, UO)
	use global_data
	implicit none
		
	real*8 :: endtime
	real*8 :: UI( 3*TBN );
	real*8 :: UO( 3*TBN );
	real*8 :: X( 3*(TBN-1) ), Xd( 3*(TBN-1) ), X2d( 3*(TBN-1) )

	real*8 :: KK( 3*(TBN-1), 3*(TBN-1) ), Kfg(3*(TBN-1),3), RP(3*(TBN-1)), Fp(3*(TBN-1));
	real*8 :: KKc(3*(TBN-1), 3*(TBN-1) );  
	real*8 :: CC( 3*(TBN-1), 3*(TBN-1) );
	real*8 :: dummy


	integer :: k, tstepNp, j;
	

	X = 0;
	Xd = 0;
	X2d = 0;

	tstepNp = int(endtime/dt);

	do k=0, tstepNp

		call Contact_search( UI );
		
		! 새로운 contact이 있으면 damping 계수(a1)를 다시 계산한다.
		if( new_contact ) then
			
			call Construct_Stiffness_Damping(UI, KK, Kfg, KKc);

			cpN = cpN;
		
			do j=1, TBN-1
				if( j== 3 ) then
					KK(:,j) = 0;
					KK(j,:) = 0;				
				end if
			end do
			
			call Get_a1(KK, CC);
			new_contact = 0;
		
		end if
			
		CC = a1*KKc;
		
		call Get_Spring_Force( UI, .0, Fp, dummy ); ! it contains the body force
																								
		RP = FP; 
		
		do j=1, TBN-1
			RP(3*j) = 0;
			RP(3*j-2) = 0;
		end do
			
		! Newmark Algorithm			
		call Newmark( KK, RP, CC, X, Xd, X2d, X, Xd, X2d );
		
		UI(4:3*TBN) = UI(4:3*TBN) + X;
		X = 0;
			
	end do

	UO = UI;

end subroutine

! SUBROUTINE OF DATA READING
subroutine Read_Data
	use global_data
	implicit none

	integer :: i, j, aN
	real*8 :: endtime
	character :: ainput*12

	print *,"Data Input file>" 
	read *,in_f
		
	open (unit=3, file=in_f, action="read" )

	! READ TOTAL BLOCK NUMBER
	read (3,*)
	read (3,*) TBN;
	read (3,*)	! blank line

	call Initial_Processing; 
			
	! READ MASS & INERTIA DATA
	read (3,*) 
	do i=1,TBN
		read (3,*) Block(i)%M, Block(i)%I  
	end do
	read (3,*) ! blank line

	! READ MASS CENTER DATA
	read (3,*) 
	do i=1,TBN
		read (3,*) Block(i)%Mc(1), Block(i)%Mc(2)  
	end do
	read (3,*) ! blank line
			
	! READ VERTEX DATA
	read (3,*) 
	do i=1,TBN
		read (3,*) Block(i)%vN;
		do j=1,Block(i)%vN
			read(3,*) Block(i)%vtx(j,1), Block(i)%vtx(j,2)
		end do
		read (3,*) ! blank line	
	end do
		
	! READ TIME DATA
	read (3,*)	
	read (3,*) endtime, org_dt
	dt = org_dt;
	tstepN = int((endtime/dt+1)/STN);
	allocate( t( 0:tstepN ) );
	allocate( U(3*TBN,0:tstepN), Ud(3*TBN,0:tstepN), U2d(3*TBN,0:tstepN) );
	allocate( PO(0:tstepN) );
	PO = 0;

	read (3,*)
	
	! READ PENALTY DATA & DAMPING RATIO
	read (3,*)
	read (3,*) kp, xi
	read (3,*)
	
	! READ FRICTION COEFFICIENT
	read (3,*)
	read (3,*) mu
	read (3,*)

	! READ GROUND MOTION TYPE
	read (3,*)
	read (3,*) input_type
	
	!==> harmonic input motion 
	if( input_type == 1 ) then	
		read (3,*) Amp, frq
		Amp = Amp*g;
	end if

	!==> earthquake input motion
	if( input_type == 2 ) then
		read (3,*) adt;			! delta t of ground input acceleration
		read (3,*) u_i, ud_i;	! initial condtion
		read (3,*) ainput;		! earthquake input data file
		read (3,*) pga_factor;	! pga factor
		
		aN = 1 + int(endtime/adt);
		allocate( dacc(aN) );
				
		! acc. data read	
		open (unit=33, file=ainput, action="read" )
		
		do i=1,aN
			read (33,*) dacc(i);  
		end do
		
		dacc = pga_factor * dacc;
			
		close(33)		
						
	end if

	close(3);

	
			
end subroutine


! WRITE THE RESULT
subroutine Print_Result
	use global_data
	implicit none

	integer :: n, i

	out_f = "Disp_out.dat";
	open (unit=3, file=out_f, action="write" )

	out_f = "Vel_out.dat";
	open (unit=5, file=out_f, action="write" )

	!out_f = "Pot_out.dat";
	!open (unit=7, file=out_f, action="write" )

	print *," ";
	print *, "Now writing the result...";
	
	do n=0, tstepN
		write (3,'(100E14.7)')  t(n),  (U(i,n), i=1,3*TBN);	
		write (5,'(100E14.7)')  t(n),  (Ud(i,n), i=1,3*TBN);	
		!write (7,'(100E14.7)')  t(n),  PO(n);
	end do
	
	close(3);	close(5);

	n = tstepN;
	
!	print *, "t = ", t(n)					
!	print *, U(:,n)		

end subroutine

! INITIAL CONDTION SETTING
subroutine Initial_Condition_Setting
	use global_data
	implicit none

	integer :: i
	real*8 :: xg, xgd, xg2d

	U = 0;
	Ud = 0;
	
	!! Displacement		
	do i=1,TBN
		U(3*(i-1)+1:3*i-1,0) = Block(i)%Mc(1:2);
	end do

	!! Velocity
	call getGround(0d0, xg, xgd, xg2d)
	do i=1,TBN
		Ud(3*(i-1)+1,0) =  xgd;
	end do

	t(0) = 0d0;	

	vertex_in = 0;

end subroutine

! ININTIAL PROCESSING 
subroutine Initial_Processing
	use global_data
	implicit none

	cpN = 0;	
	allocate( Block(TBN), nBlock(TBN), BCST(TBN,TBN,3) );
	BCST = 0;
		
	allocate( cpdata(TBN*MVN) );
	allocate( FCS_data(TBN*2,7) );
	fcsN = 0;

end subroutine

					   

! VARIABLE DEALLOCATION
subroutine Var_deallocate
	
	use global_data
	implicit none

	deallocate( Block, nBlock, BCST, cpdata ); 
	deallocate( FCS_data );
	deallocate( t );
	deallocate( U, Ud, U2d );
	deallocate( MM, Fc  );
	deallocate( PO );

	if( input_type == 2 ) then
		deallocate( dacc );
	end if

end subroutine

