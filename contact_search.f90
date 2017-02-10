! 2 Dimensional cross product
subroutine cross_product(a, b, c)
	implicit none

	real*8 :: a(2), b(2), c

	c = a(1)*b(2) - a(2)*b(1);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_disp_projection(UU, VE_data, disp, ppt, cbit )
	use global_data
	implicit none

	real*8 :: UU( 3*TBN );		! Input : Block Displacement at the center
	integer :: VE_data(4);	! Input : (vertex block N, vertex N, edge block N, edge N)
	real*8 :: ppt(2)	! Output : projection point 
	real*8 :: disp		! Output : Disp. between coner point and projection point	
	byte :: cbit		! Output : 안쪽에 있으면 1, 바깥쪽에 있으면 0

	integer :: vbn, ebn	! vertex block number, edge block number
	integer :: bv, ev	! begin & end vertex number
	real*8 :: ept1(2), ept2(2), cpt(2)	! Edge points & coner point	
	real*8 :: mevector		! edge vector의 크기
	real*8 :: evector(2)	! unit edge vector
	real*8 :: pvector(2)	! location vector of coner point
	real*8 :: mmv
		
	vbn = VE_data(1);		! vertex block number
	ebn = VE_data(3);		! edge block number

	cpt = nBlock(vbn)%vtx(VE_data(2),:) + UU(3*vbn-2:3*vbn-1);
	
	bv = VE_data(4);	ev = bv + 1;

	if( bv == Block(VE_data(3))%vN ) then
		ev = 1;
	end if

	ept1 = nBlock(ebn)%vtx(bv,:) + UU(3*ebn-2:3*ebn-1); 
	ept2 = nBlock(ebn)%vtx(ev,:) + UU(3*ebn-2:3*ebn-1);
	
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
		
	! projection point의 좌표를 구한다.
	ppt = mmv*evector + ept1;

	! projection point 와 coner point와의 거리를 구한다.
	disp = sqrt( (ppt(1)-cpt(1))**2 + (ppt(2)-cpt(2))**2 );

end subroutine


! INITIAL CONTACT DETACTION SUBROUTINE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
subroutine Contact_Search( UU )

	use global_data
	implicit none
	
	real*8 :: UU(3*TBN);
	integer :: i,j, k, overlapped
	integer :: ivN(2,0:MVN)				! Inside vertex number(index)
	integer :: p_edge(2,0:MVN)
	real*8 ::  Rt(2,2), c, s;	! required for rotated vertex points
	real*8 :: OutBox(TBN,4), OvlBox(4)	! Outline Box of Each Block
		
	! 1st Step =================================================================
	!! Making the Outline Box of Each Block
	do i=1,TBN

		c = dcos( UU(3*i) );		s = dsin( UU(3*i) );
		Rt(1,1) = c;		Rt(1,2) = -s;
		Rt(2,1) = s;		Rt(2,2) = c;
		
		j = Block(i)%vN;
		do k=1,j
			nBlock(i)%vtx(k,:) = matmul(Rt,Block(i)%vtx(k,:));
		end do

		!vtxp = transpose( Block(i)%vtx );
		!vtxp = matmul( Rt, vtxp );

		!nBlock(i)%vtx = transpose( vtxp );   ! 새로 계산된 vertex data를 nBlock 에 저장한다.
		
		OutBox(i,1) = minval( nBlock(i)%vtx(1:j,1) ) + UU(3*i-2);	! X_min
		OutBox(i,2) = maxval( nBlock(i)%vtx(1:j,1) ) + UU(3*i-2);	! X_max
		OutBox(i,3) = minval( nBlock(i)%vtx(1:j,2) ) + UU(3*i-1);	! Y_min
		OutBox(i,4) = maxval( nBlock(i)%vtx(1:j,2) ) + UU(3*i-1);	! Y_max

		!OutBox(i,1) = minval( vtxp(1,1:j) ) + UU(3*i-2);	! X_min
		!OutBox(i,2) = maxval( vtxp(1,1:j) ) + UU(3*i-2);	! X_max
		!OutBox(i,3) = minval( vtxp(2,1:j) ) + UU(3*i-1);	! Y_min
		!OutBox(i,4) = maxval( vtxp(2,1:j) ) + UU(3*i-1);	! Y_max

	end do

	! Next Step ================================================================= 
	do i=1,TBN-1
		do j=i+1,TBN


			!! Checking the Overlapped Region
			call Check_Overlapped_Region( OutBox(i,:), OutBox(j,:), OvlBox, overlapped );

			!! Find the vertice which belong to the overlapped box
			if( overlapped ) then
				call Find_Vertex_Inside(UU, i, j, OvlBox, ivN)				
			end if

			!! Find the edges which pass the overlapped box
			if( ivN(1,1) /= 0 .or. ivN(2,1) /= 0 )	then
				call Find_Feasible_Edges(UU, i, j, ivN, OvlBox, p_edge)
			end if
			
			!! Find the crossing edges and return contact point data
			if( p_edge(1,0) /= 0 .and. p_edge(2,0) /= 0 )	then
				call Determine_ContactType(UU, i, j, p_edge)
				ivN = 0;	p_edge = 0;
			end if

		end do
	end do

end subroutine

! 2개의 edge가 서로 교차하는지 조사하고, 교차여부와 cross product data를 돌려준다.
subroutine Get_Crossing_Data(UU, bn1, bn2, eN, cbit, cv)
	use global_data
	implicit none
	
	real*8 :: UU( 3*TBN );	! Input : Block Displacement at the center
	integer :: bn1, bn2		! Input : Block number	
	integer :: eN(2)		! Input : edge number data 
	byte :: cbit			! Output : 교차여부를 나타냄
	real*8 :: cv(4)			! Output : cross product 결과를 저장함. (과거의 aa, bb, cc, dd)
			
	integer :: bi, ei, bj, ej;		
	real*8 :: p(2), q(2), a(2), b(2), c(2), d(2)
	real*8 :: p1(2), p2(2), q1(2), q2(2)

	! begin vertex 와 end vertex number를 구한다.				
	bi = eN(1);	ei = bi+1;
	if( bi == Block(bn1)%vN ) then
		ei = 1;
	end if
			
	bj = eN(2);	ej = bj+1;			
	if( bj == Block(bn2)%vN ) then
		ej = 1;
	end if

	! block i와 j의 edge를 cross product를 통해 교점을 조사한다.
	p1 = nBlock(bn1)%vtx(bi,:) + UU(3*bn1-2:3*bn1-1); 
	p2 = nBlock(bn1)%vtx(ei,:) + UU(3*bn1-2:3*bn1-1);
	q1 = nBlock(bn2)%vtx(bj,:) + UU(3*bn2-2:3*bn2-1); 
	q2 = nBlock(bn2)%vtx(ej,:) + UU(3*bn2-2:3*bn2-1);
			
	p = p2 - p1;		q = q2 - q1;

	a = q1 - p1;		b = q2 - p1;
	c = p1 - q1;		d = p2 - q1;

	call cross_product(p,a,cv(1));	call cross_product(p,b,cv(2));
	call cross_product(q,c,cv(3));	call cross_product(q,d,cv(4));

	! control the numerical error
	if( dabs(cv(1)) < err ) then;	cv(1) = 0; end if
	if( dabs(cv(2)) < err ) then;	cv(2) = 0; end if
	if( dabs(cv(3)) < err ) then;	cv(3) = 0; end if
	if( dabs(cv(4)) < err ) then;	cv(4) = 0; end if
	
	cbit = 0;
	if( cv(1)*cv(2) <= 0 .and. cv(3)*cv(4) <= 0 ) then
		if( cv(1)*cv(2) < 0 .or. cv(3)*cv(4) < 0 ) then
			cbit = 1;
		end if
	end if

end subroutine Get_Crossing_Data

! Get the distance between vertices
subroutine get_disp_vertices(UU, bn1, vn1, bn2, vn2, disp )
	use global_data
	implicit none
	
	real*8 :: UU( 3*TBN );		! Input : Block Displacement at the center	
	integer :: bn1, vn1, bn2, vn2	! input : block N, vertex N, ...
	real*8 :: disp					! output : displacement(distance)
	real*8 :: X1(2), X2(2)
		
	X1 = nBlock(bn1)%vtx(vn1,:) + UU(3*bn1-2:3*bn1-1); 
	X2 = nBlock(bn2)%vtx(vn2,:) + UU(3*bn2-2:3*bn2-1); 

	disp = dsqrt( (X1(1)-X2(1))**2 + (X1(2)-X2(2))**2 );
		
end subroutine

! Determine light mass
subroutine Light_Mass(M1, M2, lmass)
	implicit none
	
	real*8 :: M1, M2, lmass
	
	if(M1 > M2) then
		lmass = M2;
	else
		lmass = M1;
	end if

	if( M1 == 0 ) then
		lmass = M2;
	end if

	if( M2 == 0 ) then
		lmass = M1;
	end if

end subroutine

! Determine the contact type
subroutine Determine_ContactType(UU, bn1, bn2, edge)
	use global_data
	implicit none
	
	real*8 :: UU( 3*TBN );		! Input : Block Displacement at the center
	integer :: bn1, bn2			! Input : Block number	
	integer :: edge(2,0:MVN)	! Input : edge data : en1, en2, ...
	
	byte :: cbit, s_bit(2)
	integer :: bn(2), temp_bn;
	integer :: edgeN(4,4), cedgeN(4,4), eNN
	integer :: ct_data(2,6), temp(6);	! (/en1, en2, bn, normal_en, vn1, vn2/,...)
	integer :: diff, hiN;
	integer :: i, j, bi, ei, bj, ej;	
	real*8 :: cv(4)
	real*8 :: disp1, disp2, ppt1(2), ppt2(2);	
	real*8 :: theta, c, s, Rt1(2,2), Rt2(2,2);
	integer :: ve_data1(4), ve_data2(4), c_vnd(2), t_vnd(2);

	type(CONTACT_DATA) :: temp1_cpdata, temp2_cpdata
	
	
	bn = (/bn1, bn2/);
	eNN = 0; edgeN = 0;
	do i=1, edge(1,0)
			
		! begin vertex 와 end vertex number를 구한다.				
		bi = edge(1,i);	ei = bi+1;
		if( bi == Block(bn1)%vN ) then
			ei = 1;
		end if
			
		do j=1,edge(2,0)
			
			bj = edge(2,j);	ej = bj+1;			
			if( bj == Block(bn2)%vN ) then
				ej = 1;
			end if
				
			! 교점이 존재하는지 조사한다. 
			call Get_Crossing_Data(UU, bn1, bn2, (/bi, bj/), cbit, cv)

				
			! 교점이 존재할 경우 각 edge 번호를 edgeN에 저장한다.				
			if( cbit ) then
			
				eNN = eNN + 1;
				if( eNN > 4 ) then
					!!print *, " Index overflow ( Determine_ContactType, edgeN ) "
				
				else  
					edgeN(1,eNN) = bi;		! edge number 저장
					edgeN(2,eNN) = bj;		! edge number 저장
					edgeN(3,eNN) = i;		! edge index 저장
					edgeN(4,eNN) = j;		! edge index 저장
				end if
			
			end if

		end do ! j
	end do ! i
	
	if( eNN == 0 ) then
		goto 111;
	end if

	cedgeN = edgeN;
	! 만약 교차하는 edge의 집합이 3개 이상일 때는 같은 번호의 edge를 찾아서 2개만 선택한다. 
	if( eNN > 2 ) then
		do i=1, eNN-1
			if( edgeN(1,i) == edgeN(1,i+1) .or. edgeN(2,i) == edgeN(2,i+1) ) then
				edgeN(:,1) = cedgeN(:,i);
				edgeN(:,2) = cedgeN(:,i+1);  
			end if
		end do
	end if

	! edge numbers
	ct_data(:,1:2) = edgeN(1:2,1:2);
	! block numbers
	ct_data(:,3) = (/bn1, bn2/);
	! vertex number를 모두 0으로 초기화 한다.
	ct_data(:,5:6) = 0;
	
	!! 기본 database 구축 ( ct_data )
	do i=1,2
		
		diff = abs(edgeN(2+i,2) - edgeN(2+i,1));		
		
		if( diff == 0 ) then
			ct_data(i,4) = edgeN(i,1);	! pair block의 normal or tangential vector edge number				
		end if

		if( diff == 1 ) then
			ct_data(i,4) = 0;
			! 높은 index의 edge번호를 찾는다.
			if( edgeN(2+i,1) > edgeN(2+i,2) ) then
				hiN = 1;
			else 
				hiN = 2;
			end if

			ct_data(i,5) = edgeN(i,hiN);
		end if

		if( diff == 2 ) then
			hiN = (edgeN(2+i,1) + edgeN(2+i,2))/2;	! edge number index
			ct_data(i,4) = edge(i,hiN);		! pair block의 normal or tangential vector edge number 
			ct_data(i,5) = edge(i,hiN);
			ct_data(i,6) = edge(i,hiN) + 1;
			if( ct_data(i,6) > Block(bn(i))%vN ) then
				ct_data(i,6) = 1;
			end if
		end if

	end do

	!! Determine Contact Type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if( ct_data(1,4) == 0 ) then
	!! 만약 edge가 없으면 즉, 0 일경우 0이 되는 것을 뒤로 보낸다.
		temp = ct_data(1,:);
		ct_data(1,:) = ct_data(2,:);
		ct_data(2,:) = temp;
	
	else if( ct_data(2,5) == 0 ) then
	!! 만약 edge가 있고, vertex가 없으면 0이 되는 것을 앞으로 보낸다.
		temp = ct_data(2,:);
		ct_data(2,:) = ct_data(1,:);
		ct_data(1,:) = temp;
	
	end if
	
	! 편의를 위해 bn1과 bn2를 다시 설정함..
	bn1 = ct_data(1,3);
	bn2 = ct_data(2,3);

	! 역변환 행렬을 만듦  => local frame 에서의 위치벡터를 구하기 위해서...
	theta = -UU(3*bn1);		
		
	c = dcos( theta );		s = dsin( theta );
	Rt1(1,1) = c;			Rt1(1,2) = -s;
	Rt1(2,1) = s;			Rt1(2,2) = c;  
		
	theta = -UU(3*bn2);		
		
	c = dcos( theta );		s = dsin( theta );
	Rt2(1,1) = c;			Rt2(1,2) = -s;
	Rt2(2,1) = s;			Rt2(2,2) = c;  


	!! case 별로 조사를 함...
	! ###########################################################################################
	if( ct_data(1,4) == 0 ) then
		
		! 둘다 0 일 때.. ==> Edge to Edge contact
		if( ct_data(2,4) == 0 ) then
			
			! i block 		
			!===========================================================================
			temp1_cpdata%bn = (/bn1, bn2/);

			! i block의 vertex에 대하여...
			!! vertex block N, vertex N, edge block N, edge N ( 1st edge of j block )
			ve_data1 = (/bn1, ct_data(1,5), bn2, ct_data(2,1)/);
			call get_disp_projection(UU, ve_data1, disp1, ppt1, cbit );

			!! vertex block N, vertex N, edge block N, edge N ( 2nd edge of j block )
			ve_data2 = (/bn1, ct_data(1,5), bn2, ct_data(2,2)/);
			call get_disp_projection(UU, ve_data2, disp2, ppt2, cbit );

			if( disp1 > disp2 ) then
				ve_data1 = ve_data2;
				ppt1 = ppt2;				
			end if

			!! block 중심에 대한 위치 벡터로 바꿔준다. (global frame)
			temp1_cpdata%p(1,:) =  nBlock(bn1)%vtx(ve_data1(2),:);
			temp1_cpdata%p(2,:) =  ppt1 - UU(3*bn2-2:3*bn2-1);

			!! Local frame 의 위치 벡터로 바꿔준다. 
			temp1_cpdata%p(1,:) = matmul(Rt1, temp1_cpdata%p(1,:));  
			temp1_cpdata%p(2,:) = matmul(Rt2, temp1_cpdata%p(2,:));

			!! vertex number
			temp1_cpdata%vn = ve_data1(1:2);

			!! normal & tangential vector edge
			temp1_cpdata%en = ve_data1(3:4);	 ! (/edge block, edge number/);
			
			! j block			
			!===========================================================================
			temp2_cpdata%bn = temp1_cpdata%bn;
			
			! j block의 vertex에 대하여...
			!! vertex block N, vertex N, edge block N, edge N ( 1st edge of j block )
			ve_data1 = (/bn2, ct_data(2,5), bn1, ct_data(1,1)/);
			call get_disp_projection(UU, ve_data1, disp1, ppt1, cbit );

			!! vertex block N, vertex N, edge block N, edge N ( 2nd edge of j block )
			ve_data2 = (/bn2, ct_data(2,5), bn1, ct_data(1,2)/);
			call get_disp_projection(UU, ve_data2, disp2, ppt2, cbit );

			if( disp1 > disp2 ) then
				ve_data1 = ve_data2;
				ppt1 = ppt2;				
			end if

			!! block 중심에 대한 위치 벡터로 바꿔준다. (global frame)
			temp2_cpdata%p(1,:) = ppt1 - UU(3*bn1-2:3*bn1-1);			
			temp2_cpdata%p(2,:) = nBlock(bn2)%vtx(ve_data1(2),:);
			
			!! Local frame 의 위치 벡터로 바꿔준다. 
			temp2_cpdata%p(1,:) = matmul(Rt1, temp2_cpdata%p(1,:));  
			temp2_cpdata%p(2,:) = matmul(Rt2, temp2_cpdata%p(2,:));

			!! vertex number
			temp2_cpdata%vn = ve_data1(1:2);

			!! normal & tangential vector edge
			temp2_cpdata%en = ve_data1(3:4);	 ! (/edge block, edge number/);
			
			!! 얻어진 contact data 가 새로운 contact 인지 조사한다. ------------------------
			if( bn1 > bn2 ) then;	! 블록 번호를 작은걸 앞에...
				temp_bn = bn1;
				bn1 = bn2; 
				bn2 = temp_bn;
			end if
			
			!!!		기존에 아무런 contact data가 없을 경우 
			if( BCST(bn1,bn2,1) == 0 ) then
				BCST(bn1,bn2,1) = 2;		! edge to edge contact
				BCST(bn1,bn2,2) = cpN + 1;	! cpdata index number
				BCST(bn1,bn2,3) = cpN + 2;	! cpdata index number

				cpdata(cpN+1) = temp1_cpdata;
				cpdata(cpN+1)%state = 0;	! 초기에는 무조건 stick 으로 해준다.
				cpdata(cpN+2) = temp2_cpdata;
				cpdata(cpN+2)%state = 0;	! 초기에는 무조건 stick 으로 해준다.

				cpdata(cpN+1)%prev_a = -99999;	! 초기 상대 가속도 (첫번째 값임을 의미하기 위해)
				cpdata(cpN+2)%prev_a = -99999;	! 초기 상대 가속도 

				cpdata(cpN+1)%sep_bit = 0;
				cpdata(cpN+2)%sep_bit = 0;

				cpdata(cpN+1)%d_stable = 0.0;
				cpdata(cpN+2)%d_stable = 0.0;

				cpdata(cpN+1)%stable = 0;
				cpdata(cpN+2)%stable = 0;

				call Light_Mass(Block(cpdata(cpN+1)%bn(1))%M, Block(cpdata(cpN+1)%bn(2))%M, cpdata(cpN+1)%light_mass );
				call Light_Mass(Block(cpdata(cpN+2)%bn(1))%M, Block(cpdata(cpN+2)%bn(2))%M, cpdata(cpN+2)%light_mass );

				cpN = cpN + 2;								
				new_contact = 1;



				!!print *, " Edge to Edge contact => type 1, cpN =", cpN
			
			!!!		기존에 edge to corner contact 인 경우 
			else if( BCST(bn1,bn2,1) == 1 ) then
				
				! vertex가 관입했을 때만 조사를 한다.
				if( 1 ) then ! vertex_in(1) == bn1 .and. vertex_in(2) == bn2 ) then

					BCST(bn1,bn2,1) = 2;		! edge to edge contact 으로 바꿔줌
					BCST(bn1,bn2,3) = cpN + 1;	! cpdata index number (2nd index)

					! 거리를 바탕으로 어떤 것이 기존에 있는 contact point(vertex)인지 찾는다..
					c_vnd = cpdata( BCST(bn1,bn2,2) )%vn;
				
					t_vnd = temp1_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp1 );
				
					t_vnd = temp2_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp2 );
				
					if( disp1 < disp2  ) then
						cpdata(cpN+1) = temp2_cpdata;
						cpdata(cpN+1)%state = 0;	! 초기에는 무조건 stick 으로 해준다.
					else					
						cpdata(cpN+1) = temp1_cpdata;
						cpdata(cpN+1)%state = 0;	! 초기에는 무조건 stick 으로 해준다.
					end if
					
					cpdata(cpN+1)%prev_a = -99999;	
					cpdata(cpN+1)%sep_bit = 0;		
					
					cpdata(cpN+1)%d_stable = 0.0;
					cpdata(cpN+1)%stable = 0;

					call Light_Mass(Block(cpdata(cpN+1)%bn(1))%M, Block(cpdata(cpN+1)%bn(2))%M, cpdata(cpN+1)%light_mass );
								
					cpN = cpN + 1;
								
					new_contact = 1;

					print *, " Edge to Edge contact => type 2, cpN =", cpN, " Block=", cpdata(cpN)%vn

				end if
				
			!!!		기존에 edge to edge contact 		
			else if( BCST(bn1,bn2,1) == 2 ) then
				
				s_bit(1) = cpdata( BCST(bn1,bn2,2) )%state;
				s_bit(2) = cpdata( BCST(bn1,bn2,3) )%state;

				!!!! 둘 다 slip 일 경우	
				if( s_bit(1) .and. s_bit(2) ) then
				
					cpdata( BCST(bn1,bn2,2) )%bn = temp1_cpdata%bn;
					cpdata( BCST(bn1,bn2,2) )%p = temp1_cpdata%p; 
					cpdata( BCST(bn1,bn2,2) )%vn = temp1_cpdata%vn;
					cpdata( BCST(bn1,bn2,2) )%en = temp1_cpdata%en;
					
					cpdata( BCST(bn1,bn2,3) )%bn = temp2_cpdata%bn;
					cpdata( BCST(bn1,bn2,3) )%p = temp2_cpdata%p;
					cpdata( BCST(bn1,bn2,3) )%vn = temp2_cpdata%vn;
					cpdata( BCST(bn1,bn2,3) )%en = temp2_cpdata%en;
					
					!temp1_cpdata%state = 1;	temp2_cpdata%state = 1;
					!temp1_cpdata%sign = cpdata( BCST(bn1,bn2,2) )%sign;
					!temp2_cpdata%sign = cpdata( BCST(bn1,bn2,3) )%sign;
					
					!cpdata( BCST(bn1,bn2,2) ) = temp1_cpdata;					
					!cpdata( BCST(bn1,bn2,3) ) = temp2_cpdata;
				
				!!!! 첫번째 것만 slip 일 경우
				else if( s_bit(1) .and. s_bit(2)==0 ) then
					
					c_vnd = cpdata( BCST(bn1,bn2,2) )%vn;
				
					t_vnd = temp1_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp1 );
				
					t_vnd = temp2_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp2 );
				
					if( disp1 < disp2  ) then						
						cpdata( BCST(bn1,bn2,2) )%bn = temp1_cpdata%bn;
						cpdata( BCST(bn1,bn2,2) )%p = temp1_cpdata%p;
						cpdata( BCST(bn1,bn2,2) )%vn = temp1_cpdata%vn;
						cpdata( BCST(bn1,bn2,2) )%en = temp1_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,2) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,2) )%sign;
						!cpdata( BCST(bn1,bn2,2) ) = temp1_cpdata;
						!cpdata( BCST(bn1,bn2,2) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,2) )%sign = s_bit(2); 
					else					
						cpdata( BCST(bn1,bn2,2) )%bn = temp2_cpdata%bn;
						cpdata( BCST(bn1,bn2,2) )%p = temp2_cpdata%p;
						cpdata( BCST(bn1,bn2,2) )%vn = temp2_cpdata%vn;
						cpdata( BCST(bn1,bn2,2) )%en = temp2_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,2) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,2) )%sign;
						!cpdata( BCST(bn1,bn2,2) ) = temp2_cpdata;
						!cpdata( BCST(bn1,bn2,2) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,2) )%sign = s_bit(2); 
					end if

				!!!! 두번째 것만 slip 일 경우
				else if( s_bit(1)==0 .and. s_bit(2) ) then
					
					c_vnd = cpdata( BCST(bn1,bn2,3) )%vn;
				
					t_vnd = temp1_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp1 );
				
					t_vnd = temp2_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp2 );
				
					if( disp1 < disp2  ) then
						cpdata( BCST(bn1,bn2,3) )%bn = temp1_cpdata%bn;
						cpdata( BCST(bn1,bn2,3) )%p = temp1_cpdata%p;
						cpdata( BCST(bn1,bn2,3) )%vn = temp1_cpdata%vn;
						cpdata( BCST(bn1,bn2,3) )%en = temp1_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,3) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,3) )%sign;
						!cpdata( BCST(bn1,bn2,3) ) = temp1_cpdata;
						!cpdata( BCST(bn1,bn2,3) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,3) )%sign = s_bit(2); 
					else					
						cpdata( BCST(bn1,bn2,3) )%bn = temp2_cpdata%bn;
						cpdata( BCST(bn1,bn2,3) )%p = temp2_cpdata%p;
						cpdata( BCST(bn1,bn2,3) )%vn = temp2_cpdata%vn;
						cpdata( BCST(bn1,bn2,3) )%en = temp2_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,3) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,3) )%sign;
						!cpdata( BCST(bn1,bn2,3) ) = temp2_cpdata;
						!cpdata( BCST(bn1,bn2,3) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,3) )%sign = s_bit(2); 
					end if
				
				end if			
			
			end if			
		
		end if ! ( 2번째 edge 번호도 0 인것의 끝 )
		
	! ###########################################################################################
	! 첫번째 edge 번호가 0 이 아닐 경우 
	else
	
		! Edge to Corner Contact ************************************************************ 
		if( ct_data(1,5) == 0 .and. ct_data(2,6) == 0 ) then
			
			! block data
			temp1_cpdata%bn = (/bn1, bn2/);
			
			! vertex data
			temp1_cpdata%vn(1) = bn2;
			temp1_cpdata%vn(2) = ct_data(2,5);

			! edge data
			temp1_cpdata%en(1) = bn1;
			temp1_cpdata%en(2) = ct_data(1,4);

			! vertex block N, vertex N, edge block N, edge N 
			ve_data1 = (/bn2, ct_data(2,5), bn1, ct_data(1,4)/);
			call get_disp_projection(UU, ve_data1, disp1, ppt1, cbit );

			! block 중심에 대한 위치 벡터로 바꿔준다. (global frame)
			temp1_cpdata%p(1,:) = ppt1 - UU(3*bn1-2:3*bn1-1);			
			temp1_cpdata%p(2,:) = nBlock(bn2)%vtx(ve_data1(2),:);
			
			! Local frame 의 위치 벡터로 바꿔준다. 
			temp1_cpdata%p(1,:) = matmul(Rt1, temp1_cpdata%p(1,:));  
			temp1_cpdata%p(2,:) = matmul(Rt2, temp1_cpdata%p(2,:));

			! 얻어진 contact data 가 새로운 contact 인지 조사한다. ------------------------
			if( bn1 > bn2 ) then;	! 블록 번호를 작은걸 앞에... (편의상)
				temp_bn = bn1;
				bn1 = bn2; 
				bn2 = temp_bn;
			end if

			! 기존에 아무런 contact data가 없을 경우
			if( BCST(bn1,bn2,1) == 0 ) then
				
				BCST(bn1,bn2,1) = 1
				BCST(bn1,bn2,2) = cpN + 1
				
				cpdata(cpN+1) = temp1_cpdata;
				cpdata(cpN+1)%state = 0;
				cpdata(cpN+1)%prev_a = -99999;
				cpdata(cpN+1)%sep_bit = 0;
				cpdata(cpN+1)%d_stable = 0.0;
				cpdata(cpN+1)%stable = 0;
				
				call Light_Mass(Block(cpdata(cpN+1)%bn(1))%M, Block(cpdata(cpN+1)%bn(2))%M, cpdata(cpN+1)%light_mass );
								
				cpN = cpN + 1;				
				new_contact = 1;
				
				!!print *, " Edge to Corner contact => type 3, cpN =", cpN
				!!print *, " cpN =", cpN
			
			! 기존에 Edge to Corner Contact인 경우
			else if( BCST(bn1,bn2,1) == 1 ) then
				! vertex data 와 edge data만 update 해준다.
				!cpdata( BCST(bn1,bn2,2) )%vn = temp1_cpdata%vn;
				!cpdata( BCST(bn1,bn2,2) )%en = temp1_cpdata%en;
				! dddd
				
			! 기존에 Edge to Edge Contact인 경우
			else if( BCST(bn1,bn2,1) == 2 ) then
									
				!BCST(bn1,bn2,1) = 1
				!BCST(bn1,bn2,2) = cpN + 1
				! 어떤 것이 남아 있는 vertex인지 검색한다.
				! ##$$$ 이부분은 미흡한 점이 있으니 나중에 다시 보완한다.... $$$##
				!if( cpdata( BCST(bn1,bn2,2) )%vn(1) = temp1_cpdata%vn(1) .and. cpdata( BCST(bn1,bn2,2) )%vn(2) = temp1_cpdata%vn(2) ) then
			
			end if
			
		
		! Edge to Edge contact **************************************************************
		else if( ct_data(1,5) == 0 .and. ct_data(2,6) /= 0 ) then	

			! block data
			temp1_cpdata%bn = (/bn1, bn2/);
			temp2_cpdata%bn = (/bn1, bn2/);
			
			! vertex data
			temp1_cpdata%vn(1) = bn2; 	temp1_cpdata%vn(2) = ct_data(2,5);
			temp2_cpdata%vn(1) = bn2; 	temp2_cpdata%vn(2) = ct_data(2,6);

			! edge data
			temp1_cpdata%en(1) = bn1;	temp1_cpdata%en(2) = ct_data(1,4);
			temp2_cpdata%en(1) = bn1;	temp2_cpdata%en(2) = ct_data(1,4);

			! vertex block N, vertex N, edge block N, edge N 
			ve_data1 = (/bn2, ct_data(2,5), bn1, ct_data(1,4)/);
			call get_disp_projection(UU, ve_data1, disp1, ppt1, cbit );

			! vertex block N, vertex N, edge block N, edge N 
			ve_data2 = (/bn2, ct_data(2,6), bn1, ct_data(1,4)/);
			call get_disp_projection(UU, ve_data2, disp2, ppt2, cbit );

			! block 중심에 대한 위치 벡터로 바꿔준다. (global frame)
			temp1_cpdata%p(1,:) = ppt1 - UU(3*bn1-2:3*bn1-1);			
			temp1_cpdata%p(2,:) = nBlock(bn2)%vtx(ve_data1(2),:);

			temp2_cpdata%p(1,:) = ppt2 - UU(3*bn1-2:3*bn1-1);			
			temp2_cpdata%p(2,:) = nBlock(bn2)%vtx(ve_data2(2),:);
			
			! Local frame 의 위치 벡터로 바꿔준다. 
			temp1_cpdata%p(1,:) = matmul(Rt1, temp1_cpdata%p(1,:));  
			temp1_cpdata%p(2,:) = matmul(Rt2, temp1_cpdata%p(2,:));

			temp2_cpdata%p(1,:) = matmul(Rt1, temp2_cpdata%p(1,:));  
			temp2_cpdata%p(2,:) = matmul(Rt2, temp2_cpdata%p(2,:));

			! 얻어진 contact data 가 새로운 contact 인지 조사한다. ------------------------
			if( bn1 > bn2 ) then;	! 블록 번호를 작은걸 앞에... (편의상)
				temp_bn = bn1;
				bn1 = bn2; 
				bn2 = temp_bn;
			end if

			! 기존에 아무런 contact data가 없을 경우
			if( BCST(bn1,bn2,1) == 0 ) then
				
				BCST(bn1,bn2,1) = 2;
				BCST(bn1,bn2,2) = cpN + 1; 
				BCST(bn1,bn2,3) = cpN + 2;

				cpdata(cpN+1) = temp1_cpdata;
				cpdata(cpN+2) = temp2_cpdata;
				cpdata(cpN+1)%state = 0;
				cpdata(cpN+2)%state = 0;
				
				cpdata(cpN+1)%prev_a = -99999;
				cpdata(cpN+2)%prev_a = -99999;

				cpdata(cpN+1)%sep_bit = 0;
				cpdata(cpN+2)%sep_bit = 0;

				cpdata(cpN+1)%d_stable = 0.0;
				cpdata(cpN+2)%d_stable = 0.0;

				cpdata(cpN+1)%stable = 0;
				cpdata(cpN+2)%stable = 0;

				call Light_Mass(Block(cpdata(cpN+1)%bn(1))%M, Block(cpdata(cpN+1)%bn(2))%M, cpdata(cpN+1)%light_mass );
				call Light_Mass(Block(cpdata(cpN+2)%bn(1))%M, Block(cpdata(cpN+2)%bn(2))%M, cpdata(cpN+2)%light_mass );

				cpN = cpN + 2;
				new_contact = 1;

				!!print *, "disp1 =", disp1, " disp2=",disp2
				!!print *, " Edge to Edge contact => type 4, cpN =", cpN				

			
			! 기존에 Edge to Corner Contact인 경우
			else if( BCST(bn1,bn2,1) == 1 ) then

				! vertex가 관입했을 때만 조사를 한다.
				if( 1 ) then ! vertex_in(1) == bn1 .and. vertex_in(2) == bn2 ) then

					BCST(bn1,bn2,1) = 2;		! edge to edge contact 으로 바꿔줌
					BCST(bn1,bn2,3) = cpN + 1;	! cpdata index number (2nd index)

					! 거리를 바탕으로 어떤 것이 기존에 있는 contact point(vertex)인지 찾는다..
					c_vnd = cpdata( BCST(bn1,bn2,2) )%vn;
				
					t_vnd = temp1_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp1 );
				
					t_vnd = temp2_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp2 );
				
					if( disp1 < disp2  ) then
						cpdata(cpN+1) = temp2_cpdata;
						cpdata(cpN+1)%state = 0;	! 초기에는 무조건 stick 으로 해준다.
					else					
						cpdata(cpN+1) = temp1_cpdata;
						cpdata(cpN+1)%state = 0;	! 초기에는 무조건 stick 으로 해준다.
					end if
					
					cpdata(cpN+1)%prev_a = -99999;
					cpdata(cpN+1)%sep_bit = 0;
					cpdata(cpN+1)%d_stable = 0.0;
					cpdata(cpN+1)%stable = 0;

					call Light_Mass(Block(cpdata(cpN+1)%bn(1))%M, Block(cpdata(cpN+1)%bn(2))%M, cpdata(cpN+1)%light_mass );
																									
					cpN = cpN + 1;
					new_contact = 1;

					!!print *, " Edge to Edge contact => type 5, cpN =", cpN

				end if


			! 기존에 Edge to Edge Contact인 경우
			else if( BCST(bn1,bn2,1) == 2 ) then
				
				s_bit(1) = cpdata( BCST(bn1,bn2,2) )%state;
				s_bit(2) = cpdata( BCST(bn1,bn2,3) )%state;

				!!!! 둘 다 slip 일 경우	
				if( s_bit(1) .and. s_bit(2) ) then
				
					!temp1_cpdata%state = 1;	temp2_cpdata%state = 1;
					!temp1_cpdata%sign = cpdata( BCST(bn1,bn2,2) )%sign;
					!temp2_cpdata%sign = cpdata( BCST(bn1,bn2,3) )%sign;
					
					cpdata( BCST(bn1,bn2,2) )%bn = temp1_cpdata%bn;
					cpdata( BCST(bn1,bn2,2) )%p = temp1_cpdata%p;					
					cpdata( BCST(bn1,bn2,2) )%vn = temp1_cpdata%vn;
					cpdata( BCST(bn1,bn2,2) )%en = temp1_cpdata%en;
					
					cpdata( BCST(bn1,bn2,3) )%bn = temp2_cpdata%bn;
					cpdata( BCST(bn1,bn2,3) )%p = temp2_cpdata%p;
					cpdata( BCST(bn1,bn2,3) )%vn = temp2_cpdata%vn;
					cpdata( BCST(bn1,bn2,3) )%en = temp2_cpdata%en;
				
				!!!! 첫번째 것만 slip 일 경우
				else if( s_bit(1) .and. s_bit(2)==0 ) then
					
					c_vnd = cpdata( BCST(bn1,bn2,2) )%vn;
				
					t_vnd = temp1_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp1 );
				
					t_vnd = temp2_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp2 );
				
					if( disp1 < disp2  ) then						
						cpdata( BCST(bn1,bn2,2) )%bn = temp1_cpdata%bn;
						cpdata( BCST(bn1,bn2,2) )%p = temp1_cpdata%p;
						cpdata( BCST(bn1,bn2,2) )%vn = temp1_cpdata%vn;
						cpdata( BCST(bn1,bn2,2) )%en = temp1_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,2) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,2) )%sign;
						!cpdata( BCST(bn1,bn2,2) ) = temp1_cpdata;
						!cpdata( BCST(bn1,bn2,2) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,2) )%sign = s_bit(2); 
					else					
						cpdata( BCST(bn1,bn2,2) )%bn = temp2_cpdata%bn;
						cpdata( BCST(bn1,bn2,2) )%p = temp2_cpdata%p;
						cpdata( BCST(bn1,bn2,2) )%vn = temp2_cpdata%vn;
						cpdata( BCST(bn1,bn2,2) )%en = temp2_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,2) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,2) )%sign;
						!cpdata( BCST(bn1,bn2,2) ) = temp2_cpdata;
						!cpdata( BCST(bn1,bn2,2) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,2) )%sign = s_bit(2); 
					end if

				!!!! 두번째 것만 slip 일 경우
				else if( s_bit(1)==0 .and. s_bit(2) ) then
					
					c_vnd = cpdata( BCST(bn1,bn2,3) )%vn;
				
					t_vnd = temp1_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp1 );
				
					t_vnd = temp2_cpdata%vn;
					call get_disp_vertices(UU, c_vnd(1), c_vnd(2), t_vnd(1), t_vnd(2), disp2 );
				
					if( disp1 < disp2  ) then
						cpdata( BCST(bn1,bn2,3) )%bn = temp1_cpdata%bn;
						cpdata( BCST(bn1,bn2,3) )%p = temp1_cpdata%p;
						cpdata( BCST(bn1,bn2,3) )%vn = temp1_cpdata%vn;
						cpdata( BCST(bn1,bn2,3) )%en = temp1_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,3) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,3) )%sign;
						!cpdata( BCST(bn1,bn2,3) ) = temp1_cpdata;
						!cpdata( BCST(bn1,bn2,3) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,3) )%sign = s_bit(2); 
					else					
						cpdata( BCST(bn1,bn2,3) )%bn = temp2_cpdata%bn;
						cpdata( BCST(bn1,bn2,3) )%p = temp2_cpdata%p;
						cpdata( BCST(bn1,bn2,3) )%vn = temp2_cpdata%vn;
						cpdata( BCST(bn1,bn2,3) )%en = temp2_cpdata%en;
						!s_bit(1) = cpdata( BCST(bn1,bn2,3) )%state;
						!s_bit(2) = cpdata( BCST(bn1,bn2,3) )%sign;
						!cpdata( BCST(bn1,bn2,3) ) = temp2_cpdata;
						!cpdata( BCST(bn1,bn2,3) )%state = s_bit(1); 
						!cpdata( BCST(bn1,bn2,3) )%sign = s_bit(2); 
					end if
				
				end if
				
			end if
			
		end if

	end if

	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! Fast Contact Data 만들기
	if(0) then

	fcsN = 0;

	do i=1, TBN
		do j=i+1, TBN
			
			if( BCST(i,j,1) == 1 ) then		! Edge to corner contact 일 경우에만 FCS_data 에 들어감
				
				fcsN = fcsN + 1;
				
				! Edge block data
				FCS_data(fcsN,1) = cpdata( BCST(i,j,2) )%en(1);		! edge block number
				FCS_data(fcsN,2) = cpdata( BCST(i,j,2) )%en(2);		!   vertex number 1
				FCS_data(fcsN,3) = cpdata( BCST(i,j,2) )%en(2) + 1;	!   vertex number 2

				if( FCS_data(fcsN,3) > Block( FCS_data(fcsN,1) )%vN	) then
					FCS_data(fcsN,3) = 1;
				end if

				! Vertex block data
				FCS_data(fcsN,4) = cpdata( BCST(i,j,2) )%vn(1);		! Vertex block number
				FCS_data(fcsN,5) = cpdata( BCST(i,j,2) )%vn(2) - 1;	!   vertex number1
				FCS_data(fcsN,6) = cpdata( BCST(i,j,2) )%vn(2);		!   vertex number2
				FCS_data(fcsN,7) = cpdata( BCST(i,j,2) )%vn(2) + 1;	!   vertex number3

				if( FCS_data(fcsN,5) == 0 ) then
					FCS_data(fcsN,5) = Block( FCS_data(fcsN,4) )%vN;
				end if

				if( FCS_data(fcsN,7) > Block( FCS_data(fcsN,4) )%vN	) then
					FCS_data(fcsN,7) = 1;
				end if

			end if

		end do
	end do
	end if
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

111	continue;
			
end subroutine Determine_ContactType


! Find the edges which pass the overlapped box 
subroutine Find_Feasible_Edges(UU, BN1, BN2, ivtxN, Ovlbox, edge)
	use global_data
	implicit none
	
	real*8 :: UU( 3*TBN );		! Input : Block Displacement at the center
	integer :: BN1, BN2			! Input : Block number
	integer :: ivtxN(2,0:MVN)	! Input : Inside Vertex Number 
	real*8 :: OvlBox(4)			! Input : Overlapped box
	integer :: edge(2,0:MVN)	! Output : edge data : en1, en2, ...
	
	byte :: edge_bit(MVN);
	integer :: k1, k2, temp_edge(MVN);
	integer :: i, j, k, nn, bn(2), bi, ei, n;	
	integer :: overlapped 	
	real*8 :: xx(2), yy(2), xb(4), yb(4)
	real*8 :: edge_box(4), a(2), b(2), c1, c2
	
	edge = 0;
	bn = (/BN1, BN2/);
	do i=1,2
				
		if( ivtxN(i,0) /= 0 ) then	! vertexs are in the ovlbox
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11	
			! edge_bit 를 구성한다.			
			edge_bit = 0;
			do j=1,ivtxN(i,0)	
				
				nn = ivtxN(i,j);				
				if( nn == 1 ) then
					edge_bit( Block(bn(i))%vN ) = 1;
					edge_bit(1) = 1;
				else
					edge_bit(nn-1) = 1;
					edge_bit(nn) = 1;
				end if
			end do
			
			! feasible edge를 선택한다.
			k=1;
			do j=1,Block(bn(i))%vN
				
				if( edge_bit(j) == 1 ) then
					edge(i,k) = j;
					k = k + 1;
				end if

			end do
			edge(i,0) = k-1;
			
			! feasible edge를 정렬한다.
			do j=1, edge(i,0)-1
				
				if( (edge(i,j+1)-edge(i,j)) > 1 ) then
					
					k1 = j;		k2 = edge(i,0);				
					temp_edge = edge(i,1:MVN);
					edge(i,1:k2-k1) = temp_edge(k1+1:k2);
					edge(i,k2-k1+1:k2) = temp_edge(1:k1);
					goto 10;
					
				end if

			end do

10			continue;

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! vertexs are outside of ovlbox
		else 
			n=1;
			do j=1,Block(bn(i))%vN
				bi = j; ei = j+1;
				if( j == Block(bn(i))%vN ) then
					ei = 1;
				end	if
				
				xx = (/nBlock(bn(i))%vtx(bi,1), nBlock(bn(i))%vtx(ei,1)/) + UU(3*bn(i)-2);
				yy = (/nBlock(bn(i))%vtx(bi,2), nBlock(bn(i))%vtx(ei,2)/) + UU(3*bn(i)-1);
								
				edge_box =  (/dmin1(xx(1), xx(2)), dmax1(xx(1), xx(2)), dmin1(yy(1), yy(2)), dmax1(yy(1), yy(2))/);
			
				call Check_Overlapped_Region( edge_box, OvlBox, edge_box, overlapped )
				
				if( overlapped ) then
					
					xb = (/OvlBox(1), OvlBox(1), OvlBox(2), OvlBox(2)/);
					yb = (/OvlBox(3), OvlBox(4), OvlBox(3), OvlBox(4)/);
				
					a = (/xx(2)-xx(1), yy(2)-yy(1)/);	
					b = (/xb(1)-xx(1), yb(1)-yy(1)/);
					call cross_product(a, b, c1);
					if( dabs(c1) < err ) then; c1 = 0; end if
					
					do k=2,4
						b = (/xb(k)-xx(1), yb(k)-yy(1)/);
						
						call cross_product(a, b, c2);
						if( dabs(c2) < err ) then; c2 = 0; end if
						
						if( c1*c2 <= 0 ) then						
							edge(i,n) = j;
							n = n+1;
							goto 20;		 						
						end if
					end do
				
				end if

20				continue;

			end do

			edge(i,0) = n-1;

		end if

	end do	

end subroutine

! Check whether the vertex is inside of the Overlapped Box
subroutine Find_Vertex_Inside(UU, BN1, BN2, OvlBox, ivtxN)
	use global_data
	implicit none

	real*8 :: UU( 3*TBN );		! Input : Block Displacement at the center
	integer :: i, j, BN1, BN2, check
	real*8 :: OvlBox(4)			! Common Overlapped Box : input
	integer :: ivtxN(2,0:MVN)	! Output : Inside Vertex Number ( i block TN, n1, n2, ... ; j block )
	real*8 :: x, y, x_min, y_min, x_max, y_max

	x_min = OvlBox(1);		x_max = OvlBox(2);
	y_min = OvlBox(3);		y_max = OvlBox(4);

	ivtxN = 0;

	! Block 1
	j=1;
	do i=1,Block(BN1)%vN
		
		x = nBlock(BN1)%vtx(i,1) + UU(3*BN1-2);	
		y = nBlock(BN1)%vtx(i,2) + UU(3*BN1-1);
		
		check = 0;
		if( (x >= x_min .and. x <= x_max) .and. (y >= y_min .and. y <= y_max) ) then
			ivtxN(1,j) = i;
			j = j+1;
			check = 1;
		end if
		
		!! Control the numerical errors
		if( check == 0 ) then
			if( dabs(x-x_min) < err .or. dabs(x-x_max) < err ) then
				if( dabs(y-y_min) < err .or. dabs(y-y_max) < err ) then 
					ivtxN(1,j) = i;
					j = j+1;
				end if
			end if
		end if

	end do
	ivtxN(1,0) = j-1;

	! Block 2
	j=1;
	do i=1,Block(BN2)%vN
		
		x = nBlock(BN2)%vtx(i,1) + UU(3*BN2-2);	
		y = nBlock(BN2)%vtx(i,2) + UU(3*BN2-1);
		
		check = 0;
		if( (x >= x_min .and. x <= x_max) .and. (y >= y_min .and. y <= y_max) ) then
			ivtxN(2,j) = i;
			j = j+1;
			check = 1;
		end if
		
		!! Control the numerical errors
		if( check == 0 ) then
			if( dabs(x-x_min) < err .or. dabs(x-x_max) < err ) then
				if( dabs(y-y_min) < err .or. dabs(y-y_max) < err ) then 				
					ivtxN(2,j) = i;
					j = j+1;
				end if
			end if
		end if

	end do
	ivtxN(2,0) = j-1;

end subroutine

subroutine Check_Overlapped_Region( OutBox1, OutBox2, OverlappedBox, overlapped )
	implicit none

	integer :: i, j
	integer :: overlapped_X, overlapped_Y, overlapped
	real*8 :: OutBox1(4), OutBox2(4), OverlappedBox(4)	! Outline & Overlapped Box of Block
	real*8 :: ox_min, ox_max, oy_min, oy_max	! Overlapped Box detail
	real*8 :: xi_min, xi_max, xj_min, xj_max, yi_min, yi_max, yj_min, yj_max
	real*8 :: err = 1e-13	

	!! controll the numerical errors
	do i=1,4
		do j=1,4
			if( abs(OutBox1(i)-OutBox2(j)) < err ) then
				OutBox1(i) = OutBox2(j);
			end if
		end do
	end do
	!! 

	xi_min = OutBox1(1);		xj_min = OutBox2(1);
	xi_max = OutBox1(2);		xj_max = OutBox2(2);			
			
	yi_min = OutBox1(3);		yj_min = OutBox2(3);
	yi_max = OutBox1(4);		yj_max = OutBox2(4);			
			
	overlapped = 0;
	overlapped_X = 0;
	overlapped_Y = 0;
			
	! checking x *******************************************
	if(	xi_max <= xj_max ) then
				
		if( xi_min <= xj_min .and. xj_min <= xi_max ) then
			overlapped_X = 1;
			ox_min = xj_min;		ox_max = xi_max;
		end if

		if( xi_min >= xj_min ) then
			overlapped_X = 1;
			ox_min = xi_min;		ox_max = xi_max;
		end if
			
	else
				
		if( xj_max >= xi_min .and. xj_min <= xi_min ) then
			overlapped_X = 1;
			ox_min = xi_min;		ox_max = xj_max;
		end if

		if( xj_min >= xi_min ) then
			overlapped_X = 1;
			ox_min = xj_min;		ox_max = xj_max;
		end if
			
	end	if
						
	! checking y *******************************************
	if(	yi_max <= yj_max ) then
		
		if( yi_min <= yj_min .and. yj_min <= yi_max ) then
			overlapped_Y = 1;
			oy_min = yj_min;		oy_max = yi_max;
		end if

		if( yi_min >= yj_min ) then
			overlapped_Y = 1;
			oy_min = yi_min;		oy_max = yi_max;
		end if
			
	else
			
		if( yj_max >= yi_min .and. yj_min <= yi_min ) then
			overlapped_Y = 1;
			oy_min = yi_min;		oy_max = yj_max;
		end if

		if( yj_min >= yi_min ) then
			overlapped_Y = 1;
			oy_min = yj_min;		oy_max = yj_max;
		end if
			
	end	if
	! ******************************************************		
	
	if( overlapped_X .and. overlapped_Y ) then
		overlapped = 1;
		OverlappedBox=(/ox_min, ox_max, oy_min, oy_max/);
	end	if

end subroutine


subroutine Fast_contact_search( UU )

	use global_data
	implicit none
	
	real*8 :: UU(3*TBN);
	integer :: i
	integer :: bsn		
			
	integer :: bi, ei;		
	real*8 :: p(2), q(2)
	real*8 :: p1(2), p2(2), p3(2), p4(2);
	real*8 :: c, s, a, b, a2, b2, a3, b3, a4, b4;


	do i=1,fcsN
		
		!! 1st block (edge block)
		bsn = FCS_data(i,1);	! block number
		bi = FCS_data(i,2);		! begin vertex	
		ei = FCS_data(i,3);		! end vertex

		c = dcos( UU(3*bsn) );		s = dsin( UU(3*bsn) );

		p1 = UU(3*bsn-2:3*bsn-1);	! block center
		p2 = p1;

		p1(1) = p1(1) + c*Block(bsn)%vtx(bi,1) - s*Block(bsn)%vtx(bi,2);
		p1(2) = p1(2) + s*Block(bsn)%vtx(bi,1) + c*Block(bsn)%vtx(bi,2);

		p2(1) = p2(1) + c*Block(bsn)%vtx(ei,1) - s*Block(bsn)%vtx(ei,2);
		p2(2) = p2(2) + s*Block(bsn)%vtx(ei,1) + c*Block(bsn)%vtx(ei,2);

		!! 2nd block (vertex block): 5번, 6번 vertex 에 의한 vector
		bsn = FCS_data(i,4);	! block number
		bi = FCS_data(i,5);		! begin vertex	
		ei = FCS_data(i,6);		! end vertex

		c = dcos( UU(3*bsn) );		s = dsin( UU(3*bsn) );

		p3 = UU(3*bsn-2:3*bsn-1);	! block center
		p4 = p3;

		p3(1) = p3(1) + c*Block(bsn)%vtx(bi,1) - s*Block(bsn)%vtx(bi,2);
		p3(2) = p3(2) + s*Block(bsn)%vtx(bi,1) + c*Block(bsn)%vtx(bi,2);

		p4(1) = p4(1) + c*Block(bsn)%vtx(ei,1) - s*Block(bsn)%vtx(ei,2);
		p4(2) = p4(2) + s*Block(bsn)%vtx(ei,1) + c*Block(bsn)%vtx(ei,2);
		
		
		p = p2 - p1;
		q = p3 - p1;		
		call cross_product(p,q,a);
		q = p4 - p1;		
		call cross_product(p,q,b);

		p = p4 - p3;
		q = p2 - p3;		
		call cross_product(p,q,a2);
		q = p1 - p3;		
		call cross_product(p,q,b2);


		!! 2nd block (vertex block): 6번, 7번 vertex 에 의한 vector
		bi = FCS_data(i,6);		! begin vertex	
		ei = FCS_data(i,7);		! end vertex

		p3 = UU(3*bsn-2:3*bsn-1);	! block center
		p4 = p3;

		p3(1) = p3(1) + c*Block(bsn)%vtx(bi,1) - s*Block(bsn)%vtx(bi,2);
		p3(2) = p3(2) + s*Block(bsn)%vtx(bi,1) + c*Block(bsn)%vtx(bi,2);

		p4(1) = p4(1) + c*Block(bsn)%vtx(ei,1) - s*Block(bsn)%vtx(ei,2);
		p4(2) = p4(2) + s*Block(bsn)%vtx(ei,1) + c*Block(bsn)%vtx(ei,2);
		
		p = p2 - p1;
		q = p3 - p1;		
		call cross_product(p,q,a3);
		q = p4 - p1;		
		call cross_product(p,q,b3);


		p = p4 - p3;
		q = p2 - p3;		
		call cross_product(p,q,a4);
		q = p1 - p3;		
		call cross_product(p,q,b4);

		! Crossing 여부 검사
		if( a*b > 0 .or. a2*b2 > 0 .or. a3*b3 > 0 .or. a4*b4 > 0) then 
				Contact_search_bit = 1;
		end if

	end do
	

end subroutine





