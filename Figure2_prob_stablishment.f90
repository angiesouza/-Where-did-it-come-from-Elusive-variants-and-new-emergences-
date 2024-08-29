program dois_hosts_assexuada

	! Declarando as variáveis 

	IMPLICIT REAL*8 (A-H,O-Z)
	integer, ALLOCATABLE :: M(:,:), Maumentada(:,:), Msub(:,:), Mauxrep (:,:), Mvar(:,:), MiniRep(:,:)
	real, dimension(3)  :: ntp
	real, dimension(2)  :: ntinv

   	! Parameters

	npar = 200 ! carrying capacity
	ngen = 1000 ! genome size 
	tmut = 0.001d0 ! repication rate
	b = 15 ! replication rate
	dp = 0.01 ! standart deviation
	fenotimo1in = 0.500d0 ! ph2(minimum value)
	fenotimo1f = 0.541d0 ! ph2 (maximum value)
	fenotimo2 = 0.500d0 ! ph1 (fixed value)
	deltafenotimo = 0.001d0! ph2 variation between simulations
	ntlag = 50 ! time interval between the migration event and the results report
	nrep = 1000 ! number of replicates
	ntinv = (/ 1,200/)	! tint (time to migration)
	ntp = (/ 1, 10, 100/) ! propagule size
	
	! Output
	
	OPEN(UNIT=17, FILE='Output_prob_stablishment.txt', STATUS='UNKNOWN')
	
	!
	
	NG=nint(ngen*dmax)
	nl=ngen+1
	
	call random_seed()
		
	ALLOCATE (M(npar,ngen+1), Maumentada(npar,ngen+1), Msub(npar,ngen+1), Mauxrep (npar,ngen+1),&
	Mvar(nind,ngen+2),MiniRep(1,ngen+1))
	
	do jyyy = 1, size(ntinv)
	intro = ntinv(jyyy)

		do jy = 1, size(ntp)
		nind = ntp(jy)
		
			DEALLOCATE (Mvar)
			ALLOCATE (Mvar(nind,ngen+2))
		
			fenotimo1 = fenotimo1in
			Do while (fenotimo1.le.fenotimo1f)
				
				nspfacum=0
				ngrupoacum=0
				ncrep=0 
				
				Do while (ncrep.lt.nrep)
				
					ncrep=ncrep+1 

					! CRIAR ESPAÇO INICIAL
			
					M=0 ! New host
					Mvar=0 
					Mauxrep=0
					mesp=1
					
					Msub=0 ! Donor host
					nuns=nint(ngen*fenotimo2)
					Msub(:,1:nuns+1)=1
						
					ntempo= intro+ntlag
					
					DO n=1, ntempo
				
						! INTRODUZIR PARASITOS DE ACORDO COM A PRESSÃO DE PROPÁGULOS

						If(n.eq.intro) then
			
							!Selecionando os que serão introduzidos
						
							nsubs=0 
							Mvar=0 
							
							Do i=1,nind	
							
								30 continue
								call random_number(h)
								h3=h*real(npar)
								nsort=1+int(h3)
						
								If (Msub(nsort,1).ne.0) then
									Mvar(i,1:ngen+1)=Msub(nsort,1:ngen+1) 
									Msub(nsort,1) = 0 
								else
									go to 30
								Endif
							
							ENDDO
								
							!Verificando a compatibilidade
							
							Do ik=1, nind
								nsomaloci=count(Mvar(ik,2:ngen+1).eq.1)
								fennormal= real(nsomaloci)/real(ngen) !normalizando o fenótipo
								dist= fennormal-fenotimo1
								probsob= exp(-dist**2/(2*dp**2))
								call random_number(h4)
								If (h4.le.probsob) Mvar(ik,ngen+2) = 1 !Marcando os compatíveis
							Enddo
							
							nintroduzidos=count(Mvar(:,ngen+2).eq.1)
						
							!Agora sim vamos introduzi-los na populaçao
							
							ntipoderep=0
							nespacos= count(M(:,1) .eq.0)
							nlinha = npar-nespacos+1
								
							If (nespacos.ge.nintroduzidos) then! Se eu tenho espaço de sobra eu uso aqui
								
								Do i=1, nind
									If (Mvar(i,ngen+2).eq.1) then
										M(nlinha,1:ngen+1)=Mvar(i,1:ngen+1)
										nlinha=nlinha+1
									Endif
								enddo
								
							else

								ntipoderep=1
								naumento=nintroduzidos-nespacos
								
								DEALLOCATE (Maumentada)
								ALLOCATE (Maumentada(npar+naumento,ngen+1))
								
								Maumentada=0
							
								Do iii=1, npar
									Maumentada(iii,1:ngen+1) = M(iii,1:ngen+1)
								Enddo
					
								Do i=1, nind
									If (Mvar(i,ngen+2).eq.1) then
										Maumentada(nlinha,1:ngen+1)=Mvar(i,1:ngen+1)
										nlinha=nlinha+1
									Endif
								enddo

							Endif
							
						Endif			

					
						!REPLICAÇÃO POP 1
					
						nvazios=count(M(:,1).eq.0)
						npop=npar-nvazios
						
						nprole=nint(npop*b)
						
						If (nprole.gt.npar) nprole=npar
						
						nfilhostotal=0
						Mauxrep=0
						nsobrev=0
						
						Do while (nfilhostotal.lt.nprole)
						
							IF (ntipoderep.eq.0) then
								
								call random_number(h1)
								nescolhido=1+int(h1*npop) !na população fonta é H1*npar
								MiniRep(1,:)= M(nescolhido,:)
								
							else
								
								call random_number(h1)
								nescolhido=1+int(h1*(npar+naumento)) !na população fonta é H1*npar
								MiniRep(1,:)= Maumentada(nescolhido,:)

							Endif
									
							!Mutação
								
							do j=2, ngen+1
								call random_number(h1)			
								if(h1.le.tmut) then
									if(MiniRep(1,j).eq.0) then 
										MiniRep(1,j)=1
									else
										MiniRep(1,j)=0
									endif				
								endif	
							enddo
							
							!Compatibilidade
							
							nsomaloci=count(MiniRep(1,2:ngen+1).eq.1)
							fennormal= real(nsomaloci)/real(ngen) !normalizando o fenótipo
							dist= fennormal-fenotimo1
							probsob= exp(-dist**2/(2*dp**2))
							call random_number(hii)
							nfilhostotal=nfilhostotal+1
							
							If (hii.le.probsob) then
								nsobrev=nsobrev+1
								Mauxrep(nsobrev,:)=MiniRep(1,:)
							Endif
								
						Enddo

						ntipoderep=0
						
						M=Mauxrep


						! REPLICAÇÃO POP 2 
						
						nvazios=count(Msub(:,1).eq.0)
						
						npop=npar-nvazios
						nprole=nint(npop*b)
						if(nprole.gt.npar) nprole=npar
						nfilhostotal=0
						Mauxrep=0
						nsobrev=0
					
						Do while (nfilhostotal.lt.nprole)
						
							call random_number(h1)
							nescolhido=1+int(h1*npop)
							MiniRep(1,:)= Msub(nescolhido,:)
							
							!Mutação
							
							do j=2, ngen+1
								call random_number(h1)			
								if(h1.le.tmut) then
									if(MiniRep(1,j).eq.0) then 
										MiniRep(1,j)=1
									else
										MiniRep(1,j)=0
									endif				
								endif	
							enddo
							
							!Compatibilidade
							
							nsomaloci=count(MiniRep(1,2:ngen+1).eq.1)
							fennormal= real(nsomaloci)/real(ngen) !normalizando o fenótipo
							dist= fennormal-fenotimo2
							probsob= exp(-dist**2/(2*dp**2))
							call random_number(hii)
							nfilhostotal=nfilhostotal+1
						
							If (hii.le.probsob) then
								nsobrev=nsobrev+1
								Mauxrep(nsobrev,:)=MiniRep(1,:)
							Endif
								
						Enddo
				
						Msub=Mauxrep

						If (n.eq.ntempo) then	
						
							If (M(1,1).ne.0) then
								ntemalguemnofinal=1
							else
								ntemalguemnofinal=0
							End if
					
							nsobrevfinal=count(M(:,1).ne.0)
				
							dist2= fenotimo1-fenotimo2
							
							propintro = real(nintroduzidos)/real(nind)
						
							OPEN(UNIT=17, FILE='Output_prob_stablishment.txt', STATUS='OLD', access="append")
								write(17,103) intro, b, nind, dist2/dp , propintro, ntemalguemnofinal, nsobrevfinal						
							close(17)
							
						Endif 
					 
					end do !tempo

				enddo !repetições

			fenotimo1=fenotimo1+deltafenotimo
			
			Enddo !fenotimo

		Enddo

	Enddo

print*, "Fim!"

103 FORMAT(1(2x,i5),1(3x,f5.2),1(3x,i5),2(7x,f5.2), 2(7x,i5))

end program dois_hosts_assexuada
