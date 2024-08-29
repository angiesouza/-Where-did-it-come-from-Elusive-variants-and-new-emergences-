program dois_hosts_assexuada

	! Declarando as variáveis 

	IMPLICIT REAL*8 (A-H,O-Z)
	integer, ALLOCATABLE :: M(:,:), Maumentada(:,:), Msub(:,:), MS(:,:), Mauxrep (:,:), Mvar(:,:), MiniRep(:,:)
	integer :: ngrupo, mesp
	common /c3/  ncrep, nind, intro
	character*6 nomeRep,numInd, nomeIntro, nnt, nnpar, nnind, txmut, fenotipo, nomeDist


   	! Parameters

	npar = 200 ! carrying capacity
	ngen = 1000 ! genome size 
	tmut = 0.001d0 ! repication rate
	b = 15 ! replication rate
	dmax=0.05 ! maxium dissimilarity
	dp = 0.01 ! standart deviation
	fenotimo1in = 0.500d0 ! ph2(minimum value)
	fenotimo1f = 0.531d0 ! ph2 (maximum value)
	fenotimo2 = 0.500d0 ! ph1 (fixed value)
	deltafenotimo = 0.005d0! ph2 variation between simulations
	nrep = 200 ! number of replicates
	ntempo=1000 ! number of generations
	nind = 100 ! propagule size
	intro = 200 !time to migration
	
	
	
	! Output
	
	OPEN(UNIT=18, FILE='propagule_id.txt', STATUS='UNKNOWN')
	OPEN(UNIT=20, FILE='phenotypic_variability.txt', STATUS='UNKNOWN') 

	write(18,*)  "  nind", "    distancia" , "   repeticao ", "   tempo","   id_linhagem", &
	"   fenotipo", "   sobreviveu"
	write(20,*)    "  individuos ", "  distancia ", "  repeticao ", " tempo ", &
	" id_linhagem ", "  fenotipo ", "    host",  "    sucesso"

	close(18)
	close(20)
	
	!
	
	NG=nint(ngen*dmax)
	nl=ngen+1
	
	call random_seed()

	ALLOCATE (M(npar,ngen+1), Maumentada(npar,ngen+1),  Msub(npar,ngen+1), MS(npar*10000,4), &
	Mauxrep (npar,ngen+1), Mvar(nind,ngen+2),MiniRep(1,ngen+1))
	
	DEALLOCATE (Mvar)
	ALLOCATE (Mvar(nind,ngen+2))
		
			fenotimo1 = fenotimo1in
			Do while (fenotimo1.le.fenotimo1f)
			
				ncrep=0 
				
				print*, "     fenotimo", fenotimo1
				
				Do while (ncrep.lt.nrep)
				
					ncrep=ncrep+1 
					
					print*,  "     rep", ncrep

					! INITIAL SPACE
			
					M=0 ! New host
					Mvar=0 
					Mauxrep=0
					mesp=1
					
					Msub=0 ! Donor host
					nuns=nint(ngen*fenotimo2)
					Msub(:,1:nuns+1)=1
					
					
					DO n=1, ntempo
				
						! INTRODUCE PARASITES ACCORDING TO PROPAGULE PRESSURE

						If(n.eq.intro) then
			
							!Selecting individuals
						
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
								
							! Checking compatibility
							
							Do ik=1, nind
								nsomaloci=count(Mvar(ik,2:ngen+1).eq.1)
								fennormal= real(nsomaloci)/real(ngen) !normalizando o fenótipo
								dist= fennormal-fenotimo1
								probsob= exp(-dist**2/(2*dp**2))
								call random_number(h4)
								If (h4.le.probsob) Mvar(ik,ngen+2) = 1 ! Surviving individuals
							Enddo
							
							nintroduzidos=count(Mvar(:,ngen+2).eq.1)
							
							ndist = int(((fenotimo1-fenotimo2)/dp)*10)
							
							OPEN(UNIT=18, FILE='propagule_id.txt', STATUS='OLD', access="append")
								do j=1, nind	
									nfen=count(Mvar(j,2:ngen+1).eq.1)
									write(18,100)  nind, ndist, ncrep, n, Mvar(j,1), nfen, Mvar(j,ngen+2)
								enddo
							Close(18)
							
							
							!Now we will introduce them into the new host
							
							ntipoderep=0
							nespacos= count(M(:,1) .eq.0)
							nlinha = npar-nespacos+1
								
							If (nespacos.ge.nintroduzidos) then
								
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

					
						! PARASITE REPLICATION - NEW HOST
					
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
								nescolhido=1+int(h1*npop) 
								MiniRep(1,:)= M(nescolhido,:)
								
							else
								
								call random_number(h1)
								nescolhido=1+int(h1*(npar+naumento)) 
								MiniRep(1,:)= Maumentada(nescolhido,:)

							Endif
							
							!Mutation
							
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
							
							!Compatibility
							
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


						! PARASITE REPLICATION - DONOR HOST 
						
						nvazios=count(Msub(:,1).eq.0)
						
						npop=npar-nvazios
						nprole=nint(npop*b)
						if(nprole.gt.npar) nprole=npar
						nfilhostotal=0
						Mauxrep=0
						nsobrev=0
					
						Do while (nfilhostotal.lt.nprole)
						
							22 continue
							call random_number(h1)
							nescolhido=1+int(h1*npop)
							If (Msub(nescolhido,1).eq.0) go to 22
							MiniRep(1,:)= Msub(nescolhido,:)
							
							!Mutation
							
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
							
							!Compatibility
							
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
						
						dist2= (fenotimo1-fenotimo2)/dp
						
						ntemalguem=0
						
						if (n.eq.199 .or. n.eq.ntempo)	call phenotypic_variability (M, Msub, ngen,  mesp, npar,&
							ncrep, n, nind, invin, b, dist2)
						
						
					end do !tempo

				enddo !repetições

				fenotimo1=fenotimo1+deltafenotimo
			
			Enddo 

print*, "Fim!"

100 FORMAT(7(6x,i5))
102 FORMAT(1(1x,i5),1(10x,f5.2), 4(7x,i5))
103 FORMAT(2(2x,i5),2(7x,f5.2), 2(7x,i5))

end program dois_hosts_assexuada


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  SUB-ROTINA  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE phenotypic_variability (MaV, MaV2, ngen,  mesp2, npar, ncrep2, nt, nind, invin, b, dist2)
IMPLICIT REAL*8 (A-H,O-Z)
integer, ALLOCATABLE:: MaV4(:,:)   
dimension:: MaV(npar,ngen+1), MaV2(npar,ngen+1)

	ntamanho1=count(MaV(:,1).ne.0)
	ntamanho2=count(MaV2(:,1).ne.0)
	nsomatamanho=ntamanho1+ntamanho2
	
	ALLOCATE(MaV4(nsomatamanho,ngen+1)) 
	
	ntem=0
	if(MaV(1,1).ne.0) ntem=1
	
	MaV4=0
	
	MaV4(1:ntamanho1, :) = MaV (1:ntamanho1, :)
	MaV4(ntamanho1+1:nsomatamanho, :) = MaV2 (1:ntamanho2, :)

OPEN(UNIT=20, FILE='phenotypic_variability.txt', STATUS='old', access="append") 


Do iv=1,nsomatamanho
	nsomafenotimo=count(MaV4(iv,2:ngen+1).eq.1)
	If(iv.le.ntamanho1) then
		nh=1
	else
		nh=2
	endif
	write(20,108)  nind, dist2, ncrep2, nt, MaV4(iv,1), nsomafenotimo, nh, ntem
ENDDO

close(20)

108 FORMAT(1(2x,i5),1(3x,f5.2), 6(6x,i5)) 

END SUBROUTINE phenotypic_variability
