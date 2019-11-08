MODULE LearningSimulation
!
USE globals
USE QL_routines
USE omp_lib
USE ifport
!
! Computes Monte Carlo Q-Learning simulations
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeModel ( iModel, codModel, alpha, ExplorationParameters, delta )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel, codModel
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: alpha, delta
    REAL(8), DIMENSION(numExplorationParameters) :: ExplorationParameters
    !
    ! Declaring local variable
    !
    INTEGER :: idumIP, ivIP(32), iyIP, idum2IP, idum, iv(32), iy, idum2, idumQ, ivQ(32), iyQ, idum2Q
    INTEGER :: iIters, iItersFix, i, j, h, l, iGame, iItersInStrategy, convergedGame, numGamesConverged
    INTEGER :: state, statePrime, stateFix, actionPrime, market, marketPrime
    INTEGER, DIMENSION(numStates,numAgents) :: strategy, strategyPrime, strategyFix
    INTEGER :: pPrime(numAgents), p(DepthState,numAgents)
    INTEGER :: iAgent, iState, iPrice, jAgent
    INTEGER :: minIndexStrategies, maxIndexStrategies
    INTEGER(8) :: numGames_I8
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: Q
    REAL(8) :: uIniPrice(DepthState,numAgents,numGames), uExploration(2,numAgents), &
        uIniMarket, uMarket
    REAL(8) :: u(2), eps(numAgents)
    REAL(8) :: newq, oldq
    REAL(8) :: meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
    CHARACTER(len = 25) :: QFileName
    CHARACTER(len = LengthFormatTotModelsPrint) :: iGamesChar, codModelChar
    CHARACTER(len = 200) :: PTrajectoryFileName
    LOGICAL :: maskConverged(numGames)
    !
    ! Beginning execution
    !
    ! Initializing various quantities
    !
    converged = 0    
    indexStrategies = 0
    indexLastState = 0
    indexLastMarket = 0
    timeToConvergence = 0.d0
    WRITE(codModelChar,'(I0.<LengthFormatTotModelsPrint>)') codModel
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over numGames
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    ! Generating uIniPrice
    !
    idumIP = -1
    idum2IP = 123456789
    ivIP = 0
    iyIP = 0
    CALL generate_uIniPrice(uIniPrice,uIniMarket,idumIP,ivIP,iyIP,idum2IP)  
    !
    ! Starting loop over games
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !$omp parallel do &
    !$omp private(idum,iv,iy,idum2,idumQ,ivQ,iyQ,idum2Q,Q,maxValQ, &
    !$omp   strategyPrime,pPrime,p,statePrime,actionPrime,iIters,iItersFix,iItersInStrategy,convergedGame, &
    !$omp   state,stateFix,strategy,strategyFix,market,marketPrime, &
    !$omp   eps,uExploration,uMarket,u,oldq,newq,iAgent,iState,iPrice,jAgent, &
    !$omp   QFileName,iGamesChar) &
    !$omp firstprivate(numGames,PI,delta,uIniPrice,uIniMarket,ExplorationParameters,itersPerYear,alpha, &
    !$omp   itersInPerfMeasPeriod,maxIters,printQ,codModelChar)
    DO iGame = 1, numGames
        !
        PRINT*, 'Game = ', iGame, ' started'
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Learning phase
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Initializing random number generators
        !
        idum = -iGame
        idum2 = 123456789
        iv = 0
        iy = 0
        !
        idumQ = -iGame
        idum2Q = 123456789
        ivQ = 0
        iyQ = 0
        !
        ! Initializing Q matrices
        !
        !$omp critical
        CALL initQMatrices(iGame,idumQ,ivQ,iyQ,idum2Q,PI,delta,Q,maxValQ,strategyPrime)
        !$omp end critical
        strategy = strategyPrime
        !
        ! Randomly initializing prices and state
        !
        CALL initState(uIniPrice(:,:,iGame),uIniMarket,p,statePrime,actionPrime,marketPrime)
        state = statePrime
        market = marketPrime
        !
        ! Loop
        !
        iIters = 0
        iItersInStrategy = 0
        convergedGame = -1
        IF ((typeExplorationMechanism .EQ. 2) .OR. (typeExplorationMechanism .EQ. 3)) eps = ExplorationParameters(:numAgents)
        IF (typeExplorationMechanism .EQ. 4) eps = ExplorationParameters(:numAgents)
        !
        DO 
            !
            ! Iterations counter
            !
            iIters = iIters+1
            !
            ! Generating exploration random numbers
            !            
            CALL generateUExploration(uExploration,idum,iv,iy,idum2)  
            CALL generateUMarket(uMarket,idum,iv,iy,idum2)  
            !
            ! Compute pPrime by balancing exploration vs. exploitation
            !
            CALL computePPrime(ExplorationParameters,uExploration,strategyPrime,state,iIters,pPrime,Q,eps)
            !
            ! Compute marketPrime 
            !
            CALL computeMarketPrime(uMarket,market,marketPrime)
            !
            ! Defining the new state
            !
            IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
            p(1,:) = pPrime
            statePrime = computeStateNumber(p)
            actionPrime = computeActionNumber(pPrime)
            !
            ! Each agent collects his payoff and updates
            !
            DO iAgent = 1, numAgents
                !
                ! Q matrices and strategies update
                !
                oldq = Q(state,pPrime(iAgent),iAgent)
                newq = oldq+alpha(iAgent)*(PI(actionPrime,iAgent,market)+delta(iAgent)*maxValQ(statePrime,iAgent)-oldq)
                Q(state,pPrime(iAgent),iAgent) = newq
                IF (newq .GT. maxValQ(state,iAgent)) THEN
                    !
                    maxValQ(state,iAgent) = newq
                    IF (strategyPrime(state,iAgent) .NE. pPrime(iAgent)) strategyPrime(state,iAgent) = pPrime(iAgent)
                    !
                END IF
                IF ((newq .LT. maxValQ(state,iAgent)) .AND. (strategyPrime(state,iAgent) .EQ. pPrime(iAgent))) THEN
                    !
                    CALL MaxLocBreakTies(numPrices,Q(state,:,iAgent),idumQ,ivQ,iyQ,idum2Q, &
                        maxValQ(state,iAgent),strategyPrime(state,iAgent))
                    !
                END IF
                !
            END DO
            !
            ! Measuring performance
            !
            IF ((PerfMeasPeriodTime .LE. 0) .AND. ((iIters-1)/itersPerYear .GE. ABS(PerfMeasPeriodTime))) THEN
                !
                ! PerfMeasPeriodTime < 0: at convergence after mandatory training of PerfMeasPeriodTime years
                !
                IF (ALL(strategyPrime(state,:) .EQ. strategy(state,:))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                ELSE
                    !
                    iItersInStrategy = 1
                    !
                END IF
                !
            ELSE 
                !
                ! PerfMeasPeriodTime > 1: Starts right after the completion of year PerfMeasPeriodTime - 1
                !
                IF ((PerfMeasPeriodTime .GE. 1) .AND. ((iIters-1)/itersPerYear .GE. (PerfMeasPeriodTime-1))) THEN
                    !
                    iItersInStrategy = iItersInStrategy+1
                    !
                END IF
                !
            END IF
            !
            ! Check for convergence in strategy
            !
            IF (convergedGame .EQ. -1) THEN
                !
                ! Maximum number of iterations exceeded
                IF (iIters .GT. maxIters) THEN
                    !
                    convergedGame = 0
                    strategyFix = strategy
                    stateFix = state
                    iItersFix = iIters
                    !
                END IF
                !
                ! Prescribed number of iterations reached
                IF ((iItersInStrategy .EQ. itersInPerfMeasPeriod) .AND. (PerfMeasPeriodTime .GE. 1)) THEN
                    !
                    convergedGame = 0
                    strategyFix = strategy
                    stateFix = state
                    iItersFix = iIters
                    !
                END IF
                !
                ! Convergence in strategy reached
                IF ((iItersInStrategy .EQ. itersInPerfMeasPeriod) .AND. (PerfMeasPeriodTime .LE. 0)) THEN
                    !
                    convergedGame = 1
                    strategyFix = strategy
                    stateFix = state
                    iItersFix = iIters
                    !
                END IF
                !
            END IF
            !
            ! Check for loop exit criterion
            !
            IF (convergedGame .NE. -1) EXIT
            !
            ! If no convergence yet, update and iterate
            !
            strategy(state,:) = strategyPrime(state,:)
            state = statePrime
            market = marketPrime            
            !
            ! End of loop over iterations
            !
        END DO          
        !
        ! Write Q matrices to file
        !
        IF (printQ .EQ. 1) THEN
            !
            ! Open Q matrices output file
            !
            !$omp critical
            WRITE(iGamesChar,'(I0.5)') iGame
            QFileName = 'Q_' // TRIM(codModelChar) // '_' // iGamesChar // '.txt'
            !
            ! Write on Q matrices to file
            !
            OPEN(UNIT = iGame,FILE = QFileName,RECL = 10000)
            DO iAgent = 1, numAgents
                !
                DO iState = 1, numStates
                    !
                    WRITE(iGame,*) Q(iState,:,iAgent)
                    !
                END DO
                !
            END DO
            CLOSE(UNIT = iGame)
            !$omp end critical
            !
        END IF
        !
        ! Record results at convergence
        !
        converged(iGame) = convergedGame
        timeToConvergence(iGame) = DBLE(iItersFix-itersInPerfMeasPeriod)/itersPerYear
        indexLastState(:,iGame) = convertNumberBase(stateFix-1,numPrices,LengthStates)
        indexStrategies(:,iGame) = computeStrategyNumber(strategyFix)
        indexLastMarket(iGame) = market
        !
        IF (convergedGame .EQ. 1) PRINT*, 'Game = ', iGame, ' converged'
        IF (convergedGame .EQ. 0) PRINT*, 'Game = ', iGame, ' did not converge'
        !
        ! End of loop over games
        !
    END DO
    !$omp end parallel do
    !
    ! Print InfoModel file
    !
    OPEN(UNIT = 996,FILE = FileNameInfoModel,STATUS = "REPLACE")
    DO iGame = 1, numGames
        !
        WRITE(996,*) iGame
        WRITE(996,*) converged(iGame)
        WRITE(996,*) timeToConvergence(iGame)
        WRITE(996,*) indexLastState(:,iGame)
        WRITE(996,*) indexLastMarket(iGame)
        DO iState = 1, numStates
            !
            WRITE(996,*) (indexStrategies((iAgent-1)*numStates+iState,iGame), iAgent = 1, numAgents)
            !
        END DO
        !
    END DO
    CLOSE(UNIT = 996)
    !
    ! Prints the RES output file
    !
    numGamesConverged = SUM(converged)
    maskConverged = (converged .EQ. 1)
    meanNashProfit = SUM(NashProfits)/numAgents
    meanCoopProfit = SUM(CoopProfits)/numAgents
    !
    ! Time to convergence
    !
    meanTimeToConvergence = SUM(timeToConvergence,MASK = maskConverged)/numGamesConverged
    seTimeToConvergence = &
        SQRT(SUM(timeToConvergence**2,MASK = maskConverged)/numGamesConverged-meanTimeToConvergence**2)
    numGames_I8 = numGames
    CALL SORTQQ(LOC(timeToConvergence),numGames_I8,SRT$REAL8)
    medianTimeToConvergence = timeToConvergence(NINT(0.5d0*numGames))
    !
    ! Print output
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10002,891) &
            (i, i = 1, numAgents), &
            (i, i = 1, numExplorationParameters), (i, i = 1, numAgents), &
            (i, (j, i, j = 1, numAgents), i = 1, numAgents), &
            (i, i = 1, numDemandParameters), &
            ((i, j, j = 1, numMarkets), i = 1, numAgents), &
            ((i, j, j = 1, numMarkets), i = 1, numAgents), &
            ((i, j, j = 1, numMarkets), i = 1, numAgents), &
            ((i, j, j = 1, numMarkets), i = 1, numAgents), &
            ((i, j, j = 1, numMarkets), i = 1, numAgents), &
            ((i, j, j = 1, numMarkets), i = 1, numAgents), &
            ((i, j, j = 1, numPrices), i = 1, numAgents)
891     FORMAT('Model ', &
            <numAgents>('    alpha', I1, ' '), &
            <numExplorationParameters>('     beta', I1, ' '), <numAgents>('    delta', I1, ' '), &
            <numAgents>('typeQini', I1, ' ', <numAgents>('par', I1, 'Qini', I1, ' ')), &
            <numDemandParameters>('  DemPar', I0.2, ' '), &
            <numAgents>(<numMarkets>('NashPrice', I1, 'Mkt', I1, ' ')), &
            <numAgents>(<numMarkets>('CoopPrice', I1, 'Mkt', I1, ' ')), &
            <numAgents>(<numMarkets>('NashProft', I1, 'Mkt', I1, ' ')), &
            <numAgents>(<numMarkets>('CoopProft', I1, 'Mkt', I1, ' ')), &
            <numAgents>(<numMarkets>('NashMktSh', I1, 'Mkt', I1, ' ')), &
            <numAgents>(<numMarkets>('CoopMktSh', I1, 'Mkt', I1, ' ')), &
            <numAgents>(<numPrices>('Ag', I1, 'Price', I0.2, ' ')), &
            '   numConv     avgTTC      seTTC     medTTC ')
        !
    END IF
    !
    WRITE(10002,9911) codModel, &
        alpha, MExpl, delta, &
        (typeQInitialization(i), parQInitialization(i, :), i = 1, numAgents), &
        DemandParameters, &
        ((NashPrices(i,j), j = 1, numMarkets), i = 1, numAgents), &
        ((CoopPrices(i,j), j = 1, numMarkets), i = 1, numAgents), &
        ((NashProfits(i,j), j = 1, numMarkets), i = 1, numAgents), &
        ((CoopProfits(i,j), j = 1, numMarkets), i = 1, numAgents), &
        ((NashMarketShares(i,j), j = 1, numMarkets), i = 1, numAgents), &
        ((CoopMarketShares(i,j), j = 1, numMarkets), i = 1, numAgents), &
        (PricesGrids(:,i), i = 1, numAgents), &
        numGamesConverged, meanTimeToConvergence, seTimeToConvergence, medianTimeToConvergence
9911 FORMAT(I5, 1X, &
        <numAgents>(F10.5, 1X), <numExplorationParameters>(F10.5, 1X), <numAgents>(F10.5, 1X), &
        <numAgents>(A9, 1X, <numAgents>(F9.2, 1X)), &
        <numDemandParameters>(F10.5, 1X), &
        <6*numAgents*numMarkets>(F14.7, 1X), &
        <numPrices*numAgents>(F10.7, 1X), &
        I10, 1X, <3>(F10.2, 1X))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeModel
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePPrime ( ExplorationParameters, uExploration, strategyPrime, state, iIters, &
        pPrime, Q, eps )
    !
    ! Computes pPrime by balancing exploration vs. exploitation
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: ExplorationParameters(numExplorationParameters)
    REAL(8), INTENT(IN) :: uExploration(2,numAgents)
    INTEGER, INTENT(IN) :: strategyPrime(numStates,numAgents)
    INTEGER, INTENT(IN) :: state, iIters
    INTEGER, INTENT(OUT) :: pPrime(numAgents)
    REAL(8), INTENT(IN) :: Q(numStates,numPrices,numAgents)
    REAL(8), INTENT(INOUT) :: eps(numAgents)
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, iPrice
    REAL(8) :: u(2), maxQ
    REAL(8) :: probs(numPrices)
    !
    ! Beginning execution
    !
    ! 1. Greedy with probability 1-epsilon, with constant epsilon
    !
    IF (typeExplorationMechanism .EQ. 1) THEN
        !
        DO iAgent = 1, numAgents
            !
            eps(iAgent) = ExplorationParameters(iAgent)
            u = uExploration(:,iAgent)
            IF (u(1) .LE. eps(iAgent)) THEN
                !
                pPrime(iAgent) = 1+INT(numPrices*u(2))
                !
            ELSE
                !
                pPrime(iAgent) = strategyPrime(state,iAgent)
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! 2. & 3. Greedy with probability 1-epsilon, with exponentially decreasing epsilon
    !
    IF ((typeExplorationMechanism .EQ. 2) .OR. (typeExplorationMechanism .EQ. 3)) THEN
        !
        DO iAgent = 1, numAgents
            !
            IF (MExpl(numAgents+iAgent) .LT. 0.d0) THEN
                !
                pPrime(iAgent) = strategyPrime(state,iAgent)
                !
            ELSE
                !
                u = uExploration(:,iAgent)
                IF (u(1) .LE. eps(iAgent)) THEN
                    !
                    pPrime(iAgent) = 1+INT(numPrices*u(2))
                    !
                ELSE
                    !
                    pPrime(iAgent) = strategyPrime(state,iAgent)
                    !
                END IF
                eps(iAgent) = eps(iAgent)*ExplorationParameters(numAgents+iAgent)
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! 4. Boltzmann exploration
    !
    IF (typeExplorationMechanism .EQ. 4) THEN
        !
        DO iAgent = 1, numAgents
            !
            maxQ = MAXVAL(Q(state,:,iAgent))
            probs = EXP((Q(state,:,iAgent)-maxQ)/eps(iAgent))
            u = uExploration(:,iAgent)*SUM(probs)
            DO iPrice = 1, numPrices
                !
                IF (iPrice .GT. 1) probs(iPrice) = probs(iPrice)+probs(iPrice-1)
                IF (u(1) .LE. probs(iPrice)) THEN
                    !
                    pPrime(iAgent) = iPrice
                    EXIT
                    !
                END IF
                !
            END DO
            eps(iAgent) = eps(iAgent)*ExplorationParameters(numAgents+iAgent)
            !
        END DO
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computePPrime
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE LearningSimulation
