MODULE QL_routines
!
USE globals
!
! Various routines used in exp2
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE initQMatrices ( iGame, idumQ, ivQ, iyQ, idum2Q, PI, delta, Q, maxValQ, maxLocQ )
    !
    ! Initializing Q matrices
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iGame
    INTEGER, INTENT(INOUT) :: idumQ, ivQ(32), iyQ, idum2Q
    REAL(8), DIMENSION(numActions,numAgents,numMarkets), INTENT(IN) :: PI
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    REAL(8), DIMENSION(numStates,numPrices,numAgents), INTENT(OUT) :: Q
    INTEGER, DIMENSION(numStates,numAgents), INTENT(OUT) :: maxLocQ
    REAL(8), DIMENSION(numStates,numAgents), INTENT(OUT) :: maxValQ
    !
    ! Declaring local variables
    !
    INTEGER :: iAgent, jAgent, iPrice, iState, i, h, status
    INTEGER :: tied(numPrices), Strategy(numStates,numAgents)
    INTEGER :: VisitedStates(numPeriods), PreCycleLength, CycleLength
    INTEGER, DIMENSION(DepthState,numAgents) :: p
    REAL(8) :: den, u
    CHARACTER(len = 225) :: QFileName
    CHARACTER(len = 5) :: iChar
    CHARACTER(len = 5) :: codModelChar
    CHARACTER(len = 200) :: QFileFolderNameAgent
    !
    ! Beginning execution
    !
    DO iAgent = 1, numAgents
        !
        IF (typeQInitialization(iAgent) .EQ. 'F') THEN
            !
            ! Agent iAgent assumes agent jAgent plays "parQInitialization(iAgent,jAgent)"
            !
            DO jAgent = 1, numAgents
                !
                Strategy(:,jAgent) = NINT(parQInitialization(iAgent,jAgent))
                !
            END DO
            DO iState = 1, numStates            ! Start of loop over states
                !
                ! Compute state value function for Strategy in iState, for all prices
                !
                DO iPrice = 1, numPrices            ! Start of loop over prices to compute a row of Q
                    !
                    CALL computeQCell(Strategy,iState,iPrice,iAgent,delta, &
                        Q(iState,iPrice,iAgent),VisitedStates,PreCycleLength,CycleLength)
                    !
                END DO                              ! End of loop over prices to compute a row of Q
                !
            END DO                              ! End of loop over states
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'G') THEN
            !
            ! Assume Grim Trigger strategies: every agent behaves as:
            ! - If all agents played "parQInitialization(iAgent,1)", keep playing "parQInitialization(iAgent,1)"
            ! - Otherwise, play "parQInitialization(iAgent,2)"
            ! 
            DO jAgent = 1, numAgents
                !
                Strategy(:,jAgent) = NINT(parQInitialization(iAgent,2))
                !
            END DO
            p = parQInitialization(iAgent,1)
            Strategy(computeStateNumber(p),:) = parQInitialization(iAgent,1)
            DO iState = 1, numStates            ! Start of loop over states
                !
                ! Compute state value function for Strategy in iState, for all prices
                !
                DO iPrice = 1, numPrices            ! Start of loop over prices to compute a row of Q
                    !
                    CALL computeQCell(Strategy,iState,iPrice,iAgent,delta, &
                        Q(iState,iPrice,iAgent),VisitedStates,PreCycleLength,CycleLength)
                    !
                END DO                              ! End of loop over prices to compute a row of Q
                !
            END DO                              ! End of loop over states
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'O') THEN
            !
            ! Randomize over the opponents decisions
            !
            DO iPrice = 1, numPrices
                !
                den = COUNT(indexActions(:,iAgent) .EQ. iPrice)*(1.d0-delta(iAgent))
                Q(:,iPrice,iAgent) = &
                    SUM(ExpectedPI(:,iAgent),MASK = indexActions(:,iAgent) .EQ. iPrice)/DBLE(den)
                !
            END DO
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'T') THEN
            !
            ! Start from a randomly drawn Q matrix at convergence 
            ! on model "parQInitialization(iAgent,1)"
            !
            
            WRITE(codModelChar,'(I0.<LengthFormatTotModelsPrint>)') NINT(parQInitialization(iAgent,1))
            i = 1+INT(DBLE(numGames)*ran2(idumQ,ivQ,iyQ,idum2Q))
            WRITE(iChar,'(I0.5)') i
            QFileName = 'Q_' // TRIM(codModelChar) // '_' // iChar // '.txt'
            QFileFolderNameAgent = QFileFolderName(iAgent)
            QFileName = TRIM(QFileFolderNameAgent) // TRIM(QFileName)
            !
            ! Read Q matrices from file
            !
            OPEN(UNIT = iGame,FILE = QFileName,READONLY,RECL = 10000,IOSTAT = status)
            IF (iAgent .GT. 1) READ(iGame,100)
100         FORMAT(<(iAgent-1)*numStates-1>(/))
            DO iState = 1, numStates
                !
                READ(iGame,*) Q(iState,:,iAgent)
                !
            END DO
            CLOSE(UNIT = iGame)
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'R') THEN
            !
            ! Randomly initialized Q matrix using a uniform distribution between 
            ! parQInitialization(iAgent,1) and parQInitialization(iAgent,2)
            !
            DO iState = 1, numStates
                !
                DO iPrice = 1, numPrices
                    !
                    Q(iState,iPrice,iAgent) = ran2(idumQ,ivQ,iyQ,idum2Q)
                    !
                END DO
                !
            END DO
            Q(:,:,iAgent) = parQInitialization(iAgent,1)+ &
                (parQInitialization(iAgent,2)-parQInitialization(iAgent,1))*Q(:,:,iAgent)
            !
        ELSE IF (typeQInitialization(iAgent) .EQ. 'U') THEN
            !
            ! Constant Q matrix with all elements set to parQInitialization(iAgent,1)
            !
            Q(:,:,iAgent) = parQInitialization(iAgent,1)
            !
        END IF
        !
    END DO
    !
    ! Find initial optimal strategy
    !
    DO iAgent = 1, numAgents
        !
        DO iState = 1, numStates
            !
            CALL MaxLocBreakTies(numPrices,Q(iState,:,iAgent),idumQ,ivQ,iyQ,idum2Q, &
                maxValQ(iState,iAgent),maxLocQ(iState,iAgent))
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initQMatrices
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE initState ( uPrice, uMarket, p, stateNumber, actionNumber, marketNumber )
    !
    ! Randomly initializing prices 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: uPrice(DepthState,numAgents)
    REAL(8), INTENT(IN) :: uMarket
    INTEGER, DIMENSION(DepthState,numAgents), INTENT(OUT) :: p
    INTEGER, INTENT(OUT) :: stateNumber, actionNumber, marketNumber
    !
    ! Beginning execution
    !
    p = 1+INT(numPrices*uPrice)
    stateNumber = computeStateNumber(p)
    actionNumber = computeActionNumber(p(1,:))
    !
    ! Computing the initial market number
    !
    IF (numMarkets .EQ. 1) THEN
        !
        ! One market size
        marketNumber = 1
        !
    ELSE
        !
        ! Any other number of market sizes
        marketNumber = 1+INT(DBLE(numMarkets)*uMarket)
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE initState
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generate_uIniPrice ( uIniPrice, uIniMarket, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(OUT) :: uIniPrice(DepthState,numAgents,numGames)
    REAL(8), INTENT(OUT) :: uIniMarket
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Declaring local variables
    !
    INTEGER :: iGame, iAgent, iDepth
    !
    ! Beginning execution
    !
    ! Generate U(0,1) draws for price initialization
    !
    DO iGame = 1, numGames
        !
        DO iDepth = 1, DepthState
            !
            DO iAgent = 1, numAgents
                !
                uIniPrice(iDepth,iAgent,iGame) = ran2(idum,iv,iy,idum2)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    ! Generate U(0,1) draws for market initialization
    !
    uIniMarket = ran2(idum,iv,iy,idum2)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE generate_uIniPrice
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generateUExploration ( uExploration, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(OUT) :: uExploration(2,numAgents)
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Declaring local variables
    !
    INTEGER :: iDecision, iAgent
    !
    ! Beginning execution
    !
    ! Generate U(0,1) draws for price initialization
    !
    DO iDecision = 1, 2
        !
        DO iAgent = 1, numAgents
            !
            uExploration(iDecision,iAgent) = ran2(idum,iv,iy,idum2)
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE generateUExploration
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE generateUMarket ( uMarket, idum, iv, iy, idum2 )   
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(OUT) :: uMarket
    INTEGER, INTENT(INOUT) :: idum
    INTEGER, INTENT(INOUT) :: iv(32)
    INTEGER, INTENT(INOUT) :: iy
    INTEGER, INTENT(INOUT) :: idum2
    !
    ! Beginning execution
    !
    ! Generate U(0,1) draw for market transition
    !
    uMarket = ran2(idum,iv,iy,idum2)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE generateUMarket
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeMarketPrime ( uMarket, market, marketPrime )
    !
    ! Computes marketPrime from market and transition probability probA0
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: uMarket
    INTEGER, INTENT(IN) :: market
    INTEGER, INTENT(OUT) :: marketPrime
    !
    ! Declaring local variables
    !
    INTEGER :: i
    REAL(8) :: sumProbA0
    !
    ! Beginning execution
    !
    IF (numMarkets .EQ. 1) THEN
        !
        ! One market size
        marketPrime = 1
        !
    ELSE 
        !
        ! Any other number of market sizes
        sumProbA0 = 0.d0
        DO i = 1, numMarkets
            !
            sumProbA0 = sumProbA0+probA0(market,i)
            IF (uMarket .LE. sumProbA0) THEN
                !
                marketPrime = i
                EXIT
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeMarketPrime
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStateNumber ( p )
    !
    ! Given the price vectors, computes the state number
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(DepthState,numAgents), INTENT(IN) :: p
    !
    ! Declaring function's type
    !
    INTEGER :: computeStateNumber
    !
    ! Declaring local variables
    !
    INTEGER, DIMENSION(LengthStates) :: stateVector
    !
    ! Beginning execution
    !
    IF (DepthState0 .GT. 0) THEN
        !
        stateVector = RESHAPE(TRANSPOSE(p),(/ LengthStates /))
        computeStateNumber = 1+SUM(cStates*(stateVector-1))
        !
    ELSE IF (DepthState0 .EQ. 0) THEN
        !
        computeStateNumber = 1
        !
    END IF
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStateNumber
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeActionNumber ( p )
    !
    ! Given the prices, computes the action number
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(numAgents), INTENT(IN) :: p
    !
    ! Declaring function's type
    !
    INTEGER :: computeActionNumber
    !
    ! Declaring local variables
    !
    INTEGER, DIMENSION(numAgents) :: tmp
    !
    ! Beginning execution
    !
    tmp = cActions*(p-1)
    computeActionNumber = 1+SUM(tmp)
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeActionNumber
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStatesCodePrint ( )
    !
    ! Compute the states code in printable format (with '.')
    !
    IMPLICIT NONE
    !
    ! Declaring function's type
    !
    CHARACTER(len = LengthFormatStatesPrint) :: computeStatesCodePrint(numStates)
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, indexState(LengthStates)
    CHARACTER(len = lengthFormatActionPrint) :: tmp
    CHARACTER(len = LengthFormatStatesPrint) :: labelState
    !
    ! Beginning execution
    !
    DO i = 1, numStates
        !
        indexState = convertNumberBase(i-1,numPrices,LengthStates)
        !
        DO j = 1, LengthStates
            !
            WRITE(tmp,'(I0.<lengthFormatActionPrint>)') indexState(j)
            IF (j .EQ. 1) THEN 
                !
                labelState = TRIM(tmp)   
                !
            ELSE IF (MOD(j,numAgents) .NE. 1) THEN
                !
                labelState = TRIM(labelState) // '.' // TRIM(tmp)   
                !
            ELSE
                !
                labelState = TRIM(labelState) // '-' // TRIM(tmp)   
                !
            END IF
            !
        END DO
        !
        computeStatesCodePrint(i) = labelState
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStatesCodePrint
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    FUNCTION computeStrategyNumber ( maxLocQ )
    !
    ! Given the maxLocQ vectors, computes the lengthStrategies-digit strategy number
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, DIMENSION(numStates,numAgents), INTENT(IN) :: maxLocQ
    !
    ! Declaring function's type
    !
    INTEGER :: computeStrategyNumber(lengthStrategies)
    !
    ! Declaring local variables
    !
    INTEGER :: i, il, iu
    !
    ! Beginning execution
    !
    iu = 0
    DO i = 1, numAgents
        !
        il = iu+1
        iu = iu+numStates
        computeStrategyNumber(il:iu) = maxLocQ(:,i)
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END FUNCTION computeStrategyNumber
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeQCell ( OptimalStrategy, iState, iPrice, iAgent, delta, &
        QCell, VisitedStates, PreCycleLength, CycleLength )
    !
    ! Computes a cell of the 'true' (i.e., theoretical) Q matrix
    !
    ! INPUT:
    !
    ! - OptimalStrategy     : strategy for all agents
    ! - iState              : current state
    ! - iPrice              : price (i.e., action) index
    ! - iAgent              : agent index
    ! - delta               : discount factors
    !
    ! OUTPUT:
    !
    ! - QCell               : 'theoretical'/'true' Q(iState,iPrice,iAgent)
    ! - VisitedStates       : numPeriods array of states visited (0 after start of cycling)
    ! - PreCycleLength      : number of periods in the pre-cycle phase
    ! - CycleLength         : number of periods in the cycle phase
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: OptimalStrategy(numStates,numAgents)
    INTEGER, INTENT(IN) :: iState
    INTEGER, INTENT(IN) :: iPrice
    INTEGER, INTENT(IN) :: iAgent
    REAL(8), DIMENSION(numAgents), INTENT(IN) :: delta
    REAL(8), INTENT(OUT) :: QCell
    INTEGER, DIMENSION(numPeriods), INTENT(OUT) :: VisitedStates
    INTEGER, INTENT(OUT) :: PreCycleLength, CycleLength
    !
    ! Declaring local variable
    !
    INTEGER :: iPeriod, p(DepthState,numAgents), pPrime(numAgents)
    REAL(8) :: VisitedProfits(numPeriods), PreCycleProfit, CycleProfit
    !
    ! Beginning execution
    !
    ! Initial p and pPrime, including deviation to iPrice
    !
    p = RESHAPE(convertNumberBase(iState-1,numPrices,numAgents*DepthState),(/ DepthState,numAgents /))
    pPrime = OptimalStrategy(iState,:)
    pPrime(iAgent) = iPrice
    !
    ! Loop over deviation period
    !
    VisitedStates = 0
    VisitedProfits = 0.d0
    DO iPeriod = 1, numPeriods
        !
        IF (DepthState .GT. 1) p(2:DepthState,:) = p(1:DepthState-1,:)
        p(1,:) = pPrime
        VisitedStates(iPeriod) = computeStateNumber(p)
        VisitedProfits(iPeriod) = ExpectedPI(computeActionNumber(pPrime),iAgent)
        !
        ! Check if the state has already been visited
        !
        IF ((iPeriod .GE. 2) .AND. (ANY(VisitedStates(:iPeriod-1) .EQ. VisitedStates(iPeriod)))) THEN
            !
            PreCycleLength = MINVAL(MINLOC((VisitedStates(:iPeriod-1)-VisitedStates(iPeriod))**2))
            CycleLength = iPeriod-PreCycleLength
            EXIT
            !
        END IF
        !
        ! After period 1, every agent follows the optimal strategy
        !
        pPrime = OptimalStrategy(VisitedStates(iPeriod),:)
        !
    END DO
    !
    ! 2. Compute state value function for the optimal strategy
    !
    PreCycleProfit = SUM(DiscountFactors(0:PreCycleLength-1,iAgent)*VisitedProfits(1:PreCycleLength))
    CycleProfit = SUM(DiscountFactors(0:CycleLength-1,iAgent)*VisitedProfits(PreCycleLength+1:iPeriod))
    Qcell = PreCycleProfit+delta(iAgent)**PreCycleLength*CycleProfit/(1.d0-delta(iAgent)**CycleLength)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeQCell
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE ReadInfoModel ( )
    !
    ! Reads the InfoModel txt file
    !
    ! INPUT:
    !
    ! None
    !
    ! OUTPUT (via global variables):
    !
    ! - converged           : numGames array, = 1 if replication converged, = 0 otherwise
    ! - timeToConvergence   : numGames array of number of iterations to convergence (/ItersPerYear)
    ! - CycleLength         : numGames array of the length of the cycles at convergence
    ! - CycleStates         : numPeriods x numGames array of states in the cycles at convergence
    ! - CyclePrices         : numAgents x numPeriods x numGames array of prices in the cycles at convergence
    ! - CycleProfits        : numAgents x numPeriods x numGames array of profits in the cycles at convergence
    ! - indexStrategies     : lengthStrategies x numGames array of strategies at convergence 
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: iGame, rGame, iCycle, iState, iAgent
    !
    ! Beginning execution
    !
    OPEN(UNIT = 998,FILE = FileNameInfoModel,STATUS = "OLD")    ! Open InfoModel file    
    DO iGame = 1, numGames
        !
        IF (MOD(iGame,100) .EQ. 0) PRINT*, 'Read ', iGame, ' strategies'
        READ(998,*) rGame
        READ(998,*) converged(iGame)
        READ(998,*) timeToConvergence(iGame)
        READ(998,*) CycleLength(iGame)
        READ(998,*) CycleStates(:CycleLength(iGame),iGame)
        READ(998,*) ((CyclePrices(iAgent,iCycle,iGame), iCycle = 1, CycleLength(iGame)), iAgent = 1, numAgents)
        READ(998,*) ((CycleProfits(iAgent,iCycle,iGame), iCycle = 1, CycleLength(iGame)), iAgent = 1, numAgents)
        DO iState = 1, numStates
            !
            READ(998,*) (indexStrategies((iAgent-1)*numStates+iState,iGame), iAgent = 1, numAgents)
            !
        END DO
        !
    END DO
    CLOSE(UNIT = 998)                   ! Close indexStrategies txt file
    PRINT*, 'Finished reading InfoModel'
    !
    ! Ending execution and returning control
    !
END SUBROUTINE ReadInfoModel
! 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
! End of execution
!
END MODULE QL_routines