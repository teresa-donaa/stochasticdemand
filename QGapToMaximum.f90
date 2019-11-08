MODULE QGapToMaximum
!
USE globals
USE QL_routines
USE generic_routines
USE EquilibriumCheck
!
! Computes gap in Q function values w.r.t. maximum
!
IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeQGapToMax ( iModel )
    !
    ! Computes statistics for one model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: iModel
    !
    ! Declaring local variable
    !
    INTEGER, PARAMETER :: numThresPathCycleLength = 10
    INTEGER, PARAMETER :: ThresPathCycleLength(numThresPathCycleLength) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
    INTEGER :: iGame, iState, iAgent, iThres, i, j, CycleLengthGame, CycleStatesGame(numPeriods), &
        OptimalStrategy(numStates,numAgents), OptimalStrategyVec(lengthStrategies)
    INTEGER, DIMENSION(0:numThresPathCycleLength,0:numAgents) :: NumQGapTot, NumQGapOnPath, &
        NumQGapNotOnPath, NumQGapNotBRAllStates, NumQGapNotBRonPath, & 
        NumQGapNotEqAllStates, NumQGapNotEqonPath
    REAL(8), DIMENSION(0:numThresPathCycleLength,0:numAgents) :: SumQGapTot, SumQGapOnPath, &
        SumQGapNotOnPath, SumQGapNotBRAllStates, SumQGapNotBRonPath, SumQGapNotEqAllStates, SumQGapNotEqonPath
    REAL(8), DIMENSION(0:numAgents) :: QGapTotGame,QGapOnPathGame,QGapNotOnPathGame,QGapNotBRAllStatesGame, &
            QGapNotBRonPathGame,QGapNotEqAllStatesGame,QGapNotEqonPathGame
    !
    ! Beginning execution
    !
    PRINT*, 'Computing Q gaps'
    !
    ! Initializing variables
    !
    !$ CALL OMP_SET_NUM_THREADS(numCores)
    !
    SumQGapTot = 0.d0
    SumQGapOnPath = 0.d0
    SumQGapNotOnPath = 0.d0
    SumQGapNotBRAllStates = 0.d0
    SumQGapNotBRonPath = 0.d0
    SumQGapNotEqAllStates = 0.d0
    SumQGapNotEqonPath = 0.d0
    NumQGapTot = 0
    NumQGapOnPath = 0
    NumQGapNotOnPath = 0
    NumQGapNotBRAllStates = 0
    NumQGapNotBRonPath = 0
    NumQGapNotEqAllStates = 0
    NumQGapNotEqonPath = 0
    !
    ! Reading strategies and states at convergence from file
    !
    CALL ReadInfoModel()
    !
    ! Beginning loop over games
    !
    !$omp parallel do &
    !$omp private(OptimalStrategyVec,CycleLengthGame,CycleStatesGame,OptimalStrategy, &
    !$omp   QGapTotGame,QGapOnPathGame,QGapNotOnPathGame,QGapNotBRAllStatesGame, &
    !$omp   QGapNotBRonPathGame,QGapNotEqAllStatesGame,QGapNotEqonPathGame,iThres) &
    !$omp reduction(+ : SumQGapTot,SumQGapOnPath,SumQGapNotOnPath,SumQGapNotBRAllStates,SumQGapNotBRonPath, &
    !$omp   SumQGapNotEqAllStates,SumQGapNotEqonPath,NumQGapTot,NumQGapOnPath,NumQGapNotOnPath, &
    !$omp   NumQGapNotBRAllStates,NumQGapNotBRonPath,NumQGapNotEqAllStates,NumQGapNotEqonPath)
    !
    DO iGame = 1, numGames                  ! Start of loop aver games
        !
        PRINT*, 'iGame = ', iGame
        !
        ! Read strategy and last observed state at convergence from file
        !
        !$omp critical
        OptimalStrategyVec = indexStrategies(:,iGame)
        CycleLengthGame = CycleLength(iGame)
        CycleStatesGame(:CycleLengthGame) = CycleStates(:CycleLengthGame,iGame)
        !$omp end critical
        !
        OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
        !
        ! Compute Q gap for the optimal strategy for all agents, in all states and actions
        !
        CALL computeQGapToMaxGame(OptimalStrategy,CycleLengthGame,CycleStatesGame(:CycleLengthGame), &
            QGapTotGame,QGapOnPathGame,QGapNotOnPathGame,QGapNotBRAllStatesGame, &
            QGapNotBRonPathGame,QGapNotEqAllStatesGame,QGapNotEqonPathGame)
        !
        ! Summing by agent and threshold
        !
        iThres = MIN(CycleLength(iGame),ThresPathCycleLength(numThresPathCycleLength))
        IF (NOT(ANY(QGapTotGame .LE. 0.d0))) THEN
            !
            SumQGapTot(0,:) = SumQGapTot(0,:)+QGapTotGame
            SumQGapTot(iThres,:) = SumQGapTot(iThres,:)+QGapTotGame
            NumQGapTot(0,:) = NumQGapTot(0,:)+1
            NumQGapTot(iThres,:) = NumQGapTot(iThres,:)+1
            !
        END IF
        IF (NOT(ANY(QGapOnPathGame .LE. 0.d0))) THEN
            !
            SumQGapOnPath(0,:) = SumQGapOnPath(0,:)+QGapOnPathGame
            SumQGapOnPath(iThres,:) = SumQGapOnPath(iThres,:)+QGapOnPathGame
            NumQGapOnPath(0,:) = NumQGapOnPath(0,:)+1
            NumQGapOnPath(iThres,:) = NumQGapOnPath(iThres,:)+1
            !
        END IF
        IF (NOT(ANY(QGapNotOnPathGame .LE. 0.d0))) THEN
            !
            SumQGapNotOnPath(0,:) = SumQGapNotOnPath(0,:)+QGapNotOnPathGame
            SumQGapNotOnPath(iThres,:) = SumQGapNotOnPath(iThres,:)+QGapNotOnPathGame
            NumQGapNotOnPath(0,:) = NumQGapNotOnPath(0,:)+1
            NumQGapNotOnPath(iThres,:) = NumQGapNotOnPath(iThres,:)+1
            !
        END IF
        IF (NOT(ANY(QGapNotBRAllStatesGame .LE. 0.d0))) THEN
            !
            SumQGapNotBRAllStates(0,:) = SumQGapNotBRAllStates(0,:)+QGapNotBRAllStatesGame
            SumQGapNotBRAllStates(iThres,:) = SumQGapNotBRAllStates(iThres,:)+QGapNotBRAllStatesGame
            NumQGapNotBRAllStates(0,:) = NumQGapNotBRAllStates(0,:)+1
            NumQGapNotBRAllStates(iThres,:) = NumQGapNotBRAllStates(iThres,:)+1
            !
        END IF
        IF (NOT(ANY(QGapNotBRonPathGame .LE. 0.d0))) THEN
            !
            SumQGapNotBRonPath(0,:) = SumQGapNotBRonPath(0,:)+QGapNotBRonPathGame
            SumQGapNotBRonPath(iThres,:) = SumQGapNotBRonPath(iThres,:)+QGapNotBRonPathGame
            NumQGapNotBRonPath(0,:) = NumQGapNotBRonPath(0,:)+1
            NumQGapNotBRonPath(iThres,:) = NumQGapNotBRonPath(iThres,:)+1
            !
        END IF
        IF (NOT(ANY(QGapNotEqAllStatesGame .LE. 0.d0))) THEN
            !
            SumQGapNotEqAllStates(0,:) = SumQGapNotEqAllStates(0,:)+QGapNotEqAllStatesGame
            SumQGapNotEqAllStates(iThres,:) = SumQGapNotEqAllStates(iThres,:)+QGapNotEqAllStatesGame
            NumQGapNotEqAllStates(0,:) = NumQGapNotEqAllStates(0,:)+1
            NumQGapNotEqAllStates(iThres,:) = NumQGapNotEqAllStates(iThres,:)+1
            !
        END IF
        IF (NOT(ANY(QGapNotEqonPathGame .LE. 0.d0))) THEN
            !
            SumQGapNotEqonPath(0,:) = SumQGapNotEqonPath(0,:)+QGapNotEqonPathGame
            SumQGapNotEqonPath(iThres,:) = SumQGapNotEqonPath(iThres,:)+QGapNotEqonPathGame
            NumQGapNotEqonPath(0,:) = NumQGapNotEqonPath(0,:)+1
            NumQGapNotEqonPath(iThres,:) = NumQGapNotEqonPath(iThres,:)+1
            !
        END IF
        !
    END DO                                  ! End of loop over games
    !$omp end parallel do
    !
    ! Averaging
    !
    SumQGapTot = SumQGapTot/DBLE(NumQGapTot)
    WHERE (ISNAN(SumQGapTot)) SumQGapTot = -999.999d0
    SumQGapOnPath = SumQGapOnPath/DBLE(NumQGapOnPath)
    WHERE (ISNAN(SumQGapOnPath)) SumQGapOnPath = -999.999d0
    SumQGapNotOnPath = SumQGapNotOnPath/DBLE(NumQGapNotOnPath)
    WHERE (ISNAN(SumQGapNotOnPath)) SumQGapNotOnPath = -999.999d0
    SumQGapNotBRAllStates = SumQGapNotBRAllStates/DBLE(NumQGapNotBRAllStates)
    WHERE (ISNAN(SumQGapNotBRAllStates)) SumQGapNotBRAllStates = -999.999d0
    SumQGapNotBRonPath = SumQGapNotBRonPath/DBLE(NumQGapNotBRonPath)
    WHERE (ISNAN(SumQGapNotBRonPath)) SumQGapNotBRonPath = -999.999d0
    SumQGapNotEqAllStates = SumQGapNotEqAllStates/DBLE(NumQGapNotEqAllStates)
    WHERE (ISNAN(SumQGapNotEqAllStates)) SumQGapNotEqAllStates = -999.999d0
    SumQGapNotEqonPath = SumQGapNotEqonPath/DBLE(NumQGapNotEqonPath)
    WHERE (ISNAN(SumQGapNotEqonPath)) SumQGapNotEqonPath = -999.999d0
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Printing averages and descriptive statistics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
    IF (iModel .EQ. 1) THEN
        !
        WRITE(10006,1) &
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
            ((i, j, j = 1, numPrices), i = 1, numAgents), &
            (i, i, i, i, i, i, i, i = 1, numAgents), &
            (j, j, j, j, j, j, j, &
                (j, i, j, i, j, i, j, i, j, i, j, i, j, i, i = 1, numAgents), j = 1, numThresPathCycleLength)
1       FORMAT('Model ', &
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
            <numAgents>(<numPrices>('Ag', I1, 'Price', I2.2, ' ')), &
            '   QGapTot QGaponPath QGapNotOnPath ', &
                'QGapNotBRAllSt QGapNotBROnPath QGapNotEqAllSt QGapNotEqOnPath ', &
            <numAgents>('QGapTotAg', I1, ' QGaponPathAg', I1, ' QGapNotOnPathAg', I1, ' ', &
                'QGapNotBRAllStAg', I1, ' QGapNotBROnPathAg', I1, ' QGapNotEqAllStAg', I1, ' QGapNotEqOnPathAg', I1, ' '), &
            <numThresPathCycleLength>('   PL', I2.2, 'QGapTot PL', I2.2, 'QGaponPath PL', I2.2, 'QGapNotOnPath ', &
                'PL', I2.2, 'QGapNotBRAllSt PL', I2.2, 'QGapNotBROnPath PL', I2.2, 'QGapNotEqAllSt PL', I2.2, 'QGapNotEqOnPath ', &
            <numAgents>('PL', I2.2, 'QGapTotAg', I1, ' PL', I2.2, 'QGaponPathAg', I1, ' PL', I2.2, 'QGapNotOnPathAg', I1, ' ', &
                'PL', I2.2, 'QGapNotBRAllStAg', I1, ' PL', I2.2, 'QGapNotBROnPathAg', I1, ' PL', I2.2, 'QGapNotEqAllStAg', I1, ' PL', I2.2, 'QGapNotEqOnPathAg', I1, ' ') ) &
            )
        !
    END IF
    !
    WRITE(10006,2) codModel, &
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
        SumQGapTot(0,0), SumQGapOnPath(0,0), SumQGapNotOnPath(0,0), &
            SumQGapNotBRAllStates(0,0), SumQGapNotBRonPath(0,0), SumQGapNotEqAllStates(0,0), SumQGapNotEqonPath(0,0), &
        (SumQGapTot(0,i), SumQGapOnPath(0,i), SumQGapNotOnPath(0,i), &
            SumQGapNotBRAllStates(0,i), SumQGapNotBRonPath(0,i), SumQGapNotEqAllStates(0,i), SumQGapNotEqonPath(0,i), i = 1, numAgents), &
        (SumQGapTot(j,0), SumQGapOnPath(j,0), SumQGapNotOnPath(j,0), &
            SumQGapNotBRAllStates(j,0), SumQGapNotBRonPath(j,0), SumQGapNotEqAllStates(j,0), SumQGapNotEqonPath(j,0), &
        (SumQGapTot(j,i), SumQGapOnPath(j,i), SumQGapNotOnPath(j,i), &
            SumQGapNotBRAllStates(j,i), SumQGapNotBRonPath(j,i), SumQGapNotEqAllStates(j,i), SumQGapNotEqonPath(j,i), i = 1, numAgents), j = 1, numThresPathCycleLength)
2   FORMAT(I5, 1X, &
        <numAgents>(F10.5, 1X), <numExplorationParameters>(F10.5, 1X), <numAgents>(F10.5, 1X), &
        <numAgents>(A9, 1X, <numAgents>(F9.2, 1X)), &
        <numDemandParameters>(F10.5, 1X), &
        <6*numAgents*numMarkets>(F14.7, 1X), &
        <numPrices*numAgents>(F10.5, 1X), &
        F10.7, 1X, F10.7, 1X, F13.7, 1X, F14.7, 1X, F15.7, 1X, F14.7, 1X, F15.7, 1X, &
        <numAgents>(F10.7, 1X, F13.7, 1X, F16.7, 1X, F17.7, 1X, F18.7, 1X, F17.7, 1X, F18.7, 1X), &
        <numThresPathCycleLength>(F14.7, 1X, F14.7, 1X, F17.7, 1X, F18.7, 1X, F19.7, 1X, F18.7, 1X, F19.7, 1X, &
           <numAgents>(F14.7, 1X, F17.7, 1X, F20.7, 1X, F21.7, 1X, F22.7, 1X, F21.7, 1X, F22.7, 1X)) &
        )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeQGapToMax
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computeQGapToMaxGame ( OptimalStrategy, CycleLength, CycleStates, &
        QGapTot, QGapOnPath, QGapNotOnPath, QGapNotBRAllStates, &
        QGapNotBRonPath, QGapNotEqAllStates, QGapNotEqonPath )
    !
    ! Computes Q gap w.r.t. maximum by state for an individual replication
    !
    ! INPUT:
    !
    ! - OptimalStrategy     : strategy for all agents
    ! - CycleLength         : length of the replication's equilibrium path (i.e., state cycle)
    ! - CycleStates         : replication's equilibrium path (i.e., state cycle)
    !
    ! OUTPUT:
    !
    ! - QGapTot             : Average Q gap over all states 
    ! - QGapOnPath          : Average Q gap over cycle states
    ! - QGapNotOnPath       : Average Q gap over non-cycle states
    ! - QGapNotBRAllStates  : Average Q gap over non-best responding states
    ! - QGapNotBRonPath     : Average Q gap over non-best responding, non-cycle states
    ! - QGapNotEqAllStates  : Average Q gap over non-equilibrium states
    ! - QGapNotEqonPath     : Average Q gap over non-equilibrium cycle states
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: OptimalStrategy(numStates,numAgents)
    INTEGER, INTENT(IN) :: CycleLength
    INTEGER, INTENT(IN) :: CycleStates(CycleLength)
    REAL(8), DIMENSION(0:numAgents), INTENT(OUT) :: QGapTot, QGapOnPath, QGapNotOnPath, QGapNotBRAllStates, &
        QGapNotBRonPath, QGapNotEqAllStates, QGapNotEqonPath
    !
    ! Declaring local variables
    !
    INTEGER :: iState, iAgent, iPrice
    INTEGER :: CellVisitedStates(numPeriods), CellPreCycleLength, CellCycleLength
    REAL(8), DIMENSION(numStates,numPrices,numAgents) :: QTrue
    REAL(8), DIMENSION(numStates,numAgents) :: MaxQTrue, QGap
    LOGICAL, DIMENSION(numStates,numAgents) :: IsOnPath, IsNotOnPath, &
        IsNotBRAllStates, IsNotBROnPath, IsNotEqAllStates, IsNotEqOnPath
    LOGICAL, DIMENSION(numAgents) :: IsBR
    !
    ! Beginning execution
    !
    ! 1. Compute true Q for the optimal strategy for all agents, in all states and actions
    !
    QTrue = 0.d0
    MaxQTrue = 0.d0
    QGap = 0.d0
    !
    DO iState = 1, numStates                ! Start of loop over states
        !
        DO iAgent = 1, numAgents            ! Start of loop over agents
            !
            DO iPrice = 1, numPrices        ! Start of loop over prices
                !
                ! Compute state value function of agent iAgent for the optimal strategy in (iState,iPrice)
                !
                CALL computeQCell(OptimalStrategy,iState,iPrice,iAgent,delta, &
                    QTrue(iState,iPrice,iAgent),CellVisitedStates,CellPreCycleLength,CellCycleLength)
                !
            END DO                          ! End of loop over prices
            !
            ! Compute gap in Q function values w.r.t. maximum
            !
            MaxQTrue(iState,iAgent) = MAXVAL(QTrue(iState,:,iAgent))
            QGap(iState,iAgent) = &
                (MaxQTrue(iState,iAgent)-QTrue(iState,OptimalStrategy(iState,iAgent),iAgent))/ABS(MaxQTrue(iState,iAgent))
            !
        END DO                              ! End of loop over agents
        !
    END DO                                  ! End of loop over initial states
    !
    ! 2. Compute mask matrices
    !
    IsOnPath = .FALSE.
    IsNotOnPath = .FALSE.
    IsNotBRAllStates = .FALSE.
    IsNotBROnPath = .FALSE.
    IsNotEqAllStates = .FALSE.
    IsNotEqOnPath = .FALSE.
    !
    DO iState = 1, numStates                ! Start of loop over states
        !
        IF (ANY(CycleStates .EQ. iState)) IsOnPath(iState,:) = .TRUE.
        IF (ALL(CycleStates .NE. iState)) IsNotOnPath(iState,:) = .TRUE.
        !
        IsBR = .FALSE.
        DO iAgent = 1, numAgents            ! Start of loop over agents
            !
            IF (AreEqualReals(QTrue(iState,OptimalStrategy(iState,iAgent),iAgent),MaxQTrue(iState,iAgent))) THEN
                !
                IsBR(iAgent) = .TRUE.
                !
            ELSE
                !
                IsNotBRAllStates(iState,iAgent) = .TRUE.
                IF (ANY(CycleStates .EQ. iState)) IsNotBROnPath(iState,iAgent) = .TRUE.
                !
            END IF
            !
        END DO
        !
        DO iAgent = 1, numAgents            ! Start of loop over agents
            !
            IF (NOT(ALL(IsBR))) THEN
                !
                IsNotEqAllStates(iState,iAgent) = .TRUE.
                IF (ANY(CycleStates .EQ. iState)) IsNotEqOnPath(iState,iAgent) = .TRUE.
                !
            END IF
            !
        END DO                          ! End of loop over agents
        !
    END DO                              ! End of loop over states
    !
    ! 3. Compute Q gap averages over subsets of states
    !
    QGapTot(0) = SUM(QGap)/DBLE(numAgents*numStates)
    QGapOnPath(0) = SUM(QGap,MASK = IsOnPath)/DBLE(COUNT(IsOnPath))
    QGapNotOnPath(0) = SUM(QGap,MASK = IsNotOnPath)/DBLE(COUNT(IsNotOnPath))
    QGapNotBRAllStates(0) = SUM(QGap,MASK = IsNotBRAllStates)/DBLE(COUNT(IsNotBRAllStates))
    QGapNotBRonPath(0) = SUM(QGap,MASK = IsNotBRonPath)/DBLE(COUNT(IsNotBRonPath))
    QGapNotEqAllStates(0) = SUM(QGap,MASK = IsNotEqAllStates)/DBLE(COUNT(IsNotEqAllStates))
    QGapNotEqonPath(0) = SUM(QGap,MASK = IsNotEqonPath)/DBLE(COUNT(IsNotEqonPath))
    !
    DO iAgent = 1, numAgents
        !
        QGapTot(iAgent) = SUM(QGap(:,iAgent))/DBLE(numStates)
        QGapOnPath(iAgent) = SUM(QGap(:,iAgent),MASK = IsOnPath(:,iAgent))/DBLE(COUNT(IsOnPath(:,iAgent)))
        QGapNotOnPath(iAgent) = SUM(QGap(:,iAgent),MASK = IsNotOnPath(:,iAgent))/DBLE(COUNT(IsNotOnPath(:,iAgent)))
        QGapNotBRAllStates(iAgent) = SUM(QGap(:,iAgent),MASK = IsNotBRAllStates(:,iAgent))/DBLE(COUNT(IsNotBRAllStates(:,iAgent)))
        QGapNotBRonPath(iAgent) = SUM(QGap(:,iAgent),MASK = IsNotBRonPath(:,iAgent))/DBLE(COUNT(IsNotBRonPath(:,iAgent)))
        QGapNotEqAllStates(iAgent) = SUM(QGap(:,iAgent),MASK = IsNotEqAllStates(:,iAgent))/DBLE(COUNT(IsNotEqAllStates(:,iAgent)))
        QGapNotEqonPath(iAgent) = SUM(QGap(:,iAgent),MASK = IsNotEqonPath(:,iAgent))/DBLE(COUNT(IsNotEqonPath(:,iAgent)))
        !
    END DO
    !
    WHERE (ISNAN(QGapTot)) QGapTot = -999.999d0
    WHERE (ISNAN(QGapOnPath)) QGapOnPath = -999.999d0
    WHERE (ISNAN(QGapNotOnPath)) QGapNotOnPath = -999.999d0
    WHERE (ISNAN(QGapNotBRAllStates)) QGapNotBRAllStates = -999.999d0
    WHERE (ISNAN(QGapNotBRonPath)) QGapNotBRonPath = -999.999d0
    WHERE (ISNAN(QGapNotEqAllStates)) QGapNotEqAllStates = -999.999d0
    WHERE (ISNAN(QGapNotEqonPath)) QGapNotEqonPath = -999.999d0
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE computeQGapToMaxGame
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE QGapToMaximum
