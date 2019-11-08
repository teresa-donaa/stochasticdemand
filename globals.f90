MODULE globals
!
USE generic_routines
!
! Declare global parameters and variables 
!
IMPLICIT NONE
!
! Parameters
!
INTEGER, PARAMETER :: numShockPeriodsPrint = 10
INTEGER, PARAMETER :: numThresCycleLength = 10
!
! Variables
!
INTEGER :: numModels, totModels, numCores, numGames, itersPerYear, maxNumYears, maxIters, &
    itersInPerfMeasPeriod, printQ, codModel, PerfMeasPeriodTime, numPrices, numMarkets, &
    typeExplorationMechanism, DepthState0, DepthState, LengthStates, numStates, lengthStrategies, &
    LengthFormatStatesPrint, LengthFormatActionPrint, LengthFormatTotModelsPrint, &
    typePayoffInput, numAgents, numActions, numDemandParameters, numPeriods, &
    numExplorationParameters, SwitchQLearningResults, SwitchConvergenceResults, &
    SwitchImpulseResponseToBR, SwitchImpulseResponseToNash, SwitchImpulseResponseToAll, &
    SwitchEquilibriumCheck, SwitchQGapToMaximum, SwitchDetailedAnalysis
REAL(8) :: PerfMeasPeriodLength, meanNashProfit, meanCoopProfit, gammaSinghVives, SlackOnPath, SlackOffPath
CHARACTER(len = 50) :: ModelNumber, FileNameInfoModel, ModelName
!
INTEGER, ALLOCATABLE :: converged(:), indexStrategies(:,:), indexLastState(:,:), indexLastMarket(:), &
    CycleLength(:), CycleStates(:,:), CyclePrices(:,:,:), &
    indexActions(:,:), cStates(:), cActions(:), &
    indexNashPrices(:,:), indexCoopPrices(:,:)
REAL(8), ALLOCATABLE :: timeToConvergence(:), CycleProfits(:,:,:), probA0(:,:), &
    NashProfits(:,:), CoopProfits(:,:), maxValQ(:,:), NashPrices(:,:), CoopPrices(:,:), &
    ExpectedNashPrices(:), ExpectedCoopPrices(:), ExpectedNashProfits(:), ExpectedCoopProfits(:), &
    PI(:,:,:), PIQ(:,:,:), avgPI(:,:), avgPIQ(:,:), ExpectedPI(:,:), &
    PG(:,:,:), PGQ(:,:,:), avgPG(:,:), avgPGQ(:,:), &
    alpha(:), delta(:), DiscountFactors(:,:), & 
    DemandParameters(:), MExpl(:), ExplorationParameters(:), &
    NashMarketShares(:,:), CoopMarketShares(:,:), PricesGrids(:,:), &
    parQInitialization(:,:)
CHARACTER(len = :), ALLOCATABLE :: labelStates(:)
CHARACTER(len = :), ALLOCATABLE :: QFileFolderName(:)
CHARACTER(len = 1), ALLOCATABLE :: typeQInitialization(:)
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readBatchVariables ( unitNumber )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unitNumber
    !
    ! Declaring local variables
    !
    INTEGER :: i, iAgent, iPrices, jPrices, iState, iAction, jState
    INTEGER, ALLOCATABLE :: switchedState(:), indA1(:), indA2(:)
    !
    ! Beginning execution
    !
    READ(unitNumber,'(1X)') 
    READ(unitNumber,*) numModels, totModels
    READ(unitNumber,'(1X)') 
    READ(unitNumber,*) numCores
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numGames
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) itersPerYear
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) maxNumYears
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) PerfMeasPeriodTime
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) PerfMeasPeriodLength
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numAgents
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) DepthState0
    DepthState = MAX(1,DepthState0)          ! Accomodates the DepthState = 0 case
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numPrices
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) numMarkets
    !
    ! Global variables
    !
    LengthFormatTotModelsPrint = 1+INT(LOG10(DBLE(totModels)))
    maxIters = maxNumYears*itersPerYear
    itersInPerfMeasPeriod = INT(PerfMeasPeriodLength*itersPerYear)
    LengthStates = MAX(1,numAgents*DepthState0)
    LengthFormatStatesPrint = LengthStates*(1+FLOOR(LOG10(DBLE(numPrices))))+LengthStates-1
    numStates = numPrices**(numAgents*DepthState0)
    numPeriods = numStates+1
    numActions = numPrices**numAgents   ! Actions contain combinations of prices;
                                        ! they coincide with states when DepthState == 1
    lengthStrategies = numAgents*numStates
    lengthFormatActionPrint = FLOOR(LOG10(DBLE(numPrices)))+1
    READ(unitNumber,'(1X)')
    !
    ! Read type of exploration mechanism
    !
    READ(unitNumber,*) typeExplorationMechanism
    numExplorationParameters = 2*numAgents           
    READ(unitNumber,'(1X)')
    !
    ! Read type of payoff input
    !
    READ(unitNumber,*) typePayoffInput
    IF (typePayoffInput .EQ. 2) THEN
        !
        IF (numMarkets .EQ. 1) numDemandParameters = numMarkets+2*numAgents+1+2     ! a0, ai, ci, sigma, extend
        IF (numMarkets .EQ. 3) numDemandParameters = numMarkets+2*numAgents+1+2+1   ! a0, ai, ci, sigma, extend, probA0
        !
    END IF
    IF (typePayoffInput .EQ. 3) THEN
        !
        IF (numMarkets .EQ. 1) numDemandParameters = numMarkets+2*numAgents+1+2     ! a0, ai, ci, sigma = 0, extend
        IF (numMarkets .EQ. 3) numDemandParameters = numMarkets+2*numAgents+1+2+1   ! a0, ai, ci, sigma = 0, extend, probA0
        !
    END IF
    READ(unitNumber,'(1X)')
    !
    ! Continue reading input settings
    !
    READ(unitNumber,*) SwitchQLearningResults
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchConvergenceResults
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchImpulseResponseToBR
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchImpulseResponseToNash
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchImpulseResponseToAll
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchEquilibriumCheck, SlackOnPath, SlackOffPath
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchQGapToMaximum
    READ(unitNumber,'(1X)')
    READ(unitNumber,*) SwitchDetailedAnalysis
    READ(unitNumber,'(1X)')
    !
    ! Allocating matrices and vectors
    !
    ALLOCATE(converged(numGames),timeToConvergence(numGames), &
        indexStrategies(lengthStrategies,numGames),indexLastState(lengthStates,numGames), &
        indexLastMarket(numGames), CycleLength(numGames),CycleStates(numPeriods,numGames), &
        CyclePrices(numAgents,numPeriods,numGames),CycleProfits(numAgents,numPeriods,numGames), &
        indexActions(numActions,numAgents), &
        cStates(LengthStates),cActions(numAgents),DiscountFactors(0:numStates,numAgents), &
        maxValQ(numStates,numAgents), DemandParameters(numDemandParameters), &
        ExplorationParameters(numExplorationParameters), MExpl(numExplorationParameters), &
        alpha(numAgents),delta(numAgents),NashProfits(numAgents,numMarkets),CoopProfits(numAgents,numMarkets), &
        probA0(numMarkets,numMarkets), &
        PI(numActions,numAgents,numMarkets),PIQ(numActions,numAgents,numMarkets), &
        avgPI(numActions,numMarkets),avgPIQ(numActions,numMarkets),ExpectedPI(numActions,numAgents), &
        PG(numActions,numAgents,numMarkets),PGQ(numActions,numAgents,numMarkets), &
        avgPG(numActions,numMarkets),avgPGQ(numActions,numMarkets), &
        indexNashPrices(numAgents,numMarkets),indexCoopPrices(numAgents,numMarkets), &
        NashPrices(numAgents,numMarkets),CoopPrices(numAgents,numMarkets), &
        ExpectedNashPrices(numAgents),ExpectedCoopPrices(numAgents), &
        ExpectedNashProfits(numAgents),ExpectedCoopProfits(numAgents), &
        typeQInitialization(numAgents),parQInitialization(numAgents,numAgents), &
        NashMarketShares(numAgents,numMarkets),CoopMarketShares(numAgents,numMarkets), &
        PricesGrids(numPrices,numAgents))
    ALLOCATE(CHARACTER(len = 3+LengthFormatStatesPrint) :: labelStates(numStates))
    ALLOCATE(CHARACTER(len = 200) :: QFileFolderName(numAgents))
    !
    cStates = (/ (numPrices**i, i = LengthStates-1, 0, -1) /)
    cActions = (/ (numPrices**i, i = numAgents-1, 0, -1) /)
    !
    ! Actions contain the most recent prices of all agents. Actions and States
    ! coincide when DepthState == 1
    !
    DO iAction = 1, numActions
        !
        indexActions(iAction,:) = convertNumberBase(iAction-1,numPrices,numAgents)
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readBatchVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE closeBatch ( )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    DEALLOCATE(converged,timeToConvergence,indexStrategies,indexLastState, &
    indexLastMarket,CycleLength,CycleStates,CyclePrices,CycleProfits, &
    indexActions,cStates,cActions,DiscountFactors,maxValQ,DemandParameters, &
    ExplorationParameters, MExpl,alpha,delta,NashProfits,CoopProfits, &
    probA0,PI,PIQ,avgPI,avgPIQ,ExpectedPI,PG,PGQ,avgPG,avgPGQ, &
    indexNashPrices,indexCoopPrices,NashPrices,CoopPrices, &
    ExpectedNashPrices,ExpectedCoopPrices,ExpectedNashProfits,ExpectedCoopProfits, &
    typeQInitialization,parQInitialization,NashMarketShares,CoopMarketShares, &
    PricesGrids,labelStates,QFileFolderName)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE closeBatch
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readModelVariables ( unitNumber )
    !
    ! Reads input variables
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unitNumber
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, iMarket, iAgent
    REAL(8) :: pr
    !
    ! Beginning execution
    !
    READ(unitNumber,*) codModel, printQ, alpha, MExpl, delta, &
        DemandParameters, &
        ((NashPrices(iAgent,iMarket), iMarket = 1, numMarkets), iAgent = 1, numAgents), &
        ((CoopPrices(iAgent,iMarket), iMarket = 1, numMarkets), iAgent = 1, numAgents), &
        (typeQInitialization(i), parQInitialization(i,:), i = 1, numAgents)
    !
    IF (typeExplorationMechanism .EQ. 2) THEN
        !
        ExplorationParameters(:numAgents) = MExpl(:numAgents)
        ExplorationParameters(numAgents+1:) = EXP(-MExpl(numAgents+1:)/DBLE(itersPerYear))
        !
    ELSE IF (typeExplorationMechanism .EQ. 3) THEN
        !
        ExplorationParameters(:numAgents) = MExpl(:numAgents)
        ExplorationParameters(numAgents+1:) = -DBLE(itersPerYear)/DBLE(numAgents+1)* &
            LOG(1.d0-(DBLE(numPrices-1)/DBLE(numPrices))**numAgents/(DBLE(numStates*numPrices)*MExpl(numAgents+1:)))
        ExplorationParameters(numAgents+1:) = EXP(-ExplorationParameters(numAgents+1:)/DBLE(itersPerYear))
        !
    ELSE IF (typeExplorationMechanism .EQ. 4) THEN
        !
        ExplorationParameters(:numAgents) = MExpl(:numAgents)
        ExplorationParameters(numAgents+1:) = 1.d0-1.d1**MExpl(numAgents+1:)
        !
    END IF
    DiscountFactors = TRANSPOSE(RESHAPE((/ (delta**i, i = 0, numPeriods-1) /),(/ numAgents,numPeriods /)))
    ExpectedNashPrices = SUM(NashPrices,DIM = 2)/DBLE(numMarkets)
    ExpectedCoopPrices = SUM(CoopPrices,DIM = 2)/DBLE(numMarkets)
    !
    ! Constructing matrix of transition probabilities for market sizes
    !
    pr = (DBLE(numMarkets-1)*DemandParameters(numDemandParameters)+1.d0)/DBLE(numMarkets)
    DO i = 1, numMarkets
        !
        DO j = 1, numMarkets
            !
            IF (i .EQ. j) THEN
                !
                probA0(i,i) = pr
                !
            ELSE 
                !
                probA0(i,j) = (1.d0-pr)/DBLE(numMarkets-1)
                !
            END IF
            !
        END DO
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE readModelVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE globals