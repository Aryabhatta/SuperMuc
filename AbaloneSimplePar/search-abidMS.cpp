/**
 * A real world, sequential strategy:
 * Alpha/Beta with Iterative Deepening (ABID)
 *
 * (c) 2005, Josef Weidendorfer
 */

#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <iostream>
using namespace std;

#include "search.h"
#include "board.h"
#include "eval.h"

#define print 0

class ABIDStrategy: public SearchStrategy
{
 public:
    ABIDStrategy(): SearchStrategy("ABID", 2) {}
    SearchStrategy* clone() { return new ABIDStrategy(); }

    Move& nextMove() { return _pv[1]; }

    void wait();

 private:
    void searchBestMove();
    void slaveWorks();
    /* recursive alpha/beta search */
    int alphabeta(int depth, int alpha, int beta);

    /* prinicipal variation found in last search */
    Variation _pv;
    Move _currentBestMove;
    bool _inPV;
    int _currentMaxDepth;
};

// Function wait for ranks other than 0
void ABIDStrategy::wait()
{
    printf("\nI'm waiting in ABID, Rank %d\n", _sc->getRank());    

    int rank = _sc->getRank();
    int NumTask = _sc->getNumTask();
    char boardLayout[500];
    int RankEval = 0;
    MPI_Status stat;
    
    while(1)
    {
        // Receive board
        MPI_Recv( boardLayout, 500, MPI_CHAR, 0, 10, MPI_COMM_WORLD, &stat);

        if( strncmp( boardLayout, "exit", 4 ) == 0)
        {
            int NoEvaluation;
            NoEvaluation = _ev->getNoEval();
            MPI_Send( &NoEvaluation, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
            break;
        }

        // set board
        _board->setState( boardLayout );

        _pv.clear(_maxDepth);
        _currentBestMove.type = Move::none;
        _currentMaxDepth=1;

//        this->searchBestMove();
        cout << "Rank " << rank << "going to slaveWorks" << endl;
        this->slaveWorks();

    }

    cout << "Total Evaluations for Rank: " << _sc->getRank() << " Evaluations: " << RankEval << endl;
}

/**
 * Entry point for search
 *
 * Does iterative deepening and alpha/beta width handling, and
 * calls alpha/beta search
 */
void ABIDStrategy::searchBestMove()
{    
    int alpha = -15000, beta = 15000;
    int nalpha, nbeta, currentValue = 0;

    _pv.clear(_maxDepth);
    _currentBestMove.type = Move::none;
    _currentMaxDepth=1;

    int rank = _sc->getRank();
    int NumTask = _sc->getNumTask();
    static int OtherEval = 0;
    int SelfEval = 0;

    // if Rank 0, send board
    if( rank == 0 )
    {
        char boardLayout[500];

        // get board
        int len = sprintf( boardLayout, "%s", _board->getState() );
    
        for( int k=1; k< NumTask; k++ )
        {
            // Send board
            MPI_Send( boardLayout, 500, MPI_CHAR, k, 10, MPI_COMM_WORLD );
        }
    }
    
    
    /* iterative deepening loop */
    do {

		/* searches on same level with different alpha/beta windows */
		while(1) 
		{
	
		    nalpha = alpha, nbeta = beta;
		    _inPV = (_pv[0].type != Move::none);
	
		    if (_sc && _sc->verbose()) {
			char tmp[100];
			sprintf(tmp, "Alpha/Beta [%d;%d] with max depth %d", alpha, beta, _currentMaxDepth);
			_sc->substart(tmp);
		    }
	
		    currentValue = alphabeta(0, alpha, beta);
	
		    /* stop searching if a win position is found */
		    if (currentValue > 14900 || currentValue < -14900)
			_stopSearch = true;
	
		    /* Don't break out if we haven't found a move */
		    if (_currentBestMove.type == Move::none)
			_stopSearch = false;
	
		    if (_stopSearch) break;
	
		    /* if result is outside of current alpha/beta window,
		     * the search has to be rerun with widened alpha/beta
		     */
		    if (currentValue <= nalpha) 
		    {
				alpha = -15000;
				if (beta<15000) beta = currentValue+1;
				continue;
		    }
		    if (currentValue >= nbeta) 
		    {
				if (alpha > -15000) alpha = currentValue-1;
				beta=15000;
				continue;
		    }
		    break;
		}
	
		/* Window in both directions cause of deepening */
		alpha = currentValue - 200, beta = currentValue + 200;

        cout << "Rank 0, new window" << alpha << "," << beta << endl;
	
		if (_stopSearch) break;
	
		_currentMaxDepth++;
    }
    while(_currentMaxDepth <= _maxDepth);

    _bestMove = _currentBestMove;
    cout << "Rank 0: Current best move" << _bestMove.name() << endl;

    int exitval = -32700;    
    cout << "Rank 0 exiting from searchBestMove !!!" << endl;
    MPI_Send(&exitval,1, MPI_INT,1, 10, MPI_COMM_WORLD ); // value
    MPI_Send(&exitval,1, MPI_INT,2, 10, MPI_COMM_WORLD ); // value

    if( print )
    cout << "Best Move for Rank " << rank << " Move: " << _bestMove.name() << " Value: " << _bestMove.value << endl;
}


/*
 * Alpha/Beta search
 *
 * - first, start with principal variation
 * - depending on depth, we only do depth search for some move types
 */
int ABIDStrategy::alphabeta(int depth, int alpha, int beta)
{
    int currentValue = -14999+depth, value;
    Move m;
    MoveList list;
    bool depthPhase, doDepthSearch;

    int rank = _sc->getRank();
    int NumTask = _sc->getNumTask();
    int NoMoves = 0;

    int value1, value2;
    Move m1, m2;
    int iDepthSearch1, iDepthSearch2;

    MPI_Status stat;

    /* We make a depth search for the following move types... */
    int maxType = (depth < _currentMaxDepth-1)  ? Move::maxMoveType :(depth < _currentMaxDepth)    ? Move::maxPushType :Move::maxOutType;

    _board->generateMoves(list);

    if (_sc && _sc->verbose()) 
    {
	    char tmp[100];
	    sprintf(tmp, "Alpha/Beta [%d;%d], %d moves (%d depth)", alpha, beta,
		    list.count(Move::none), list.count(maxType));
	    _sc->startedNode(depth, tmp);
    }

    /* check for an old best move in principal variation */
    if (_inPV) 
    {
		m = _pv[depth];
	
		if ((m.type != Move::none) &&
		    (!list.isElement(m, 0, true)))
		    m.type = Move::none;
	
		if (m.type == Move::none) _inPV = false;
    }

    // first, play all moves with depth search
    depthPhase = true;

    while (1) 
    {
		// get next move
		if (m.type == Move::none) 
		{
	            if (depthPhase)
        			depthPhase = list.getNext(m, maxType);
	            if (!depthPhase)
		    	if (!list.getNext(m, Move::none)) break;
		}
        
        // we could start with a non-depth move from principal variation
        doDepthSearch = depthPhase && (m.type <= maxType);
       
        if( rank == 0 && depth == 0  )
        {
//            cout << "Move at Rank 0, depth 0: " << m.name() << endl;
            if( NoMoves % 2 == 0)
            {
                m1 = m;
                iDepthSearch1 = (doDepthSearch)?1:0;
            }
            else
            {
                m2 = m;
                iDepthSearch2 = (doDepthSearch)?1:0;
            }

            if( NoMoves%2==1 )
            {
                // distribute work among slaves

                // Send to 1
                MPI_Send(&m1.value,1, MPI_INT,1, 10, MPI_COMM_WORLD ); // value
                MPI_Send(&m1.field,1, MPI_SHORT,1, 10, MPI_COMM_WORLD ); //field
                MPI_Send(&m1.direction,1, MPI_UNSIGNED_CHAR,1, 10, MPI_COMM_WORLD );// direction
                MPI_Send(&m1.type,1, MPI_INT,1, 10, MPI_COMM_WORLD ); // TYPE-INT
                MPI_Send(&iDepthSearch1, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
                MPI_Send(&alpha, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
                MPI_Send(&beta, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
                MPI_Send( &_currentMaxDepth, 1, MPI_INT, 1, 10, MPI_COMM_WORLD);
//                cout << "Move: " << m1.name() << " sent to rnak 1"<< endl;

                // Send to 2
                MPI_Send(&m2.value,1, MPI_INT,2, 10, MPI_COMM_WORLD ); // value
                MPI_Send(&m2.field,1, MPI_SHORT,2, 10, MPI_COMM_WORLD ); //field
                MPI_Send(&m2.direction,1, MPI_UNSIGNED_CHAR,2, 10, MPI_COMM_WORLD );// direction
                MPI_Send(&m2.type,1, MPI_INT,2, 10, MPI_COMM_WORLD ); // TYPE-INT
                MPI_Send(&iDepthSearch2, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
                MPI_Send(&alpha, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
                MPI_Send(&beta, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
                MPI_Send( &_currentMaxDepth, 1, MPI_INT, 2, 10, MPI_COMM_WORLD);
//                cout << "Move: " << m2.name() << " sent to rnak 2"<< endl;
            }
        }

        if( rank != 0 ) // rank 0 never does any work
        {
            _board->playMove(m);

            /* check for a win position first */
            if (!_board->isValid()) 
            {

                /* Shorter path to win position is better */
                value = 14999-depth;
            }
               else 
            {
                if (doDepthSearch) 
                {
                    /* opponent searches for its maximum; but we want the
                     * minimum: so change sign (for alpha/beta window too!)
                     */
                    value = -alphabeta(depth+1, -beta, -alpha);
                }
                else 
                {
                    value = evaluate();
                    _ev->incrementNoEval();
                }
            }

            _board->takeBack();
        }

        if( rank == 0 && depth == 0 && NoMoves%2==1 )
        {
            // get value from slaves
            MPI_Recv(&value1, 1, MPI_INT,1, 10, MPI_COMM_WORLD, &stat );
            MPI_Recv(&value2, 1, MPI_INT,2, 10, MPI_COMM_WORLD, &stat );

            cout << "Rank 1- Move: " << m1.name()  << " Value: " << value1 << " Currentvalue:" << currentValue << endl;
            cout << "Rank 2- Move: " << m2.name()  << " Value: " << value2 << " Currentvalue:" << currentValue << endl;

            if( value1 >= value2 )
            {
                value = value1;
                m = m1;
            }
            else
            {
                value = value2;
                m = m2;
            }
            /* best move so far? */
            if (value > currentValue) 
            {
                cout << "--->Chosen move: " << m.name() << " Chosen value:" << value << endl;
                currentValue = value;
                _pv.update(depth, m);

                if (_sc) _sc->foundBestMove(depth, m, currentValue);
                if (depth == 0)
                {
                    _currentBestMove = m;
                    // Insert value in currentbestmove to distiguish other best moves
                    _currentBestMove.value = currentValue;
                }

                /* alpha/beta cut off or win position ... */
                if (currentValue>14900 || currentValue >= beta) 
                {
                    if (_sc) _sc->finishedNode(depth, _pv.chain(depth));
                    return currentValue;
                }

                /* maximize alpha */
                if (currentValue > alpha) alpha = currentValue;
            }
        }

        if( rank != 0 )
        {
        /* best move so far? */
        if (value > currentValue) 
        {
            currentValue = value;
            _pv.update(depth, m);

            if (_sc) _sc->foundBestMove(depth, m, currentValue);
            if (depth == 0)
            {
                _currentBestMove = m;
                // Insert value in currentbestmove to distiguish other best moves
                _currentBestMove.value = currentValue;
            }

            /* alpha/beta cut off or win position ... */
            if (currentValue>14900 || currentValue >= beta) 
            {
                if (_sc) _sc->finishedNode(depth, _pv.chain(depth));
                return currentValue;
            }

            /* maximize alpha */
            if (currentValue > alpha) alpha = currentValue;
        }
        }

        NoMoves++; // Increment moves that I'm working on
	
		if (_stopSearch) break; // depthPhase=false;
		m.type = Move::none;
    }

    if (_sc) _sc->finishedNode(depth, _pv.chain(depth));

    return currentValue;
}

void ABIDStrategy::slaveWorks()
{
    bool doDepthSearch;
    int iDepthSearch;
    int depth = 0;
    int value;
    int alpha,beta;
    int rank = _sc->getRank();

    // Recv move to process
    int value1;
    short field1;
    unsigned char direction1;
    int type1;
    MPI_Status stat;
while(1)
{
    // Receive moves
    MPI_Recv( &value1, 1, MPI_INT, 0,10,MPI_COMM_WORLD, &stat);

    // check exti conditions
    if( value1 == -32700)
    {
        cout << "Rank:" << rank << "Exiting...." << endl;
        break;
    }

    MPI_Recv( &field1, 1, MPI_SHORT, 0,10,MPI_COMM_WORLD, &stat);
    MPI_Recv( &direction1, 1, MPI_UNSIGNED_CHAR, 0,10,MPI_COMM_WORLD, &stat);
    MPI_Recv( &type1, 1, MPI_INT, 0,10,MPI_COMM_WORLD, &stat);

    // Create Move
    Move::MoveType t = (Move::MoveType)type1;
    Move m( field1, direction1, t );
    
//    cout << "Rank: " << rank << " Move recieved: " << m.name() << endl;
    // Receive value for depthsearch
    MPI_Recv( &iDepthSearch, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );

    // Receive alpha beta
    MPI_Recv( &alpha, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );
    MPI_Recv( &beta, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );
    MPI_Recv( &_currentMaxDepth,1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );

    if( iDepthSearch == 1 )
        doDepthSearch = true;
    else
        doDepthSearch = false;

    // Dont Receive depth as it is always 0

    _board->playMove(m);

    /* check for a win position first */
    if (!_board->isValid()) 
    {

        /* Shorter path to win position is better */
        value = 14999-depth;
    }
    else 
    {
        if (doDepthSearch) 
        {
            /* opponent searches for its maximum; but we want the
             * minimum: so change sign (for alpha/beta window too!)
             */
            value = -alphabeta(depth+1, -beta, -alpha);
        }
        else 
        {
            value = evaluate();
            _ev->incrementNoEval();
        }
    }
	
    _board->takeBack();

    // Send value computed
    MPI_Send( &value, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
  }
}

// register ourselve
ABIDStrategy abidStrategy;
