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
    int alphabetaZero(int depth, int alpha, int beta);
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

    int alpha = -15000, beta = 15000;
    int rank = _sc->getRank();
    int NumTask = _sc->getNumTask();
    char boardLayout[500];
    MPI_Status stat;
    
    _pv.clear(_maxDepth);
    _currentBestMove.type = Move::none;
    _currentMaxDepth=1;
    cout << "Rank " << rank << "going to slaveWorks" << endl;
    this->slaveWorks();
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
    
    /* iterative deepening loop */
    do {
        cout << endl<< endl;
        if( print )
        cout << "CurrentMaxDepth: " << _currentMaxDepth << endl;

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
	
		    currentValue = alphabetaZero(0, alpha, beta);
	
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

        if( print )
        {
            cout << "Rank 0, new window" << alpha << "," << beta << endl;
            cout << "Rank 0: best move" << _currentBestMove.name() << endl;
        }
	
		if (_stopSearch) break;
	
		_currentMaxDepth++;
    }
    while(_currentMaxDepth <= _maxDepth);

    _bestMove = _currentBestMove;
    cout << "Rank 0: Final best move" << _bestMove.name() << endl;

    if( print )
    cout << "Best Move for Rank " << rank << " Move: " << _bestMove.name() << " Value: " << _bestMove.value << endl;
}


/*
 * Alpha/Beta search
 *
 * - first, start with principal variation
 * - depending on depth, we only do depth search for some move types
 */
int ABIDStrategy::alphabetaZero(int depth, int alpha, int beta)
{

    int currentValue = -14999+depth, value=0;
    Move m;
    MoveList list;
    bool depthPhase, doDepthSearch;
    Move PVMove;

    int rank = _sc->getRank();
    int NumTask = _sc->getNumTask();
    int NoMoves = 0;
    int value1;

    int iDepthSearch;
    bool Parallel = true;

    MPI_Status stat;

    /* We make a depth search for the following move types... */
    int maxType = (depth < _currentMaxDepth-1)  ? Move::maxMoveType :(depth < _currentMaxDepth)? Move::maxPushType :Move::maxOutType;

    _board->generateMoves(list);

    if (_sc && _sc->verbose()) 
    {
	    char tmp[100];
	    sprintf(tmp, "Alpha/Beta [%d;%d], %d moves (%d depth)", alpha, beta,
		    list.count(Move::none), list.count(maxType));
	    _sc->startedNode(depth, tmp);
    }

    PVMove = _pv[depth];

    if ((PVMove.type != Move::none) &&
        (!list.isElement(PVMove, 0, true)))
        PVMove.type = Move::none;

    if (PVMove.type == Move::none) _inPV = false;

    // Print PVMOVE
    if( print )
    cout << "PVMOVE: " << PVMove.name() << " depth: " << depth << endl;

    // first, play all moves with depth search
    depthPhase = true;
    m = PVMove;

    // get next move
    if (m.type == Move::none) 
    {
        if (depthPhase)
            depthPhase = list.getNext(m, maxType);
        if( !depthPhase)
            depthPhase = list.getNext(m, Move::none);
    }
        
    _board->playMove(m);

    /* check for a win position first */
    if (!_board->isValid()) 
    {

        /* Shorter path to win position is better */
        value = 14999-depth;
    }
    else 
    {
        if( depth < _currentMaxDepth )
        {
            /* opponent searches for its maximum; but we want the
             * minimum: so change sign (for alpha/beta window too!)
             */
            if( print )
            cout << "Goin Recusive ?" << endl;
            value = -alphabetaZero(depth+1, -beta, -alpha);
        }
        else
        {
            if( print )
            cout << "AlphabetaZero evaluate, depth: " << depth << endl; 
            value = evaluate();
            _ev->incrementNoEval();
            Parallel = false;
        }
    }

    _board->takeBack();

    /* best move so far? */
    if (value > currentValue) 
    {
        if(print)
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
        if( currentValue>14900 || currentValue >= beta ) 
        {
            if (_sc) _sc->finishedNode(depth, _pv.chain(depth));
            if( print )
            cout << "currentval At Maxdepth(SEQ): " << currentValue  <<  " > Beta: " << beta << endl;
            return currentValue;
        }

        /* maximize alpha */
        if (currentValue > alpha) alpha = currentValue;
    }
    
    if( print )   
    cout << "currentval At Maxdepth(SEQ): " << currentValue << endl;

    depthPhase = true; 
    int alphaR, alphaBest;
    int betaR, betaBest;
    Move SlaveMoves[10];
//    Move * SlaveMoves = new Move[(_pv.maxDepth-depth)];

    Move mS[NumTask-1];
    int iDepthSearchS[ NumTask -1 ];    
    NoMoves = 0;
    
    if( Parallel )
    {
        while(1)
        {
        
            if (depthPhase)
                depthPhase = list.getNext(m, maxType);
            if( !depthPhase)
            if(!list.getNext(m, Move::none))break;

            // Accumulate Moves & iDepthSearch values
            mS[ (NoMoves%(NumTask-1))] = m;          
            iDepthSearchS[ (NoMoves%(NumTask-1)) ] = (depthPhase)?1:0;

            // if sufficient Moves accumulated
            if(NoMoves%(NumTask-1)==(NumTask-2))
            {
                char boardLayout[500];

                // get board
                int len = sprintf( boardLayout, "%s", _board->getState() );

                for( int k=1; k< NumTask; k++ )
                {
                    // Send board
                    MPI_Send( boardLayout, 500, MPI_CHAR, k, 10, MPI_COMM_WORLD );
                
                    if( print )
                    cout << "Board Sent to Slave" << endl;

                    MPI_Send(&mS[k-1].value,1, MPI_INT,k, 10, MPI_COMM_WORLD ); // value
                    MPI_Send(&mS[k-1].field,1, MPI_SHORT,k, 10, MPI_COMM_WORLD ); //field
                    MPI_Send(&mS[k-1].direction,1, MPI_UNSIGNED_CHAR,k, 10, MPI_COMM_WORLD );// direction
                    MPI_Send(&mS[k-1].type,1, MPI_INT,k, 10, MPI_COMM_WORLD ); // TYPE-INT
                
//                    iDepthSearch = (doDepthSearch)?1:0;
                    MPI_Send(&iDepthSearchS[k-1], 1, MPI_INT, k, 10, MPI_COMM_WORLD);
                    MPI_Send(&alpha, 1, MPI_INT, k, 10, MPI_COMM_WORLD);
                    MPI_Send(&beta, 1, MPI_INT, k, 10, MPI_COMM_WORLD);
                    MPI_Send( &_currentMaxDepth, 1, MPI_INT, k, 10, MPI_COMM_WORLD);
                    MPI_Send( &depth,1, MPI_INT, k, 10, MPI_COMM_WORLD);
                }

                for( int k =1; k< NumTask; k++ )
                {
                    MPI_Recv(&value1, 1, MPI_INT,k, 10, MPI_COMM_WORLD, &stat );

        //            cout << "Before Recv" << endl;
                    MPI_Recv( &SlaveMoves, sizeof(Move)* _pv.maxDepth, MPI_BYTE, k, 10, MPI_COMM_WORLD, &stat);
        //            cout << "After Recv" << endl;

        //            MPI_Send( &SlaveMoves, sizeof(Move)* (_pv.maxDepth-depth), MPI_BYTE, k, 10, MPI_COMM_WORLD);

                    // recv move chain

                    if( print )
                    {
                        cout << "Value received from slave 1: " << value1 << "\tMove: " << SlaveMoves[depth].name() << endl;
                        cout << "Value now: " << value << " CurrentValue: " << currentValue <<  endl;
                    }

//                    m = SlaveMoves[depth];
                    m = mS[k-1];

                    /* best move so far? */
                    if (value > currentValue) 
                    {
                        if(print)
                        cout << "--->Chosen move: " << m.name() << " Chosen value:" << value << endl;

                        currentValue = value;
                        _pv.update(depth, m);

                        // Update chain 
                        memcpy( &(_pv.move[depth]), &SlaveMoves, sizeof(Move) * _pv.maxDepth );
        //                memcpy( &(_pv.move[depth]), SlaveMoves, sizeof(Move) * (_pv.maxDepth-depth) );
                
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
            }

            NoMoves++;

            if(_stopSearch) break;
        }
    }

    if( print )
    cout << "NoMoves: " << NoMoves << endl;

    if (_sc) _sc->finishedNode(depth, _pv.chain(depth));

    if( print )
    cout << "currentval: " << currentValue << endl;
    return currentValue;
}

void ABIDStrategy::slaveWorks()
{
    bool doDepthSearch;
    int currentValue = -14999;
    int iDepthSearch;
    int depth = 0;
    int value;
    int alpha,beta;
    int rank = _sc->getRank();
    char boardLayout[500];

    // Recv move to process
    int value1;
    short field1;
    unsigned char direction1;
    int type1;
    MPI_Status stat;

while(1)
{
    // Receive board
    MPI_Recv( boardLayout, 500, MPI_CHAR, 0, 10, MPI_COMM_WORLD, &stat);

    if( strncmp( boardLayout, "exit", 4 ) == 0)
    {
        long NoEvaluation;
        NoEvaluation = _ev->getNoEval();
        MPI_Send( &NoEvaluation, 1, MPI_LONG, 0, 10, MPI_COMM_WORLD);
        break;
    }

    // set board
    _board->setState( boardLayout );

    // Receive moves
    MPI_Recv( &value1, 1, MPI_INT, 0,10,MPI_COMM_WORLD, &stat);

    // check exit conditions
    if( value1 == -32700)
    {
        cout << "Rank:" << rank << "Exiting...." << endl;
    }

    MPI_Recv( &field1, 1, MPI_SHORT, 0,10,MPI_COMM_WORLD, &stat);
    MPI_Recv( &direction1, 1, MPI_UNSIGNED_CHAR, 0,10,MPI_COMM_WORLD, &stat);
    MPI_Recv( &type1, 1, MPI_INT, 0,10,MPI_COMM_WORLD, &stat);

    // Create Move
    Move::MoveType t = (Move::MoveType)type1;
    Move m( field1, direction1, t );

//    if( print )    
//    cout << "Rank: " << rank << " Move recieved: " << m.name() << endl;

    // Receive value for depthsearch
    MPI_Recv( &iDepthSearch, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );

    // Receive alpha beta
    MPI_Recv( &alpha, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );
    MPI_Recv( &beta, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );
    MPI_Recv( &_currentMaxDepth,1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );

    // Receive depth
    MPI_Recv( &depth,1, MPI_INT, 0, 10, MPI_COMM_WORLD, &stat );

    if( iDepthSearch == 1 )
        doDepthSearch = true;
    else
        doDepthSearch = false;

//    doDepthSearch = true;
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

    // best move so far?
    if (value > currentValue)
    {
        currentValue = value;
        _pv.update(depth, m);
    
        if (_sc)
            _sc->foundBestMove(depth, m, currentValue);

        //if (depth == 0)
        _currentBestMove = m;
     }


    // Send value computed
    MPI_Send( &value, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);

    // send move chain to update the move chain
    MPI_Send( _pv.move[depth], sizeof(Move)* _pv.maxDepth, MPI_BYTE, 0, 10, MPI_COMM_WORLD);
//    cout << "Move chain sent " << endl;
//    MPI_Send( _pv.chain(depth), sizeof(Move)* (_pv.maxDepth-depth), MPI_BYTE, 0, 10, MPI_COMM_WORLD);
  }
}

int ABIDStrategy::alphabeta(int depth, int alpha, int beta)
{
    if( print )
    cout << "Slave1 in aplhabeta" << endl;

    int currentValue = -14999+depth, value;
    Move m;
    MoveList list;
    bool depthPhase, doDepthSearch;

    int rank = _sc->getRank();
    int NumTask = _sc->getNumTask();

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
	
		if (_stopSearch) break; // depthPhase=false;
		m.type = Move::none;
    }

    if (_sc) _sc->finishedNode(depth, _pv.chain(depth));

    return currentValue;
}

// register ourselve
ABIDStrategy abidStrategy;
