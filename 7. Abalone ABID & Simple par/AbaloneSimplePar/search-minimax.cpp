/**
 * Very simple example strategy:
 * Search all possible positions reachable via one move,
 * and return the move leading to best position
 *
 * (c) 2006, Josef Weidendorfer
 */

#include <math.h>
#include <string.h>
#include <iostream>
#include <mpi.h>
using namespace std;

#include "search.h"
#include "board.h"
#include "eval.h"

#define print 0

/**
 * To create your own search strategy:
 * - copy this file into another one,
 * - change the class name one the name given in constructor,
 * - adjust clone() to return an instance of your class
 * - adjust last line of this file to create a global instance
 *   of your class
 * - adjust the Makefile to include your class in SEARCH_OBJS
 * - implement searchBestMove()
 *
 * Advises for implementation of searchBestMove():
 * - call foundBestMove() when finding a best move since search start
 * - call finishedNode() when finishing evaluation of a tree node
 * - Use _maxDepth for strength level (maximal level searched in tree)
 */
class MiniMaxStrategy: public SearchStrategy
{
 public:
    // Defines the name of the strategy
	 MiniMaxStrategy(): SearchStrategy("MiniMax") {}

    // Factory method: just return a new instance of this class
    SearchStrategy* clone() { return new MiniMaxStrategy(); }

    // Function such that the processes other than rank 0 do not open port
    void wait();
 private:

    /**
     * Implementation of the strategy.
     */
    void searchBestMove();
    int minimax( int depth);   
};

void MiniMaxStrategy::wait()
{
	printf("\nHello, I'm waiting in Minimax, goin to searchBestMove()\n");
    this->searchBestMove();
}

int MiniMaxStrategy::minimax( int depth)
{
	int value =0;
	int localdepth=_maxDepth;
	int CurrentColor = _board->actColor();
	
	int min = 40000;
	int max = -40000;
	
	Move m;
    MoveList list;
	
	if( depth >= localdepth )
	{	
		value = evaluate();

        // Increment no of evaluations
        _ev->incrementNoEval();

		return value;
	}
	
	depth++;
	
	// generate list of allowed moves, put them into <list>
    generateMoves(list);
    
    // loop over all moves
    while(list.getNext(m)) 
    {
    	// PLAY MOVE
    	_board->playMove(m);
   
        // recusive call until minimax do evaluate 
		value = minimax(depth);
    	
    	if( CurrentColor == _board->color1 ) // Assuming, O, Black, max
    	{
    		if( value > max )
    		{
    			max = value;   		
    		}   		
    	}
    	else // Assuming X, Red, min
    	{
    		if( value < min )
    		{
    			min = value;
    		}
    	}
    	
    	// TAKE BACK
    	_board->takeBack();
    	
    	// Finishing the evaluation of tree node
    	finishedNode(depth,&m);
    }
	
    if ( CurrentColor == _board->color1 ) // color O, blck, max
	{
    	return max;
	}
    else // Color X, red, min
    {
    	return min;
    }
}

void MiniMaxStrategy::searchBestMove()
{
   
    // Round Robin Scheduling
    if( _sc->getRank() == 0 )
    {
        char boardLayout[500];
        int value = 0;
        long min = 40000;
        long max = -40000;
        int currentDepth = 0;
        int NoMoves= 0;
        int CurrentColor = _board->actColor(); // stored own color
        int NumTask = _sc->getNumTask();
        int rank = _sc->getRank();
        long SelfEval = 0;
        static long OtherEval = 0;
        int i;

        
        // 2. send board to slave , along with move no
        int len = sprintf(boardLayout, "%s\n",_board->getState());
        for( i=1; i< NumTask; i++ )
        {
            MPI_Send( boardLayout, 500, MPI_CHAR, i, 10, MPI_COMM_WORLD );
        }

        Move m;
        MoveList list;

        // Generate All Moves
        // generate list of allowed moves, put them into <list>
        generateMoves(list);
    
    // loop over all moves
    while(list.getNext(m)) 
    {
        if( NoMoves % _sc->getNumTask() == _sc->getRank() )	
        {	   
    	if( print )
    	cout << "Rank:"<< rank <<  "\tMove:" << m.name();

        // Below code executes only for specific ranks

    	// PLAY MOVE
    	_board->playMove(m);
    	
    	value = minimax(currentDepth);
    	
    	if( print )
    	cout << "\tValue: " << value;
    	
    	if( CurrentColor == _board->color1 ) // Assuming, O, Blck, max
    	{
    		if( value > max )
    		{
    			if( print )
    			cout << "Color: " << CurrentColor << "-Maximising" << endl;
    			max = value;
			    m.value = value;
    			foundBestMove(0, m, value);
    		}    			
    	}
    	else // Assuming X, Red
    	{
    		if( value < min )
			{
    			if( print )
				cout << "Color: " << CurrentColor << "-Minimising" << endl;
				min = value;
				m.value = value;
				foundBestMove(0, m, value);
			}    	
    	}
    	
    	if( print )
    	cout << "\t\tMax: " << max << " Min: " << min << endl;
    	
    	// TAKE BACK
    	_board->takeBack();
    	
    	// Finishing the evaluation of tree node
    	finishedNode(0,&m);
	    }
	    else
    	{
	    	finishedNode(0,&m);
	    }
    	NoMoves++;
    }//while ends

    // calculating evaluations for rank 0
    SelfEval = _ev->getNoEval() - OtherEval;

	int k = 0;
    for( k=1; k< NumTask; k++ )
    {
    	int value1;
	    short field1;
    	unsigned char direction1;
	    int type1;
    	MPI_Status stat;
	    int NoEvaluations;

        // Receive moves
    	MPI_Recv( &value1, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);
	    MPI_Recv( &field1, 1, MPI_SHORT, k,10,MPI_COMM_WORLD, &stat); 
    	MPI_Recv( &direction1, 1, MPI_UNSIGNED_CHAR, k,10,MPI_COMM_WORLD, &stat);
	    MPI_Recv( &type1, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);	
    	MPI_Recv( &NoEvaluations, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);

	    _ev->addNoEval( NoEvaluations );
        OtherEval += NoEvaluations;

        // Choose the best move
	    if( CurrentColor == _board->color1 ) // O, black
    	{
	      if( value1 > max )
	      {
            max = value1; // resetting value
	        Move::MoveType t = (Move::MoveType)type1;
	        Move m1( field1, direction1, t );
	        foundBestMove(0, m1, value1);
	      }
	    }
    	else // X, red
	    {
	      if( value1 < min )
	      {
            min = value1;
	        Move::MoveType t = (Move::MoveType)type1;
	        Move m1( field1, direction1, t );
	        foundBestMove(0, m1, value1);
	      }
	    }
      }

      cout << "Rank0, SelfEval: " << SelfEval << endl;

    }
    else
    {
        char boardLayout[500];
        int value = 0;
        long min = 40000;
        long max = -40000;
        int currentDepth = 0;
        int NoMoves= 0;
        long RankEval=0;
        MPI_Status stat;

        while(1)// new search
        {
            // Set Evaluation counter to Nul
            _ev->resetEval();

            // Receive board
            MPI_Recv( boardLayout, 500, MPI_CHAR, 0, 10, MPI_COMM_WORLD, &stat );

            if( strncmp( boardLayout, "exit", 4 ) == 0)
            {
                break;
            }

            // Update board
            _board->setState( boardLayout );

            int CurrentColor = _board->actColor(); // stored own color
            Move m;
            MoveList list;
            int rank = _sc->getRank();

            // Generate All Moves
            // generate list of allowed moves, put them into <list>
            generateMoves(list);

          // loop over all moves
          while(list.getNext(m)) 
          {
            if( NoMoves % _sc->getNumTask() == _sc->getRank() )	
            {	   
    	        if( print )
    	        cout << "Rank:"<< rank <<  "\tMove:" << m.name();

            // Below code executes only for specific ranks

        	// PLAY MOVE
    	    _board->playMove(m);
    	
        	value = minimax(currentDepth);
        	
    	    if( print )
        	cout << "\tValue: " << value;
    	
        	if( CurrentColor == _board->color1 ) // Assuming, O, Blck, max
    	    {
    		    if( value > max )
    		    {
    			    if( print )
    			    cout << "Color: " << CurrentColor << "-Maximising" << endl;
    			    max = value;
			        m.value = value;
    		    	foundBestMove(0, m, value);
    		    }    			
        	}
        	else // Assuming X, Red
    	    {
    		    if( value < min )
    			{
        			if( print )
		    		cout << "Color: " << CurrentColor << "-Minimising" << endl;
			    	min = value;
				    m.value = value;
    				foundBestMove(0, m, value);
	    		}    	
    	    }
    	
        	if( print )
        	cout << "\t\tMax: " << max << " Min: " << min << endl;
    	
        	// TAKE BACK
        	_board->takeBack();
    	
    	    // Finishing the evaluation of tree node
        	finishedNode(0,&m);
	        }
	        else
    	    {
	    	    finishedNode(0,&m);
	        }
 
        	NoMoves++;
      }//while ends

        // send best move
	    MPI_Send(&_bestMove.value,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value
        MPI_Send(&_bestMove.field,1, MPI_SHORT,0, 10, MPI_COMM_WORLD ); //field
    	MPI_Send(&_bestMove.direction,1, MPI_UNSIGNED_CHAR,0, 10, MPI_COMM_WORLD );// direction
	    MPI_Send(&_bestMove.type,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // TYPE-INT

    	int NoEvaluations = _ev->getNoEval();
        RankEval += _ev->getNoEval();

	    // send No of eval
    	MPI_Send(&NoEvaluations ,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value
     } // till exit board (while)


      cout << "Total Evaluations for Rank: " << _sc->getRank() << " Evaluations: " << RankEval << endl; 
    }
}

// register ourselve as a search strategy
MiniMaxStrategy miniMaxStrategy;
