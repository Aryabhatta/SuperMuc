/**
 * Very simple example strategy:
 * Search all possible positions reachable via one move,
 * and return the move leading to best position
 *
 * (c) 2006, Josef Weidendorfer
 */

#include <math.h>
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
	//cout << "Color:" << _board->actColor() << endl;
	int value =0;
	int localdepth=_maxDepth;
	int CurrentColor = _board->actColor();
	
	int min = 40000;
	int max = -40000;
	
	Move m;
    MoveList list;
	
	//if( depth >= _maxDepth )
	if( depth >= localdepth )
	{	
		value = evaluate();
        // Increment no of evaluations
        _ev->incrementNoEval();

		//cout << "MaxDepth: " << localdepth << "  Aborting at this depth, value = " << value << endl;
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
    else
    {
    	return min;
    }
}

void MiniMaxStrategy::searchBestMove()
{
    if( _sc->getRank() ==  0 )
    {
    	if( print )
	    cout << "\nIn searchBestMove()\n";

        int i=0;
    
	    long min = 20000;
    	long max = -20000;
	    int currentDepth = 0;
    	int NoMoves= 0;
	    int CurrentColor = _board->actColor(); // stored own color
        char boardLayout[500];
	
    	int valueR;
    	short fieldR;
    	unsigned char directionR;
    	int typeR;
    	MPI_Status stat;
	    int NoEvaluationsR;
        int rank = _sc->getRank();

        Move m;
        MoveList list;

        // Generate All Moves
        // generate list of allowed moves, put them into <list>
        generateMoves(list);

        NoMoves = list.count();
        cout << "Total No of Moves" << NoMoves << endl;
        cout << "Rank: " << _sc->getRank() << endl;
    
        // Partition Moves
        // Assuming 1 slave process
        // loop over all moves
        while(list.getNext(m)) 
        {
            if( print )
    	    cout << "Move: " << m.name();
 
    	    // 1.PLAY MOVE
       	    _board->playMove(m);

            // 2. send board to slave , along with move no
            int len = sprintf(boardLayout, "%s\n",_board->getState());

            MPI_Send( boardLayout, 500, MPI_CHAR, 1, 10, MPI_COMM_WORLD );

        	// 3,TAKE BACK
        	_board->takeBack();
        
            // 4. Get result, evaluate value, maximise/minimise
            // Receive moves
        	MPI_Recv(&valueR, 1, MPI_INT, 1,10,MPI_COMM_WORLD, &stat);
  //      	MPI_Recv(&NoEvaluationsR, 1, MPI_INT, 1,10,MPI_COMM_WORLD, &stat);

    //    	_ev->addNoEval( NoEvaluationsR );

            // Choose the best move
	        if( CurrentColor == _board->color1 ) // O, black
	        {
	            if( valueR > max )
                {
                     max = valueR;
                     foundBestMove(0, m, valueR);
	            }
	        }
        	else // X, red
	        {
	            if( valueR < min )
	            {
                    min = valueR;
                    foundBestMove(0, m, valueR);
	            }
	        }

        } // end of while loop for moves
    }
    else
    {
        char boardLayout[500];
    	MPI_Status stat;
        int rank = _sc->getRank();
        int value=0;
        int currentDepth = 0;

        while(1)
        {
        // Receive board
        MPI_Recv( boardLayout, 500, MPI_CHAR, 0, 10, MPI_COMM_WORLD, &stat );

        if( strncmp( boardLayout, "exit", 4 ) == 0)
        {
            break;
        }

        // Update board
        _board->setState( boardLayout );

        // call minimax()
       	value = minimax(currentDepth);

        // send results
        // send best move
	    MPI_Send(&value,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value
//	    int NoEvaluations = _ev->getNoEval();

	    // send No of eval
//	    MPI_Send(&NoEvaluations ,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value
        }
    }





















/*
	if( print )
	cout << "\nIn searchBestMove()\n";

	int value = 0;
	long min = 20000;
	long max = -20000;
	int currentDepth = 0;
	int NoMoves= 0;
	int CurrentColor = _board->actColor(); // stored own color
	
    Move m;
    MoveList list;

    // generate list of allowed moves, put them into <list>
    generateMoves(list);
    
    // loop over all moves
    while(list.getNext(m)) 
    {
    	if( print )
    	cout << "Move: " << m.name();
 
    	NoMoves++;
        if( NoMoves % _sc->getNumTask() == _sc->getRank() )	
        {	   

//i	cout << "Rank: " << _sc->getRank() << " running Move: " << NoMoves << endl;
        // Below code executes only for specific ranks

    	// PLAY MOVE
    	_board->playMove(m);
    	
    	value = minimax(currentDepth);
    	
    	if( print )
    	cout << " Value: " << value;
    	
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
    	cout << " Max: " << max << " Min: " << min << endl;
    	
    	// TAKE BACK
    	_board->takeBack();
    	
    	// Finishing the evaluation of tree node
    	finishedNode(0,&m);
	}
	else
	{
	    	finishedNode(0,&m);
	}
    }
   
    if( print )
    cout << "Color:" << CurrentColor << " No of Moves: " << NoMoves << endl;

    if( print )
    cout << "Rank:" << ((CurrentColor==_board->color1)?"O":"X")<<_sc->getRank() << " No of Moves: " << NoMoves << endl;

    // print individual best moves
    cout << "Best move for rank: " << ((CurrentColor==_board->color1)?"O":"X") << \
    _sc->getRank() << " Move: " << _bestMove.name() << " Value:" << _bestMove.value << endl;

    // Send receive best moves
    if( _sc->getRank() != 0 )
    {
        // send best move
	MPI_Send(&_bestMove.value,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value
	MPI_Send(&_bestMove.field,1, MPI_SHORT,0, 10, MPI_COMM_WORLD ); //field
	MPI_Send(&_bestMove.direction,1, MPI_UNSIGNED_CHAR,0, 10, MPI_COMM_WORLD );// direction
	MPI_Send(&_bestMove.type,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // TYPE-INT

	int NoEvaluations = _ev->getNoEval();

	// send No of eval
	MPI_Send(&NoEvaluations ,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value

	// send time spent
    }
    else
    {
	int k = 0;
	int NumTask = _sc->getNumTask();
      for( k=1; k< NumTask; k++ )
      {
	int value1;
	short field1;
	unsigned char direction1;
	int type1;
	MPI_Status stat;
	int NoEvaluations;

        // Receive moves
	MPI_Recv(&value1, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);
	MPI_Recv(&field1, 1, MPI_SHORT, k,10,MPI_COMM_WORLD, &stat); 
	MPI_Recv(&direction1, 1, MPI_UNSIGNED_CHAR, k,10,MPI_COMM_WORLD, &stat);
	MPI_Recv(&type1, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);	
	MPI_Recv(&NoEvaluations, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);

	_ev->addNoEval( NoEvaluations );

        // Choose the best move
	if( CurrentColor == _board->color1 ) // O, black
	{
	  if( value1 > value )
	  {
	     Move::MoveType t = (Move::MoveType)type1;
	     Move m1( field1, direction1, t );
	     foundBestMove(0, m1, value1);
	  }
	}
	else // X, red
	{
	  if( value1 < value )
	  {
	     Move::MoveType t = (Move::MoveType)type1;
	     Move m1( field1, direction1, t );
	     foundBestMove(0, m1, value1);
	  }
	}
      }
    }*/
}

/*
int MiniMaxStrategy::minimax( int depth)
{
	//cout << "Color:" << _board->actColor() << endl;
	int value =0;
	int localdepth=_maxDepth;
	int CurrentColor = _board->actColor();
	
	int min = 40000;
	int max = -40000;
	
	Move m;
    MoveList list;
	
	//if( depth >= _maxDepth )
	if( depth >= localdepth )
	{	
		value = evaluate();
        // Increment no of evaluations
        _ev->incrementNoEval();

		//cout << "MaxDepth: " << localdepth << "  Aborting at this depth, value = " << value << endl;
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
    else
    {
    	return min;
    }
}

void MiniMaxStrategy::searchBestMove()
{
	if( print )
	cout << "\nIn searchBestMove()\n";

	int value = 0;
	long min = 20000;
	long max = -20000;
	int currentDepth = 0;
	int NoMoves= 0;
	int CurrentColor = _board->actColor(); // stored own color
	
    Move m;
    MoveList list;

    // generate list of allowed moves, put them into <list>
    generateMoves(list);
    
    // loop over all moves
    while(list.getNext(m)) 
    {
    	if( print )
    	cout << "Move: " << m.name();
 
    	NoMoves++;
        if( NoMoves % _sc->getNumTask() == _sc->getRank() )	
        {	   

//i	cout << "Rank: " << _sc->getRank() << " running Move: " << NoMoves << endl;
        // Below code executes only for specific ranks

    	// PLAY MOVE
    	_board->playMove(m);
    	
    	value = minimax(currentDepth);
    	
    	if( print )
    	cout << " Value: " << value;
    	
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
    	cout << " Max: " << max << " Min: " << min << endl;
    	
    	// TAKE BACK
    	_board->takeBack();
    	
    	// Finishing the evaluation of tree node
    	finishedNode(0,&m);
	}
	else
	{
	    	finishedNode(0,&m);
	}
    }
   
    if( print )
    cout << "Color:" << CurrentColor << " No of Moves: " << NoMoves << endl;

    if( print )
    cout << "Rank:" << ((CurrentColor==_board->color1)?"O":"X")<<_sc->getRank() << " No of Moves: " << NoMoves << endl;

    // print individual best moves
    cout << "Best move for rank: " << ((CurrentColor==_board->color1)?"O":"X") << \
    _sc->getRank() << " Move: " << _bestMove.name() << " Value:" << _bestMove.value << endl;

    // Send receive best moves
    if( _sc->getRank() != 0 )
    {
        // send best move
	MPI_Send(&_bestMove.value,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value
	MPI_Send(&_bestMove.field,1, MPI_SHORT,0, 10, MPI_COMM_WORLD ); //field
	MPI_Send(&_bestMove.direction,1, MPI_UNSIGNED_CHAR,0, 10, MPI_COMM_WORLD );// direction
	MPI_Send(&_bestMove.type,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // TYPE-INT

	int NoEvaluations = _ev->getNoEval();

	// send No of eval
	MPI_Send(&NoEvaluations ,1, MPI_INT,0, 10, MPI_COMM_WORLD ); // value

	// send time spent
    }
    else
    {
	int k = 0;
	int NumTask = _sc->getNumTask();
      for( k=1; k< NumTask; k++ )
      {
	int value1;
	short field1;
	unsigned char direction1;
	int type1;
	MPI_Status stat;
	int NoEvaluations;

        // Receive moves
	MPI_Recv(&value1, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);
	MPI_Recv(&field1, 1, MPI_SHORT, k,10,MPI_COMM_WORLD, &stat); 
	MPI_Recv(&direction1, 1, MPI_UNSIGNED_CHAR, k,10,MPI_COMM_WORLD, &stat);
	MPI_Recv(&type1, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);	
	MPI_Recv(&NoEvaluations, 1, MPI_INT, k,10,MPI_COMM_WORLD, &stat);

	_ev->addNoEval( NoEvaluations );

        // Choose the best move
	if( CurrentColor == _board->color1 ) // O, black
	{
	  if( value1 > value )
	  {
	     Move::MoveType t = (Move::MoveType)type1;
	     Move m1( field1, direction1, t );
	     foundBestMove(0, m1, value1);
	  }
	}
	else // X, red
	{
	  if( value1 < value )
	  {
	     Move::MoveType t = (Move::MoveType)type1;
	     Move m1( field1, direction1, t );
	     foundBestMove(0, m1, value1);
	  }
	}
      }
    }
   
        
    	
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
	
	
//    // we try to maximize bestEvaluation
//    int bestEval = minEvaluation();
//    int eval;
//
//    Move m;
//    MoveList list;
//
//    // generate list of allowed moves, put them into <list>
//    generateMoves(list);
//
//    // loop over all moves
//    while(list.getNext(m)) 
//    {
//
//		// draw move, evalute, and restore position
//		playMove(m);
//		eval = evaluate();
//		takeBack();
//		
//		if (eval > bestEval) 
//		{
//		    bestEval = eval;
//		    foundBestMove(0, m, eval);
//		}
//    }
//
//    finishedNode(0,0);
}
*/
// register ourselve as a search strategy
MiniMaxStrategy miniMaxStrategy;
