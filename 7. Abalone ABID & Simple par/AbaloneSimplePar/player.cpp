/**
 * Computer player
 *
 * (1) Connects to a game communication channel,
 * (2) Waits for a game position requiring to draw a move,
 * (3) Does a best move search, and broadcasts the resulting position,
 *     Jump to (2)
 *
 * (C) 2005, Josef Weidendorfer
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <mpi.h>

#include "board.h"
#include "search.h"
#include "eval.h"
#include "network.h"

/* Global, static vars */
NetworkLoop l;
Board b;
Evaluator ev;

/* Which color to play? */
int myColor = Board::color1;

/* Which search strategy to use? */
int strategyNo = 0;

/* Max search depth */
int maxDepth = 0;

/* Maximal number of moves before terminating (negative for infinity) */
int maxMoves = -1;

/* to set verbosity of NetworkLoop implementation */
extern int verbose;

/* remote channel */
char* host = 0;       /* not used on default */
int rport = 23412;

/* local channel */
int lport = 23412;

/* change evaluation after move? */
bool changeEval = true;


/**
 * MyDomain
 *
 * Class for communication handling for player:
 * - start search for best move if a position is received
 *   in which this player is about to draw
 */
class MyDomain: public NetworkDomain
{
public:
    MyDomain(int p) : NetworkDomain(p) {}

    void sendBoard();

protected:
    void received(char* str);
};

void MyDomain::sendBoard()
{
    static char tmp[500];
    sprintf(tmp, "pos %s\n", b.getState());
    if (verbose) printf(tmp+4);
    broadcast(tmp);
}

void MyDomain::received(char* str)
{
    static double Tottime=0;
    long TotNoEval=0;

    if (strncmp(str, "quit", 4)==0) {
	l.exit();
	return;
    }

    if (strncmp(str, "pos ", 4)!=0) return;

    b.setState(str+4);
    if (verbose) {
	printf("\n\n==========================================\n");
	printf(str+4);
    }

    int state = b.validState();
    int len;
    int NumTask, i;
    MPI_Status stat;
    MPI_Comm_size( MPI_COMM_WORLD, &NumTask);

    // Sum up evaluations
    TotNoEval =  ev.getNoEval();

    if ((state != Board::valid1) && 
	(state != Board::valid2)) {
	printf("%s\n", Board::stateDescription(state));

	switch(state) 
    {
	    case Board::timeout1:
	    case Board::timeout2:
	    case Board::win1:
	    case Board::win2:
    
            char boardLayout[500];
            // Send Exit Signal to slave through EXIT board
            len = sprintf(boardLayout, "exit%s\n",b.getState());
            for( i = 1; i< NumTask; i++)
            {
                MPI_Send( boardLayout, 500, MPI_CHAR, i, 10, MPI_COMM_WORLD );
                long NoEval;
                MPI_Recv(&NoEval, 1, MPI_LONG, i, 10, MPI_COMM_WORLD, &stat);
                printf("\n Evaluations for Rank: %d Evals: %ld", i, NoEval);
                TotNoEval += NoEval;
            }

        	printf("\nTotal run time for player ");
	    	printf("%s ", (myColor == Board::color1) ? "O":"X");
        	printf(" =  %lf\n", (Tottime/1000));
	        printf("Evaluations per sec = %lf \n", (TotNoEval*1000/Tottime));
            printf("Total Evaluations= %ld\n", TotNoEval);
		    l.exit();
	    default:
		break;
	}
	return;
    }

    if (b.actColor() & myColor) {
	struct timeval t1, t2;

	gettimeofday(&t1,0);
	Move m = b.bestMove();
	gettimeofday(&t2,0);

	int msecsPassed =
	    (1000* t2.tv_sec + t2.tv_usec / 1000) -
	    (1000* t1.tv_sec + t1.tv_usec / 1000);

	printf("%s ", (myColor == Board::color1) ? "O":"X");
	if (m.type == Move::none) {
	    printf(" can not draw any move ?! Sorry.\n");
	    return;
	}
	printf("draws '%s' (after %d.%03d secs)...\n",
	       m.name(), msecsPassed/1000, msecsPassed%1000);
	    Tottime += msecsPassed;

	printf("\nMoving board for Move:%s \n ", m.name());
	b.playMove(m, msecsPassed);
	sendBoard();

	if (changeEval)
	    ev.changeEvaluation();
    
    // Sum up evaluations
    TotNoEval =  ev.getNoEval();

	/* stop player at win position */
	int state = b.validState();
	if ((state != Board::valid1) && 
	    (state != Board::valid2)) {
	    printf("%s\n", Board::stateDescription(state));
	    switch(state) {
		case Board::timeout1:
		case Board::timeout2:
		case Board::win1:
		case Board::win2:

            char boardLayout[500];
            // Send Exit Signal to slave through EXIT board
            len = sprintf(boardLayout, "exit%s\n",b.getState());
            for( i = 1; i< NumTask; i++)
            {
                MPI_Send( boardLayout, 500, MPI_CHAR, i, 10, MPI_COMM_WORLD );
                long NoEval;
                MPI_Recv(&NoEval, 1, MPI_LONG, i, 10, MPI_COMM_WORLD, &stat);
                printf("\n Evaluations for Rank: %d Evals: %ld", i, NoEval);
                TotNoEval += NoEval;
            }

	        printf("\nTotal run time for player ");
    		printf("%s ", (myColor == Board::color1) ? "O":"X");
	        printf(" =  %lf\n", (Tottime/1000));
        	printf("Evaluations per sec = %lf \n", (TotNoEval*1000/Tottime));
            printf("Total Evaluations= %ld\n", TotNoEval);
		    l.exit();
		default:
		    break;
	    }
	}

		maxMoves--;
		if (maxMoves == 0) 
		{
            char boardLayout[500];
            // Send Exit Signal to slave through EXIT board
            len = sprintf(boardLayout, "exit%s\n",b.getState());
            for( i = 1; i< NumTask; i++)
            {
                MPI_Send( boardLayout, 500, MPI_CHAR, i, 10, MPI_COMM_WORLD );
                long NoEval;
                MPI_Recv(&NoEval, 1, MPI_LONG, i, 10, MPI_COMM_WORLD, &stat);
                printf("\n Evaluations for Rank: %d Evals: %ld", i, NoEval);
                TotNoEval += NoEval;
            }

			//printf("\n %s", str);
		    printf("Terminating because given number of moves drawn.\n");
		    broadcast("quit\n");
		    l.exit();
		}
    }    
}



/*
 * Main program
 */

void printHelp(char* prg, bool printHeader)
{
    if (printHeader)
	printf("Computer player V 0.1 - (C) 2005 Josef Weidendorfer\n"
	       "Search for a move on receiving a position in which we are expected to draw.\n\n");

    printf("Usage: %s [options] [X|O] [<strength>]\n\n"
	   "  X                Play side X\n"
	   "  O                Play side O (default)\n"
	   "  <strength>       Playing strength, depending on strategy\n"
	   "                   A time limit can reduce this\n\n" ,
	   prg);
    printf(" Options:\n"
	   "  -h / --help      Print this help text\n"
	   "  -v / -vv         Be verbose / more verbose\n"
	   "  -s <strategy>    Number of strategy to use for computer (see below)\n"
	   "  -n               Do not change evaluation function after own moves\n"
	   " -<integer>        Maximal number of moves before terminating\n"
	   "  -p [host:][port] Connection to broadcast channel\n"
	   "                   (default: 23412)\n\n");

    printf(" Available search strategies for option '-s':\n");

    char** strs = SearchStrategy::strategies();
    for(int i = 0; strs[i]; i++)
	printf("  %2d : Strategy '%s'%s\n", i, strs[i],
	       (i==strategyNo) ? " (default)":"");
    printf("\n");

    exit(1);
}

void parseArgs(int argc, char* argv[])
{
    int arg=0;
    while(arg+1<argc) {
	arg++;
	if (strcmp(argv[arg],"-h")==0 ||
	    strcmp(argv[arg],"--help")==0) printHelp(argv[0], true);
	if (strncmp(argv[arg],"-v",2)==0) {   
	    verbose = 1;
	    while(argv[arg][verbose+1] == 'v') verbose++;
	    continue;
	}
	if (strcmp(argv[arg],"-n")==0)	{
	    changeEval = false;
	    continue;
	}
	if ((strcmp(argv[arg],"-s")==0) && (arg+1<argc)) {
	    arg++;
	    if (argv[arg][0]>='0' && argv[arg][0]<='9')
               strategyNo = argv[arg][0] - '0';
            continue;
        }

	if ((argv[arg][0] == '-') &&
	    (argv[arg][1] >= '0') &&
	    (argv[arg][1] <= '9')) {
	    int pos = 2;

	    maxMoves = argv[arg][1] - '0';
	    while((argv[arg][pos] >= '0') &&
		  (argv[arg][pos] <= '9')) {
		maxMoves = maxMoves * 10 + argv[arg][pos] - '0';
		pos++;
	    }
	    continue;
	}

	if ((strcmp(argv[arg],"-p")==0) && (arg+1<argc)) {
	    arg++;

	    if (argv[arg][0]>'0' && argv[arg][0]<='9') {
		lport = atoi(argv[arg]);
		continue;
	    }
	    char* c = strrchr(argv[arg],':');
	    int p = 0;
	    if (c != 0) {
		*c = 0;
		p = atoi(c+1);
	    }
	    host = argv[arg];
	    if (p) rport = p;
	    continue;
	}
	
	if (argv[arg][0] == 'X') {
	    myColor = Board::color2;
	    continue; 
	}
	if (argv[arg][0] == 'O') {
	    myColor = Board::color1;
	    continue;
	}

	int strength = atoi(argv[arg]);
	if (strength == 0) {
	    printf("ERROR - Unknown option %s\n", argv[arg]);
	    printHelp(argv[0], false);
	}
	
	maxDepth = strength;
    }
}

int main(int argc, char* argv[])
{
    parseArgs(argc, argv);
    int NumTask, rank;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &NumTask);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);

    SearchStrategy* ss = SearchStrategy::create(strategyNo);

    if (verbose)
   	printf("Using strategy '%s' ...\n", ss->name());

    ss->setMaxDepth(maxDepth);

    b.setSearchStrategy( ss );
    ss->setEvaluator(&ev);

    if( verbose )
    printf("\nIn main, rank = %d\n", rank);

    ss->registerCallbacks(new SearchCallbacks(verbose, NumTask, rank));

    // Only rank 0 opens port & communicates
    if( rank == 0 )
    {
        printf("\nStarting, rank = %d\n", rank);
        MyDomain d(lport);
        l.install(&d);
        if (host) d.addConnection(host, rport);

        //l.install(&d);
	    l.run();
    }
    else // others wait for signal from Master - rank 0
    {
    	printf("\nRank %d executing wait\n", rank);
        ss->setBoard( &b );
	    ss->wait();
    }

    MPI_Finalize();

    printf("\n Rank %d exiting from Main \n", rank );
}
