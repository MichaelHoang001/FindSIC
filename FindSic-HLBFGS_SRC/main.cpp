#include <iostream>
#include <time.h>
#include "randnum.h"
#include "func.h"
#include "misc.h"
#include "HLBFGS.h"

typedef double fnum;

static int MAX_DIMENSION = 200;

funcSICff<fnum> sicfun;

funcZSICff<fnum> zsicfun;

void evalfunc(int N, double* x, double *prev_x, double* f, double* g)
{
    sicfun.getFuncGradC(x, f, g);
}

void evalfunc2(int N, double* x, double *prev_x, double* f, double* g)
{
    zsicfun.getFuncGradC(x, f, g);
}

void newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	// std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
}

int main(int argc, char** argv) {
    
    unsigned help = 0; 

    std::string FidFilename = "sicfid";
    unsigned StartDim = 3;
    fnum ETol(1E-6);
    char randtype = '0';
    std::string RandFilename = "SobolDefault";
    unsigned RandSkip = 0;
    unsigned HessianCorrections = 7;
    fnum p1 = 1E-4;
    fnum p2 = 1E-16;
    fnum p3 = 0.9;
    fnum p4 = 1E-20;
    fnum p5 = 1E+20;
    fnum p6 = 0.0;
    fnum p7 = 1E-5;
    unsigned p8 = 20; 
    unsigned p9 = 1;
    unsigned p10 = 100000;
    unsigned p11 = 10;
    unsigned p12 = 15;
    unsigned Zsym = 1;
    unsigned OutputLevel = 1;
    
    int i;
    // Assign inputs parameters
    for (i = 2; i < argc; i+=2) {
        if (argv[i-1][0] == '-' && argv[i-1][2] == '\0') {
            switch (argv[i-1][1]) {
                case 'f': FidFilename = argv[i]; break;
                case 'd': StartDim = atoi(argv[i]); break;
                case 'e': ETol = atof(argv[i]); break;
                case 'r': randtype = *argv[i]; break;
                case 's': RandSkip = atoi(argv[i]);
                case 't': RandFilename = argv[i]; break;
                case 'z': Zsym = atoi(argv[i]); break;
                case 'o': OutputLevel = atoi(argv[i]); break;
                default: help = 1;
            }
        } else if (argv[i-1][0] == '-' && argv[i-1][1] == 'p' && argv[i-1][3] == '\0') {
            switch (argv[i-1][2]) {
                case '0': HessianCorrections = atoi(argv[i]); break;
                case '1': p1 = atof(argv[i]); break;
                case '2': p2 = atof(argv[i]); break;
                case '3': p3 = atof(argv[i]); break;
                case '4': p4 = atof(argv[i]); break;
                case '5': p5 = atof(argv[i]); break;
                case '6': p6 = atof(argv[i]); break;
                case '7': p7 = atof(argv[i]); break;
                case '8': p8 = atoi(argv[i]); break;
                case '9': p9 = atoi(argv[i]); if (p9>1) p9=1; if (p9<0) p9=0; break;
                default: help = 1;
            }
        } else if (argv[i-1][0] == '-' && argv[i-1][1] == 'p' && argv[i-1][2] == '1' && argv[i-1][4] == '\0') {
            switch (argv[i-1][2]) {
                case '0': p10 = atoi(argv[i]); break;
                case '1': p11 = atoi(argv[i]); break;
                case '2': p12 = atoi(argv[i]); break;
                default: help = 1;
            }
        } else help = 1;
    }
    if (help == 1 || i == argc) {
        std::cout << std::endl << "Usage: findsics [options]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -h                           Display this information" << std::endl;
        std::cout << "  -f <output_file_prefix>      File prefix to write solutions [sicfid]" << std::endl;
        std::cout << "  -d <start_dimension>         Starting dimension to begin search [2]" << std::endl;
        std::cout << "  -e <frame_tolerance>         Allowed error tolerance in frame bound [1E-6]" << std::endl;
        std::cout << "  -r <random_generator>        Pseudorandom number generator (0: Mersenne Twister, 1: Sobol sequence) [0]" << std::endl;
        std::cout << "  -s <random_seed>             Number of random number calls for initialisation [0]" << std::endl;
        std::cout << "  -t <sobol_file>              File for Sobol sequence direction numbers [SobolDefault]" << std::endl;
        std::cout << "  -z <zauner_symmetry>         Restrict search to Zauner 1-eigenspace (0: off, 1: on) [1]" << std::endl;
        std::cout << "  -o <output_level>            Level of output (0: none, 1: each new solution, 2: each trial) [1]" << std::endl;
        std::cout << "  -p0 <p0>                     Number of LBFGS Hessian corrections (typically 3-20) [7]" << std::endl;
        std::cout << "  -p1 <p1>                     Function tolerance used in line-search [1E-4]" << std::endl;
        std::cout << "  -p2 <p2>                     Variable tolerance used in line-search [1E-16]" << std::endl;
        std::cout << "  -p3 <p3>                     Gradient tolerance used in line-search [0.9]" << std::endl;
        std::cout << "  -p4 <p4>                     Minimum step used in line-search [1E-20]" << std::endl;
        std::cout << "  -p5 <p5>                     Maximum step used in line-search [1E+20]" << std::endl;
        std::cout << "  -p6 <p6>                     Stop criterion ||G||/max(1,||X||) < p6 [0]" << std::endl;
        std::cout << "  -p7 <p7>                     Stop criterion ||G|| < p7 [1E-5]" << std::endl;
        std::cout << "  -p8 <p8>                     Maximum number of evaluations in line search [20]" << std::endl; 
        std::cout << "  -p9 <p9>                     L-BFGS strategy (0: standard, 1: M1QN3) [1]" << std::endl;
        std::cout << "  -p10 <p10>                   Maximum number of iterations [100000]" << std::endl;
        std::cout << "  -p11 <p11>                   Hessian update interval (typically 0-200) [10]" << std::endl;
        std::cout << "  -p12 <p12>                   ICFS parameter [15]" << std::endl;      
        return 0;        
    } else if (OutputLevel > 0) {
        std::cout << std::endl << "Options:" << std::endl;
        std::cout << "  File prefix to write solutions = " << FidFilename << std::endl;
        std::cout << "  Starting dimension to begin search = " << StartDim << std::endl;
        std::cout << "  Allowed error tolerance in frame bound = " << ETol << std::endl;
        std::cout << "  Pseudorandom number generator (0: Mersenne Twister, 1: Sobol sequence) = " << randtype << std::endl;
        std::cout << "  Number of random number calls for initialisation = " << RandSkip << std::endl;
        std::cout << "  File for Sobol sequence direction numbers = " << RandFilename << std::endl;
        std::cout << "  Restrict search to Zauner 1-eigenspace (0: off, 1: on) = " << Zsym << std::endl;
        std::cout << "  Level of output (0: none, 1: each new solution, 2: each trial) = " << OutputLevel << std::endl;
        std::cout << "  Number of LBFGS Hessian corrections (typically 3-20) = " << HessianCorrections << std::endl;
        std::cout << "  Function tolerance used in line-search = " << p1 << std::endl;
        std::cout << "  Variable tolerance used in line-search = " << p2 << std::endl;
        std::cout << "  Gradient tolerance used in line-search = " << p3 << std::endl;
        std::cout << "  Minimum step used in line-search = " << p4 << std::endl;
        std::cout << "  Maximum step used in line-search = " << p5 << std::endl;
        std::cout << "  Stop criterion ||G||/max(1,||X||) < " << p6 << std::endl;
        std::cout << "  Stop criterion ||G|| < " << p7 << std::endl;
        std::cout << "  Maximum number of evaluations in line search = " << p8 << std::endl; 
        std::cout << "  L-BFGS strategy (0: standard, 1: M1QN3) = " << p9 << std::endl;
        std::cout << "  Maximum number of iterations = " << p10 << std::endl;
        std::cout << "  Hessian update interval (typically 0-200) = " << p11 << std::endl;
        std::cout << "  ICFS parameter = " << p12 << std::endl;      
    }
   
    std::ofstream sicfidfile;
    sicfidfile.precision(16);
    sicfidfile << std::scientific;
    
    std::ostringstream sicfidfilename;
   
    randBase<fnum>* Prandnum;
    switch (randtype) {
        case '0': Prandnum = new randMT<fnum>(1); break;
        case '1': Prandnum = new randSobol<fnum>(1,RandFilename); break;
        default: Prandnum = new randDefault<fnum>(1);
    }
    for (unsigned k = 0; k < RandSkip; ++k) {Prandnum->getNext();}
    
    std::cout.precision(16);
    std::cout << std::scientific;

    double parameter[20];
    int info[20];
    //initialize
    INIT_HLBFGS(parameter, info);
       
    parameter[0] = p1; // function tolerance used in line-search
    parameter[1] = p2; // variable tolerance used in line-search
    parameter[2] = p3; // gradient tolerance used in line-search
    parameter[3] = p4; // stpmin used in line-search
    parameter[4] = p5; // stpmax used in line-search
    parameter[5] = p6; // the stop criterion ( ||G||/max(1,||X||) < parameter[5] )
    parameter[6] = p7; // the stop criterion ( ||G|| < parameter[6] )

    info[0] = p8; // the max number of evaluation in line-search
    info[1] = 0; // the total number of evalfunc calls
    info[2] = 0; // the current number of iterations
    info[3] = p9; // The lbfgs strategy. 0: standard, 1: M1QN3 strategy[8](recommended).
    info[4] = p10; // the max number of iterations
    info[5] = 0; // 1: print message, 0: do nothing
    info[6] = p11; // T: the update interval of Hessian. (typical choices: 0-200)
    info[7] = 0; // 0: without hessian, 1: with accurate hessian
    info[8] = p12; // icfs parameter
    info[9] = 0; // 0: classical line-search; 1: modified line-search (it is not useful in practice)
    info[10] = 0; // 0: Disable preconditioned CG; 1: Enable preconditioned CG
    info[11] = 1; // 0 or 1 defines different methods for choosing beta in CG.
    info[12] = 1; // internal usage. 0: only update the diag in USER_DEFINED_HLBFGS_UPDATE_H; 1: default.

    unsigned d, N, trials, iter;
    fnum E;
    std::vector<fnum> x, r;
    std::vector<fnum>::iterator px;
    clock_t timestart;  // create timer
    clock_t timecur;
    double sortArray[MAX_DIMENSION*2];
    
    for (d=StartDim;d<=MAX_DIMENSION;d++) {
      timestart = clock(); // timer start
      
      if (OutputLevel == 2) std::cout << "d = " << d << "\t  trials = " << 0 << std::flush; // -o 2, ready rolling console output
      
      sicfidfilename.str("");
      sicfidfilename << FidFilename << "_" << d << ".txt"; // name the file
      if (!fileexists((sicfidfilename.str()).c_str())) { // file doesn't exist, run findsic 
	
	N = 2*d;
	sicfun.setNumVars(N);
	zsicfun.setNumVars(N);
        
	x.resize(N);
	r.resize(N-1);
	Prandnum->setDim(N-1);
	E = 1;
	trials = 0;
	
	while (E > ETol && !fileexists((sicfidfilename.str()).c_str())) { // while error > tolerance and no other findsic suceeded
	  
	  Prandnum->getNextD(r);
	  rnumstox(x, r); // looks like x is the array of complex numbers

	  // use HLBFGS to process x
	  if (Zsym == 0) {
	    HLBFGS(N, HessianCorrections, &x[0], evalfunc, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info); 
	  }
	  else {
	    HLBFGS(N, HessianCorrections, &x[0], evalfunc2, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info); 
	    projectZx(x);
	  } 

	  // run x through sic function to gather E, error grading
	  sicfun.getFunc(x, E);
	  ++trials;
	  
	  if (OutputLevel == 2){
	    timecur = (clock() - timestart) / (CLOCKS_PER_SEC/1000);
	    std::cout << "\r" << "d = " << d << "\t  trials = " << trials << "\tE = " << E << "\t";
	    if (timecur < 1000*10) std::cout << "ms = " << timecur << std::flush; //ms until 10 seconds
	    else if (timecur < 1000*60*2) std::cout << "sec = " << timecur/1000 << std::flush; //seconds until 2 minutes
	    else if (timecur < 1000*60*60*2) std::cout << "mins = " << timecur/(1000*60) << std::flush; //minutes until 2 hours
	    else std::cout << "hours = " << timecur/(1000*60*60) << std::flush; //hours
	  }
	}
	
	timecur = (clock() - timestart) / (CLOCKS_PER_SEC/1000); // record time, then convert to milliseconds
	if (!fileexists((sicfidfilename.str()).c_str())) { // race condition prevention
	  sicfidfile.open((sicfidfilename.str()).c_str());
          sicfidfile << "Elapsed time:\t" << timecur << "ms\r\n\tApprox\t";

	  if (timecur < 1000*10) sicfidfile << timecur << "ms\r\n"; //ms until 10 seconds
	  else if (timecur < 1000*60*2) sicfidfile << timecur/1000 << "seconds\r\n"; //seconds until 2 minutes
	  else if (timecur < 1000*60*60*2) sicfidfile << timecur/(1000*60) << "minutes\r\n"; //minutes until 2 hours
	  else sicfidfile << timecur/(1000*60*60) << "hours\r\n"; //hours

	  sicfidfile << "Elased trials:\t" << trials << "\r\n";
	  sicfidfile << "E: \t\t" << E << "\r\n";
	  sicfidfile << "Random Seed:\t" << RandSkip << "\r\n";
	  sicfidfile << "\r\nREAL\t\t\tIMAGINARY\r\n";
	  iter = 0;
	  for (px = x.begin(); px != x.end(); ++px)
	    sortArray[iter++] = *px;
	  for (iter = 0; iter < d; iter++)
	    sicfidfile << sortArray[iter] << "\t" << sortArray[iter+d] << "\r\n";
	  sicfidfile.close();
          
	  if (OutputLevel == 1){
	    std::cout << "d = " << d << "\t  trials = " << trials << "\tE = " << E << "\t";
	    if (timecur < 1000*10) std::cout << "ms = " << timecur << std::flush; //ms until 10 seconds
	    else if (timecur < 1000*60*2) std::cout << "sec = " << timecur/1000 << std::flush; //seconds until 2 minutes
	    else if (timecur < 1000*60*60*2) std::cout << "mins = " << timecur/(1000*60) << std::flush; //minutes until 2 hours
	    else std::cout << "hours = " << timecur/(1000*60*60) << std::flush; //hours    	    std::cout << "  \tms/trial = " << timecur/trials << std::flush;
	  }
	  
	} else {// another instance of findsic beat us, skip
	  
	  if (OutputLevel == 1) std::cout << "d = " << d << "\t  trials = " << trials << std::flush;
	  else if (OutputLevel == 2) std::cout << "\r" << "d = " << d << "\t  trials = " << trials << std::flush;
	  
	}
        
      } else { // file exists, sic already found, skip
	
	if (OutputLevel == 1) std::cout << "d = " << d << "\t  trials = " << 0 << std::flush;
        
      }
      
      if (OutputLevel > 0) std::cout << std::endl;
      
    }

    delete Prandnum;

    return 0; 
}

/*
#include <ctime>

funcSIC<fnum> sicfun;
funcSICf<fnum> sicfunf;
funcSICff<fnum> sicfunff;

int main(int argc, char** argv) {
   
    unsigned dd = 100;
    unsigned NN = 2*dd;
    sicfun.setNumVars(NN);
    sicfunf.setNumVars(NN);
    sicfunff.setNumVars(NN);
    fnum ff1, ff2, ff3;
    std::vector<fnum> xx, rr, gg1, gg2, gg3;
    xx.resize(NN);
    gg1.resize(NN);gg2.resize(NN);gg3.resize(NN);
    rr.resize(NN-1);
    randMT<fnum> randnum(1);
    randnum.setDim(NN-1);
    clock_t t;
    int t1(0), t2(0), t3(0);
    for (unsigned m=0; m<10000; ++m) {
        randnum.getNextD(rr);
        rnumstox(xx, rr);

        t = clock();
        sicfun.getFuncGradC(&xx[0], &ff1, &gg1[0]);
        t1 += clock()-t;
        std::cout << double(t1)/double(m+1) << std::endl;
        
        t = clock();
        sicfunf.getFuncGradC(&xx[0], &ff2, &gg2[0]);
        t2 += clock()-t;
        std::cout << double(t2)/double(m+1) << std::endl;
        
        t = clock();
        sicfunff.getFuncGradC(&xx[0], &ff3, &gg3[0]);
        t3 += clock()-t;
        std::cout << double(t3)/double(m+1) << std::endl;
        
        std::cout << std::endl << ff1-ff2 << " " << ff1-ff3 << " " << normxmy(gg1,gg2) << " " << normxmy(gg1,gg3) << std::endl << std::endl;
    
   
    }
    
    return 0;
    
 }
 */
