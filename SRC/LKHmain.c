#include "LKH.h"
#include "Genetic.h"
/*
 * This file contains the main function of the program.
 */

//_declspec(dllexport) char* LKHmain(char* tspstr);

char* LKHmain(char* tspstr)
{
    GainType Cost, OldOptimum;
    double Time, LastTime;
    Node* N;
    Candidate* NN;
    static char candstr[1000000];
    memset(candstr,'\0',sizeof(candstr));

    ParameterFileName = "ss.par";
    StartTime = LastTime = GetTime();
    MaxMatrixDimension = 20000;
    MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
        MergeWithTourGPX2;
    ReadParameters();
    ReadProblem(tspstr);
    AllocateStructures();
    CreateCandidateSet();
    N = FirstNode;
    do {
        sprintf(candstr, "%s%d ", candstr, N->Id);
        if (N->CandidateSet)
            for (NN = N->CandidateSet; NN->To; NN++) {
		sprintf(candstr, "%s%d %lf ", candstr, NN->To->Id, (double)NN->Alpha / exmaxalpha);
		//sprintf(candstr, "%s%d %lf ", candstr, NN->To->Id, (double)NN->Alpha);
            }
        sprintf(candstr, "%s\n", candstr);
    } while ((N = N->Suc) != FirstNode);
    return candstr;
}

//int main() {
//    char* str = "TYPE : TSP\nDIMENSION: 10\nEDGE_WEIGHT_TYPE : EUC_2D\nNODE_COORD_SECTION\n1 95.63049481616359 10.660647847542737\n2 50.34023017015349 90.17911792314325\n3 83.01872387549503 20.721711683501397\n4 18.885147515632795 2.3877887454321067\n5 83.77048074520069 61.738320663186805\n6 42.031268435960726 91.43603678260021\n7 22.1884703488065 94.99523324286547\n8 68.20061249862988 91.75924210717254\n9 54.5666925016468 57.34941351098437\n10 70.83363024682373 27.44594876761376\nEOF";
//    printf("%s",LKHmain(str));
//    system("pause");
//}
