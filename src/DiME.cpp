//
//  main.cpp
//  AdaptedKLMoving
//
//  Created by Leciel on 13-3-19.
//  Copyright (c) 2013年 MSc Computer Science. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <complex>
//#include <omp.h>
//using namespace std;
using namespace Rcpp;

#define SEED (unsigned)time(NULL)+getpid()



    double randf() {
        return static_cast<double>((double)rand()/RAND_MAX);
    }
    int randBin() {
        return rand() % 2;
    }
    
    int randi(int lower, int upper) {
        if (upper > lower) {
            return rand() % (upper - lower) + lower;
        }
        else {
            return lower;
        }
    }
    
    int min(double *list, int length) {
        int min = 0;
        for (int i = 0; i < length; i++) {
            if (list[i] < list[min]) {
                min = i;
            }
        }
        return min;
    }
    
    int max(double *list, int length) {
        int max = 0;
        for (int i = 0; i < length; i++) {
            if (list[i] > list[max]) {
                max = i;
            }
        }
        return max;
    }
    
    int average(double *data, int length) {
        double sum = 0;
        for (int i = 0; i < length; i++) {
            sum += data[i];
        }
        return sum/length;
    }
    
    double stddev(double *data, int length) {
        double sd = 0;
        double mean = average(data, length);
        for (int i = 0; i < length; i++) {
            sd += (data[i] - mean) * (data[i] - mean);
        }
        sd = sqrt(sd / (length - 1));
        return sd;
    }
    
    int * randperm(int size) {
        int *a = new int[size];
        for (int k = 0; k != size; k++)
            a[k] = k;
        for (int k = size-1; k != 0; k--) {
            int j = randi(0, k+1);
            int temp = a[j];
            a[j] = a[k];
            a[k] = temp;
        }
        return a;
    }
    
    void popInit(int **pop, int popsize, int dim) {
        for (int i = 0; i < popsize; i++) {
            for (int j = 0; j < dim; j++) {
                pop[i][j] = randBin();
            }
        }
    }
    
    void structuredPopInit(int **pop, int popsize, int dim) {
        for (int i = 0; i < popsize; i++) {
            double prob = 0.5*(double)(i+1)/(double)popsize;
            for (int j = 0; j < dim; j++) {
                if (randf() < prob) {
                    pop[i][j] = 1;
                }
                else {
                    pop[i][j] = 0;
                }
            }
        }
    }
    
    void EDPopInit(int **pop, int popsize, double *estPropDist, int dim) {
        for (int i = 0; i != popsize; i++) {
            for (int j = 0; j != dim; j++) {
                if (randf() < estPropDist[j]) {
                    pop[i][j] = 1;
                }
                else {
                    pop[i][j] = 0;
                }
            }
        }
    }

    void EDPopInit2(int **pop, int popsize, double prob, int dim) {
        for (int i = 0; i != popsize; i++) {
            for (int j = 0; j != dim; j++) {
                if (randf() < prob) {
                    pop[i][j] = 1;
                }
                else {
                    pop[i][j] = 0;
                }
            }
        }
    }
    
    double modularity(double **adj_matrix, int *config, int dim, double *param, std::vector<int> selected, std::vector<int> unselected) {
        double eval = 0;
        int sum = 0;
        for (int i = 0; i < dim; i++) {
            sum += config[i];
        }
        double S = sum;
        double Sc = dim - sum;
        //        vector<int> selected;
        //        vector<int> unselected;
        //        //int *selected = new int[sum];
        //        //int *unselected = new int[dim-sum];
        //        //int count_sel = 0;
        //        //int count_unsel = 0;
        //        for (int i = 0; i < dim; i++) {
        //            if (config[i]) {
        //                selected.push_back(i);
        //                //selected[count_sel] = i;
        //                //count_sel++;
        //            }
        //            else {
        //                unselected.push_back(i);
        //                //unselected[count_unsel] = i;
        //                //count_unsel++;
        //            }
        //        }
        double OS = 0;
        double BS = 0;
        for (int m = 0; m < sum; m++) {
            for (int n = 0; n < sum; n++) {
                OS = OS + adj_matrix[selected[m]][selected[n]];
            }
        }
        for (int m = 0; m < sum; m++) {
            for (int n = 0; n < dim-sum; n++) {
                BS = BS + adj_matrix[selected[m]][unselected[n]];
            }
        }
        eval = S * Sc * (OS / (S * S) - BS / (S * Sc));
        param[0] = S;
        param[1] = OS;
        //delete [] selected;
        //delete [] unselected;
        return eval;
    }
    
    double deltaW(double **adj_matrix, int *config, int dim, int which, double *param, int flag) {
        double dW = 0;
        double S = param[0];
        double OS = param[1];
        double N = (double)dim;
        double sum1 = 0;
        double sum2 = 0;
        for (int j = 0; j < dim; j++) {
            sum1 += adj_matrix[which][j] * config[j];
            sum2 += adj_matrix[which][j];
        }
        sum1 = sum1 * 2;
        if (flag) {
            dW = OS * N / (S * (S - 1)) - 1 * N / (S - 1) * sum1 + sum2;
        }
        else {
            dW = - OS * N / (S * (S + 1)) + 1 * N / (S + 1) * sum1 - sum2;
        }
        return dW;
    }
    
    double localMoving(double **adj_matrix, std::vector< std::vector<int> > const& adj_list, int *community, int dim, int verbose) {
        //int *bestRes = new int[dim];
        double param[2];
        std::vector<int> selected;
        std::vector<int> unselected;
        for (int i = 0; i < dim; i++) {
            if (community[i]) {
                selected.push_back(i);
            }
            else {
                unselected.push_back(i);
            }
        }
        double fitness = modularity(adj_matrix, community, dim, param, selected, unselected);
        double bestfitness = fitness;
        double lastPeak = -10000;
        //int *order;
        int count = 0;
        std::vector<int> sel;
        std::vector<int> unsel;
        for (int i = 0; i < dim; i++) {
            if (community[i]) {
                sel.push_back(i);
            }
            else {
                unsel.push_back(i);
            }
        }
        while (1) {
            for (int j = 0; j != dim; j++) {
                double dW = deltaW(adj_matrix, community, dim, j, param, community[j]);
                if (dW > 0) {
                    community[j] = 1 - community[j];
                    if (community[j]) {
                        sel.push_back(j);
                        param[0] += 1;
                        double deltaOS = 0;
                        for (int k = 0; k != adj_list[j].size(); k++) {
                            if (community[adj_list[j][k]]) {
                                deltaOS += 2;
                            }
                        }
                        param[1] += deltaOS;
                    }
                    else {
                        int toErase = 0;
                        for (int l = 0; l < sel.size(); l++) {
                            if (sel[l] == j) {
                                toErase = l;
                                break;
                            }
                        }
                        sel.erase(sel.begin()+toErase);
                        unsel.push_back(j);
                        param[0] -= 1;
                        double deltaOS = 0;
                        for (int k = 0; k != adj_list[j].size(); k++) {
                            if (community[adj_list[j][k]]) {
                                deltaOS += 2;
                            }
                        }
                        param[1] -= deltaOS;
                    }
                    fitness += dW;
                    bestfitness = fitness;
                }
            }
            int sum = 0;
            for (int j = 0; j != dim; j++) {
                sum += community[j];
            }
            if (verbose) {
                //cout<<"Current best W is:\t"<<lastPeak<<"\t\tCommunity size is:\t"<<sum<<"       \r"<<flush;
            }
            if (bestfitness > lastPeak) {
                lastPeak = bestfitness;
                count = 0;
            }
            else {
                count++;
            }
            if (count == 1) {
                if (verbose) {
                    std::cout<<"\n";
                }
                break;
            }
        }
        return lastPeak;
    }
   
    // [[Rcpp::export]]
    List commextr(int d, int ps, NumericVector matrix, NumericVector res, int vb)
    {
        int verbose = vb;
        int DIM = d;
        int POPSIZE = ps;
        
        //omp_set_num_threads(*ncores);
        
        //srand(SEED);
   
        //cout<<"\nLoading adjacency matrix...";
        /* Load adjacency matrix */
        double **adj_matrix = new double*[DIM];
        for (int i = 0; i != DIM; i++) {
            adj_matrix[i] = new double[DIM];
        }
        for (unsigned int i = 0; i != DIM; i++) {
            for (unsigned int j = 0; j != DIM; j++) {
                adj_matrix[i][j] = matrix[i+DIM*j];
            }
        }
        //cout<<"Done.\n";
        
        //cout<<"\nBuilding adjacency list...";
        /* Construct adjacency list using 2-D vector */
        std::vector< std::vector<int> > adj_list;
        for (int i = 0; i != DIM; i++) {
            std::vector<int> ithEdges;
            for (int j = 0; j != DIM; j++) {
                if (adj_matrix[i][j]) {
                    ithEdges.push_back(j);
                }
            }
            adj_list.push_back(ithEdges);
        }
        //cout<<"Done.\n";
        
        //        for (int i = 0; i != DIM; i++) {
        //            cout<<i<<"\t";
        //            for (int j = 0; j < adj_list[i].size(); j++) {
        //                cout<<adj_list[i][j]<<"\t";
        //            }
        //            cout<<"\n";
        //        }
        
        /* Initialize a population for estimation of the probability distribution for community membership. */
        int initPopSize = POPSIZE; //Test using 1/10 of desired population size.
        int **estPop = new int*[initPopSize];
        for (int i = 0; i != initPopSize; i++) {
            estPop[i] = new int[DIM];
        }
        //popInit(estPop, initPopSize, DIM);
        structuredPopInit(estPop, initPopSize, DIM);
        
        /* Perform local greedy search on estimator population, similar to EDA. */
        double sum = 0;
        if (verbose) {
            std::cout<<"\nEstimating probability distribution of optimal solution(s)...\n";
        }
        //cout<<"0% Done.\r"<<flush;
        for (int i = 0; i < initPopSize; i++) {
            //int nThreads = omp_get_num_threads();
            //cout<<"Number of threads is "<<nThreads<<"\n";
            //cout<<i<<endl;
            if (verbose) {
                //cout<<"Optimizing "<<i+1<<" out of "<<POPSIZE/10<<" solutions...\n";
            }
            //double result = localMoving(adj_matrix, adj_list, estPop[i], DIM, verbose);
        }
        for (int i = 0; i != initPopSize; i++) {
            for (int j = 0; j != DIM; j++) {
                sum += estPop[i][j];
            }
        }
        if (verbose) {
            std::cout<<"\n\n";
        }
        //double averageCommSize = sum / (initPopSize*DIM);
        double *estProbDist = new double[DIM];
        double freq;
        for (int j = 0; j != DIM; j++) {
            freq = 0;
            for (int i = 0; i != initPopSize; i++) {
                freq += estPop[i][j];
            }
            //estProbDist[j] = freq / initPopSize * averageCommSize;
            estProbDist[j] = freq / sum;
        }
        //        double top = estProbDist[max(estProbDist, DIM)];
        //        for (int j = 0; j != DIM; j++) {
        //            estProbDist[j] += 1 - top; // Scale probabilities to maximum 1.
        //        }
        
        if (verbose) {
            std::cout<<"\nInitializing population with randomly seeded solutions...";
        }
        /* Initialize a population with "seeds" for community membership. */
        int newPopsize = POPSIZE * 10;
        int **pop = new int*[newPopsize];
        for (int i = 0; i < newPopsize; i++) {
            pop[i] = new int[DIM];
        }
        //structuredPopInit(pop, POPSIZE, DIM);
        EDPopInit(pop, newPopsize, estProbDist, DIM);
        //EDPopInit2(pop, POPSIZE, averageCommSize, DIM);
        
        //double **fits = new double*[newPopsize];
//        for (int i = 0; i < POPSIZE; i++) {
//            fits[i] = new double[2];
//        }
//        double **update = new double*[POPSIZE];
//        for (int i = 0; i < POPSIZE; i++) {
//            update[i] = new double[2];
//        }
        
        double *fitList = new double[newPopsize];
        if (verbose) {
            std::cout<<"Done.\n";
            std::cout<<"\nNow starting community extraction...\n";
        }
        /* Perform local on population of randomly seeded solutions and determine best result. */
        for (int i = 0; i < newPopsize; i++) {
            if (verbose) {
                std::cout<<"Now evolving individual no. "<<i+1<<"\n";
            }
            fitList[i] = localMoving(adj_matrix, adj_list, pop[i], DIM, verbose);
            //fitList[i] = adaptedKLMoving(adj_matrix, pop[i], DIM, TOL, 1);
        }
        
//        if (*followup) {
//            vector<int> selected;
//            vector<int> unselected;
//            if (verbose) {
//                cout<<"\nNow performing follow-up rounds of optimization using greedy search...\n";
//            }
//            /* Follow-up rounds of optimization by greedy search */
//            for (int i = 0; i != POPSIZE; i++) {
//                for (int j = 0; j != DIM; j++) {
//                    if (pop[i][j]) {
//                        selected.push_back(j);
//                        //selected[count_sel] = i;
//                        //count_sel++;
//                    }
//                    else {
//                        unselected.push_back(j);
//                        //unselected[count_unsel] = i;
//                        //count_unsel++;
//                    }
//                }
//                modularity(adj_matrix, pop[i], DIM, fits[i], selected, unselected);
//                selected.clear();
//                unselected.clear();
//            }
//            int count = 0;
//            int nfail = 0;
//            double lastPeak = -10000;
//            while (nfail < 1) {
//                //            int bestind = max(fitList,POPSIZE);
//                //            for (int i = 0; i < DIM; i++) {
//                //                res[i] = pop[bestind][i];
//                //            }
//                for (int i = 0; i < POPSIZE; i++) {
//                    selected.clear();
//                    unselected.clear();
//                    update[i][0] = 0;
//                    update[i][1] = deltaW(adj_matrix, pop[i], DIM, 0, fits[i], pop[i][0]);
//                    for (int j = 1; j < DIM; j++) {
//                        double now = deltaW(adj_matrix, pop[i], DIM, j, fits[i], pop[i][j]);
//                        if (now > update[i][1]) {
//                            update[i][1] = now;
//                            update[i][0] = j;
//                        }
//                    }
//                    if (update[i][1] > 0) {
//                        pop[i][(int)update[i][0]] = 1 - pop[i][(int)update[i][0]];
//                        for (int j = 0; j != DIM; j++) {
//                            if (pop[i][j]) {
//                                selected.push_back(j);
//                                //selected[count_sel] = i;
//                                //count_sel++;
//                            }
//                            else {
//                                unselected.push_back(j);
//                                //unselected[count_unsel] = i;
//                                //count_unsel++;
//                            }
//                        }
//                        fitList[i] = modularity(adj_matrix, pop[i], DIM, fits[i], selected, unselected);
//                    }
//                    //cout<<fits[i]<<"\t";
//                }
//                //            if (fits[max(fits,POPSIZE)] <= fits[bestind]) {
//                //                break;
//                //            }
//                int bestF = max(fitList, POPSIZE);
//                count++;
//                if (fitList[bestF] > lastPeak) {
//                    nfail = 0;
//                    lastPeak = fitList[bestF];
//                    for (int j = 0; j != DIM; j++) {
//                        res[j] = pop[bestF][j];
//                    }
//                }
//                else {
//                    nfail++;
//                }
//                int sum = 0;
//                for (int j = 0; j != DIM; j++) {
//                    sum += res[j];
//                }
//                //cout<<"Delta W for this round is:\t"<<update[0][1]<<"\n";
//                if (verbose) {
//                    cout<<"Iteration No.:\t"<<count<<"\n";
//                    cout<<"Current score W = "<<lastPeak<<"\n";
//                    cout<<"Optimal community size is:\t"<<sum<<"\n";
//                    cout<<"Standard deviation is:\t"<<stddev(fitList, POPSIZE)<<"\n";
//                }
//            }
//        }
        //int bestF = max(fitList,POPSIZE);
        int bestF = max(fitList, POPSIZE);
        //        for (int i = 0; i != POPSIZE; i++) {
        //            cout<<fitList[i]<<endl;
        //        }
        int commSize = 0;
        for (int j = 0; j != DIM; j++) {
            res[j] = pop[bestF][j];
            commSize += res[j];
        }
        //cout<<"Best fitness W is: "<<fitList[bestF]<<endl;
        if (verbose) {
            std::cout<<"Search process done."<<std::endl;
            std::cout<<"Best fitness W is: "<<fitList[bestF]<<std::endl;
            std::cout<<"Best community size is: "<<commSize<<std::endl;
        }
        for (int i = 0; i != DIM; i++) {
            delete [] adj_matrix[i];
        }
        delete [] adj_matrix;
        for (int i = 0; i != initPopSize; i++) {
            delete [] estPop[i];
        }
        delete [] estPop;
        delete [] estProbDist;
        for (int i = 0; i != POPSIZE; i++) {
            delete [] pop[i];
            //delete [] fits[i];
            //delete [] update[i];
        }
        delete [] pop;
        //delete [] fits;
        //delete [] update;
        delete [] fitList;
        
        return List::create(
        _["fitList"] = fitList[bestF],
        _["commSize"] = commSize,
        _["res"] = res);
    }

 