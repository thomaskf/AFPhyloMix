//
//  definitions.h
//
//  Created by Thomas Wong on 06/11/19.
//  Copyright © 2019 Thomas Wong. All rights reserved.
//

#ifndef definitions_h
#define definitions_h


// Which method? 1 - MAX or 2 - SUM?
// regarding the location-mutation variable
#define METHOD_EQU 1
#define MIN_FREQ_DIFF 0.0001
#define NUM_DIGITS 8 // for number display
#define DECIMAL_PLACE 5

//==========================
// prior parameters:
// gamma distribution for frequencies
#define FREQ_SHAPE 2
#define FREQ_RATE 0.1
// beta distribution for edge length
#define EDGE_BETA 0.005
// window sizes for move
#define FREQ_WINSIZE 0.3
#define EDGE_WINSIZE 0.01
#define ERR_WINSIZE 0.001
//==========================


//==========================
// for reads:
#define COVER_THRES 50 // the minimum coverage for potential snp positions
#define COVER_RATIO_THRES 0.2 // the minimum coverage for potential snp positions >= 0.2 * average-coverage
#define WINDOW_SIZE 50

#define LOGLIKE_THRES 20.0
#define FREQ_DIFF_THRES 0.02
#define INIT_SEQ_ERR 0.005 // initialize value of sequencing error is 0.5%
#define MAX_SEQ_ERR 0.01 // assume the maximum value of sequencing error is 1%

#define MAX_GAP_RATE 0.01 // if the column has more than 1% gap, ignore that column

#define K_FACTOR 1000

#define PAIR_COVER_RATIO_THRES 0.2 // the minimum coverage for a pair of positions
#define BAD_MIN_RATIO 0.01
#define badPairRatioTHRES 0.00
#define BAD_MIN_PAIR 2

#define TOO_LONG_TRIM 20 // alert if it has been trimmed for more than 20 bases (in any end), for identifying the "dropping" region
#define MINDROPRATIO 0.01
#define MINDROPFREQ 7
#define TOO_CLOSE_DROP 30
#define MAX_DROP_SIZE 500

// for optimzation
#define DELTA 1e-5
#define MIN_VALUE 0.01
#define MAXIT 5
#define CONVERAGE_THRES 0.1

#define MIN_COV_RATIO 0.5 // 0.5 of the average

// default method to estimate the back-mutation positions
#define DEFAULT_BACK_MUTATE_METHOD 1

#define MAPQ_THRES 20 // the alignments with mapq < MAPQ_THRES are ignored

#define QV_THRES 25 // min quality value of a base to consider

// to show the positions
// #define SHOWPOS

// to list all bad mutation positions
// #define LISTBADPOS
//==========================


#endif
