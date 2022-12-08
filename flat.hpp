#pragma once
#define AP_INT_MAX_W 4096

#include <iostream>
#define __gmp_const const
#include "ap_int.h"
#include "ap_fixed.h"
#include "hls_stream.h"
//#include "hls_math.h"

constexpr unsigned floorlog2(unsigned x)
{
#pragma HLS INLINE
    return x == 1 ? 0 : 1+floorlog2(x >> 1);
}

constexpr unsigned ceillog2(unsigned x)
{
#pragma HLS INLINE
    return x == 1 ? 0 : floorlog2(x - 1) + 1;
}

constexpr unsigned bitsNeeded(unsigned n) {
#pragma HLS INLINE
  return n <= 1 ? 0 : 1 + bitsNeeded((n + 1) / 2);
}

template <typename T>
static constexpr T ceildiv(T dividend, T divisor)
{
#pragma HLS INLINE
    return (dividend + divisor - 1) / divisor;
}

template <typename T>
static constexpr T roundup(T dividend, T divisor)
{
#pragma HLS INLINE
    return ceildiv(dividend, divisor) * divisor;
}

template <typename T>
static constexpr T ap_fixed_epsilon()
{
#pragma HLS INLINE
    return T(1.0 / (1 << (T::width - T::iwidth)));
}

// Common Functionality Definition
#define LOG_1(n) (((n) >= 2) ? 1 : 0)
#define LOG_2(n) (((n) >= 1<<2) ? (2 + LOG_1((n)>>2)) : LOG_1(n))
#define LOG_4(n) (((n) >= 1<<4) ? (4 + LOG_2((n)>>4)) : LOG_2(n))
#define LOG_8(n) (((n) >= 1<<8) ? (8 + LOG_4((n)>>8)) : LOG_4(n))
#define LOG(n)   (((n) >= 1<<16) ? (16 + LOG_8((n)>>16)) : LOG_8(n))

//typedef ap_fixed<32,8,AP_RND_CONV,AP_SAT> data_t;
typedef ap_fixed<32,8> data_t;
//Off-chip Interface Parameters
constexpr int NUM_HP_IFC = 6;
constexpr int HP_IFC_BANDWIDTH = 128;
constexpr int MAX_OFF_CHIP_BW = NUM_HP_IFC*HP_IFC_BANDWIDTH;

//Data Width Parameters. (Same as the data_t for now)
constexpr int ATTENTION_DATAWIDTH                = 32;   
constexpr int VALUE_DATAWIDTH           = 32;                                             
constexpr int OACTS_DATAWIDTH                    = 32;      

constexpr int MAX_WT_X = 2048;
constexpr int MAX_WT_Y = 1000;
constexpr int MAX_X = 1;
constexpr int MAX_Y = 2048;
constexpr int MIN_WT_X = 2048;
constexpr int MIN_WT_Y = 1000;
constexpr int MIN_X = 1;
constexpr int MIN_Y = 2048;
// constexpr int X = 1;
// constexpr int Y = 2048;

constexpr int PARALLEL_K = 32; // Given 2048 as k-dim size, process 32 at one time
constexpr int PARALLEL_N = 20; //20 //20 to get higher parallelism
constexpr int PARALLEL_M = 1;

constexpr int PROCESSING_BUFFER_DEPTH = MAX_WT_Y / PARALLEL_N;

//On-chip buffer size
constexpr int MAX_BLOCK_COUNT = MAX_WT_X * MAX_WT_Y / PARALLEL_K / PARALLEL_N;
constexpr int MIN_BLOCK_COUNT = MIN_WT_X * MIN_WT_Y / PARALLEL_K / PARALLEL_N;
constexpr int MAX_INNERLOOP_COUNT = ceildiv(PARALLEL_K*PARALLEL_N, 3*MAX_OFF_CHIP_BW/VALUE_DATAWIDTH);
constexpr int MIN_INNERLOOP_COUNT = ceildiv(PARALLEL_K*PARALLEL_N, 3*MAX_OFF_CHIP_BW/VALUE_DATAWIDTH);

constexpr int ATTEN_BRAM_WIDTH = 32; //Stores 1 ATTEN values at one row
constexpr int MAX_ATTEN_NUM = 2048;
constexpr int MAX_ATTEN_BRAM_ROW = MAX_ATTEN_NUM / (ATTEN_BRAM_WIDTH / ATTENTION_DATAWIDTH); //32 rows as the maximum -> 4 brams in need
constexpr int VALUE_URAM_WIDTH = 4 * 72;
constexpr int READ_PARALLEL_VALUE = PARALLEL_K * PARALLEL_N /(VALUE_URAM_WIDTH / VALUE_DATAWIDTH); //32*50/9, VALUE num per block

constexpr int VALUE_NUM_PER_ROW_URAM = VALUE_URAM_WIDTH / VALUE_DATAWIDTH; //9
constexpr int BLOCK_NUM = PARALLEL_K * PARALLEL_N * PARALLEL_M;
constexpr int MAX_VALUE_URAM_ROW = (MAX_WT_X * MAX_WT_Y / BLOCK_NUM) * ceildiv(BLOCK_NUM, VALUE_NUM_PER_ROW_URAM); //block num
constexpr int BLOCK_VALUE_TOTAL_IFC = PARALLEL_K * PARALLEL_N / (HP_IFC_BANDWIDTH / VALUE_DATAWIDTH); //32 * 50 / 4 = 400 ifc blocks in total
constexpr int BLOCK_VALUE_URAM_ROW = ceildiv(PARALLEL_K * PARALLEL_N, VALUE_NUM_PER_ROW_URAM) ;//URAM rows needed for one block
constexpr int MAX_VALUE_NUM = 2048000; //2048 * 1000
//constexpr int MAX_ATTENS_NUM = 2048; //1*2048
constexpr int VALUE_PER_CYCLE = HP_IFC_BANDWIDTH / VALUE_DATAWIDTH * NUM_HP_IFC; //24 VALUEs per cycle
constexpr int ATTEN_PER_CYCLE = HP_IFC_BANDWIDTH / ATTEN_DATAWIDTH * NUM_HP_IFC; //24 VALUEs per cycle
constexpr int VALUE_PER_THREE_CYCLE = 3 * VALUE_PER_CYCLE; //72 VALUEs
constexpr int BLOCK_VALUE_PER_THREE_CYCLE = 3 * VALUE_PER_CYCLE / VALUE_NUM_PER_ROW_URAM; //8 rows


constexpr int PARALLEL_ATTEN_BANDWIDTH           = PARALLEL_K*ATTEN_DATAWIDTH;
constexpr int PARALLEL_VALUE_BANDWIDTH           = PARALLEL_K*PARALLEL_N*ATTEN_DATAWIDTH;
constexpr int MAX_IFC_ENTRY = (MAX_WT_X * MAX_WT_Y * VALUE_DATAWIDTH + ATTEN_DATAWIDTH * MAX_X * MAX_Y) / HP_IFC_BANDWIDTH;

//Systolic array parameters
constexpr int SYSTOLIC_DIM_X = PARALLEL_K;
constexpr int SYSTOLIC_DIM_Y = PARALLEL_N;
//Systolic array
constexpr int Parallel_ATTEN = 1;
constexpr int Parallel_VALUE_PE = 4;
constexpr int Total_PE = Parallel_ATTEN * Parallel_VALUE_PE;
constexpr int VALUE_ENTRY_PE = MAX_VALUE_NUM / Parallel_VALUE_PE;
constexpr int ATTEN_ENTRY_PE = MAX_ATTEN_NUM / Parallel_ATTEN;

void LINEAR(
    ap_uint<HP_IFC_BANDWIDTH>               *ifc1, //ifc1[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH>               *ifc2, //ifc2[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH>               *ifc3, //ifc3[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH>               *ifc4, //ifc4[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH>               *ifc5, //ifc5[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH>               *ifc6, //ifc6[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH>               *ifc7, //ifc7[MAX_X*MAX_WT_Y],
    int X, //Input Shape X
    int Y, //Input Shape Y
    int Wt_X, //VALUE Shape X
    int Wt_Y,  //VALUE Shape Y
    int bias
);
