//!!!!!Remember to adjust the buffer size params in the header file-> Now the maximum value is for no tiling


#include "FC_Layer.hpp"

#define LOGGING_LINEAR_FUNC 0
#define LOGGIGNG_READATTENBUFF_FUNC 0
#define LOGGING_RUNDATAFLOW_FUNC 0
#define LOGGING_DPECOMPUTATION_FUNC 0
#define LOGGING_READFROMMEM_FUNC 0 

void ReadFromMem(
    ap_uint<HP_IFC_BANDWIDTH> *ifc1,
    ap_uint<HP_IFC_BANDWIDTH> *ifc2,
    ap_uint<HP_IFC_BANDWIDTH> *ifc3,
    ap_uint<HP_IFC_BANDWIDTH> *ifc4,
    ap_uint<HP_IFC_BANDWIDTH> *ifc5,
    ap_uint<HP_IFC_BANDWIDTH> *ifc6,
    ap_uint<VALUE_URAM_WIDTH> value_buffer[MAX_VALUE_URAM_ROW],
    hls::stream<ap_uint<VALUE_DATAWIDTH> > value_stream[PARALLEL_K],
    ap_uint<ATTEN_BRAM_WIDTH> atten_buffer[MAX_ATTEN_BRAM_ROW],
    hls::stream<data_t> &atten_stream,
    int X,
    int Y,
    int Wt_X,
    int Wt_Y)
{
    int block_count = Wt_X * Wt_Y / PARALLEL_K / PARALLEL_N; //compute how many blocks the program needs to load ----> 1000*2048/32/20
    int inner_loop_count = ceildiv(PARALLEL_K*PARALLEL_N, 3*MAX_OFF_CHIP_BW/VALUE_DATAWIDTH); //640/(24*3) --> burst read: how many three cycles read needed for one block
    int max_uram_row = ceildiv(PARALLEL_K*PARALLEL_N, MAX_VALUE_URAM_ROW); // compute for each block, how many lines of data it needs from URAM

    int value_offset = 0;
    //compute what is the residual after burst read
    int residual = PARALLEL_K * PARALLEL_N - inner_loop_count*72; //72 VALUEs read in 3 cycles
#ifdef LOGGING_READFROMMEM_FUNC
    std::cout << "burst read for values!!" << std::endl;
#endif
    for (size_t i = 0; i < block_count; ++i) { //iterate all blocks
#pragma HLS loop_tripcount min = MIN_BLOCK_COUNT max = MAX_BLOCK_COUNT avg = MAX_BLOCK_COUNT  

#ifdef LOGGING_READFROMMEM_FUNC
        std::cout << "Process block: " << i << std::endl;
#endif

        for (size_t j = 0; j < inner_loop_count; ++j) { //iterate cycles needed for one block transfer
#pragma HLS loop_tripcount min = MIN_INNERLOOP_COUNT max = MAX_INNERLOOP_COUNT avg = MAX_INNERLOOP_COUNT 
#pragma HLS Pipeline II=4

            //cycle1
#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Process cycle 1: " << std::endl;
#endif
            ap_uint<MAX_OFF_CHIP_BW> payload1 = 0;
            
            int addr_offset = i * BLOCK_VALUE_TOTAL_IFC;
            int addr_offset_1 = addr_offset + 3 * j * NUM_HP_IFC;
            payload1.range(6 * HP_IFC_BANDWIDTH - 1, 5 * HP_IFC_BANDWIDTH) = ifc1[addr_offset_1 + 5];
            payload1.range(5 * HP_IFC_BANDWIDTH - 1, 4 * HP_IFC_BANDWIDTH) = ifc2[addr_offset_1 + 4];
            payload1.range(4 * HP_IFC_BANDWIDTH - 1, 3 * HP_IFC_BANDWIDTH) = ifc3[addr_offset_1 + 3];
            payload1.range(3 * HP_IFC_BANDWIDTH - 1, 2 * HP_IFC_BANDWIDTH) = ifc4[addr_offset_1 + 2];
            payload1.range(2 * HP_IFC_BANDWIDTH - 1, 1 * HP_IFC_BANDWIDTH) = ifc5[addr_offset_1 + 1];
            payload1.range(1 * HP_IFC_BANDWIDTH - 1, 0) = ifc6[addr_offset_1];
#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Payload1: " << payload1 << std::endl;
#endif
            //cycle2

#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Process cycle 2: " << std::endl;
#endif
            ap_uint<MAX_OFF_CHIP_BW> payload2 = 0;
            int addr_offset_2 = addr_offset + (3 * j + 1) * NUM_HP_IFC;

            payload2.range(6 * HP_IFC_BANDWIDTH - 1, 5 * HP_IFC_BANDWIDTH) = ifc1[addr_offset_2 + 5];
            payload2.range(5 * HP_IFC_BANDWIDTH - 1, 4 * HP_IFC_BANDWIDTH) = ifc2[addr_offset_2 + 4];
            payload2.range(4 * HP_IFC_BANDWIDTH - 1, 3 * HP_IFC_BANDWIDTH) = ifc3[addr_offset_2 + 3];
            payload2.range(3 * HP_IFC_BANDWIDTH - 1, 2 * HP_IFC_BANDWIDTH) = ifc4[addr_offset_2 + 2];
            payload2.range(3 * HP_IFC_BANDWIDTH - 1, 2 * HP_IFC_BANDWIDTH) = ifc4[addr_offset_2 + 2];
            payload2.range(2 * HP_IFC_BANDWIDTH - 1, 1 * HP_IFC_BANDWIDTH) = ifc5[addr_offset_2 + 1];
            payload2.range(1 * HP_IFC_BANDWIDTH - 1, 0) = ifc6[addr_offset_2];

            // cycle 3
            ap_uint<MAX_OFF_CHIP_BW> payload3 = 0;
            int addr_offset_3 = addr_offset+(3 * j + 2) * NUM_HP_IFC;

#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Payload2: " << payload2 << std::endl;
#endif

#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Process cycle 3: " << std::endl;
#endif
            payload3.range(6 * HP_IFC_BANDWIDTH - 1, 5 * HP_IFC_BANDWIDTH) = ifc1[addr_offset_3 + 5];
            payload3.range(5 * HP_IFC_BANDWIDTH - 1, 4 * HP_IFC_BANDWIDTH) = ifc2[addr_offset_3 + 4];
            payload3.range(4 * HP_IFC_BANDWIDTH - 1, 3 * HP_IFC_BANDWIDTH) = ifc3[addr_offset_3 + 3];
            payload3.range(3 * HP_IFC_BANDWIDTH - 1, 2 * HP_IFC_BANDWIDTH) = ifc4[addr_offset_3 + 2];
            payload3.range(2 * HP_IFC_BANDWIDTH - 1, 1 * HP_IFC_BANDWIDTH) = ifc5[addr_offset_3 + 1];
            payload3.range(1 * HP_IFC_BANDWIDTH - 1, 0) = ifc6[addr_offset_3];

#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Payload3: " << payload3 << std::endl;
#endif
            //VALUE_BUFFER load data
#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Load data to buffer: " << std::endl;
#endif
            int addr_offset_4 = BLOCK_VALUE_URAM_ROW * i + j * BLOCK_VALUE_PER_THREE_CYCLE;
            ap_uint<3 *MAX_OFF_CHIP_BW> combine_payload = 0;
            combine_payload.range(MAX_OFF_CHIP_BW - 1, 0) = payload1;
            combine_payload.range(2 * MAX_OFF_CHIP_BW - 1, MAX_OFF_CHIP_BW) = payload2;
            combine_payload.range(3 * MAX_OFF_CHIP_BW - 1, 2 * MAX_OFF_CHIP_BW) = payload3;

            //ap_uint<MAX_OFF_CHIP_BW> combine_payload = 0;
            //std::cout << "check VALUE matrix" << std::endl;

#ifdef LOGGING_READFROMMEM_FUNC
            std::cout << "Check combine payload info: " << std::endl;
#endif
            for (int idx_count = 0; idx_count < VALUE_PER_THREE_CYCLE; idx_count+=9) {
                ap_uint<VALUE_URAM_WIDTH> value_load = 0;
                int offset = idx_count + VALUE_NUM_PER_ROW_URAM;
                for (int i = 0; i < VALUE_NUM_PER_ROW_URAM; ++i)
                {
                    value_load.range((i + 1) * VALUE_DATAWIDTH - 1, i * VALUE_DATAWIDTH) 
                        = combine_payload.range((idx_count + i + 1) * VALUE_DATAWIDTH - 1, (idx_count + i) * VALUE_DATAWIDTH);
#ifdef LOGGING_READFROMMEM_FUNC
                    std::cout << combine_payload.range((idx_count + i + 1) * VALUE_DATAWIDTH - 1, (idx_count + i) * VALUE_DATAWIDTH) << ", ";
#endif
                }
#ifdef LOGGING_READFROMMEM_FUNC
                std::cout << std::endl;
#endif
                value_buffer[addr_offset_4++] = value_load.range(VALUE_URAM_WIDTH - 1, 0);
            }
        }

#ifdef LOGGING_READFROMMEM_FUNC
        std::cout << "Process the remaining data for the block read" << std::endl;
#endif

        //process residual data
        //compute how many cycles needed
        int count = 0;
        int full_cycles = ceildiv(residual, VALUE_PER_CYCLE);
        ap_uint<3*MAX_OFF_CHIP_BW> payload = 0;
        for (size_t ii = 0; ii < full_cycles; ++ii) {
            ap_uint<MAX_OFF_CHIP_BW> temp = 0;
            int addr_offset_1 = BLOCK_VALUE_URAM_ROW * i + VALUE_PER_THREE_CYCLE * inner_loop_count + ii * VALUE_PER_CYCLE;
            for (int jj = 0; jj < 6; ++jj) {
                if (count >= residual) break;
                temp.range((jj+1) * HP_IFC_BANDWIDTH - 1, jj * HP_IFC_BANDWIDTH) = ifc1[addr_offset_1 + jj];
                count++;
            }
            payload.range((ii+1)*MAX_OFF_CHIP_BW-1, ii*MAX_OFF_CHIP_BW) = temp;
        }
        int value_rows = ceildiv(residual, VALUE_NUM_PER_ROW_URAM);
        int offset = VALUE_DATAWIDTH * VALUE_NUM_PER_ROW_URAM;
        int addr_offset = BLOCK_VALUE_URAM_ROW * i + inner_loop_count * BLOCK_VALUE_PER_THREE_CYCLE;
        for (int ii = 0; ii < value_rows; ++ii) {
        #pragma HLS UNROLL
            ap_uint<VALUE_URAM_WIDTH> value_load = 0;
            if (ii != value_rows - 1) {
                value_load = payload.range((ii+1)*offset-1, ii*offset); 
            }
            else {
                int w_offset = residual - ii * VALUE_NUM_PER_ROW_URAM;
                value_load.range(w_offset*VALUE_DATAWIDTH-1, 0) = payload.range(ii*offset+w_offset*VALUE_DATAWIDTH-1, ii*offset);
            }
            value_buffer[addr_offset++] = value_load.range(VALUE_URAM_WIDTH - 1, 0);
        }
        //std::cout << addr_offset << "Address Offset!!!!!" << std::endl;
    }
    
#ifdef LOGGING_READFROMMEM_FUNC
    std::cout << "Check VALUE buffer!!" << std::endl;
    int count = 0;
    for (int i = 0; i < MAX_VALUE_URAM_ROW; i+=72) {
        std::cout << "Process block: " << i/72 << std::endl;
        for (int j = 0; j < 72; ++j) {
            ap_uint<VALUE_URAM_WIDTH> load = value_buffer[i + j];
            for (int k = 0; k < 9; ++k) {
                data_t wt = 0;
                wt.range(VALUE_DATAWIDTH-1, 0) = load.range((k+1)*VALUE_DATAWIDTH-1, k*VALUE_DATAWIDTH);
                std::cout << wt << ", ";
            }
        }
        std::cout << std::endl << "********" << std::endl;
    }

    std::cout << "start reading attentions!!" << std::endl;
#endif
    int atten_count = ceildiv(X*Y, ATTEN_PER_CYCLE)- 1;
    residual = X*Y - atten_count * ATTEN_PER_CYCLE;
    for (int i = 0; i < atten_count; ++i) {
        ap_uint<MAX_OFF_CHIP_BW> payload = 0;
        int addr_offset = BLOCK_VALUE_TOTAL_IFC + i * NUM_HP_IFC;
        for (int j = 0; j < NUM_HP_IFC; ++j) {
            payload.range((j+1) * HP_IFC_BANDWIDTH - 1, j*HP_IFC_BANDWIDTH) = ifc6[addr_offset+j];
        }
        int buffer_offset = i * ATTEN_PER_CYCLE;
        for (int j = 0; j < ATTEN_PER_CYCLE; ++j) {
            atten_buffer[buffer_offset+j] = payload.range((j+1) * ATTEN_DATAWIDTH-1, j*ATTEN_DATAWIDTH);
        }
    } 
    //process remaining data -- less than one cycle
    //std::cout << "Process atten residual elements" << std::endl;
    int remaining_cycle = residual / ATTEN_PER_CYCLE + 1;
    //assert(remaining_cycle == 1);
    for (int i = 0; i < remaining_cycle; ++i) {
        ap_uint<MAX_OFF_CHIP_BW> payload = 0;
        int addr_offset = BLOCK_VALUE_TOTAL_IFC + atten_count * NUM_HP_IFC;
        for (int j = 0; j < NUM_HP_IFC; ++j) {
        #pragma HLS UNROLL
            if (j * 4 > residual) break;
            payload.range((j+1) * HP_IFC_BANDWIDTH - 1, j*HP_IFC_BANDWIDTH) = ifc6[addr_offset+j];
        }
        int buffer_offset = atten_count * ATTEN_PER_CYCLE + i * ATTEN_PER_CYCLE;
        for (int j = 0; j < residual; ++j) {
            atten_buffer[buffer_offset+j] = payload.range((j+1) * ATTEN_DATAWIDTH-1, j*ATTEN_DATAWIDTH);
        }
    }

    //std::cout << "Wegith streaming starts" << std::endl;
    // read 640 numbers to VALUEs_stream
    int block_num_x = Wt_X / PARALLEL_K;
    int block_num_y = Wt_Y / PARALLEL_N;
    for (int block_x = 0; block_x < block_num_x; ++block_x)
    {
        for (int block_y = 0; block_y < block_num_y; ++block_y)
        {
        #pragma HLS PIPELINE 
            int base_offset = (block_x * block_num_y + block_y) * BLOCK_VALUE_URAM_ROW;
            for (int i = 0; i < PARALLEL_K; ++i) //32
            {
            #pragma HLS UNROLL
                for (int j = 0; j < PARALLEL_N; ++j) //20
                {
                #pragma HLS UNROLL
                    // int idx_x = (block_x * PARALLEL_K * PARALLEL_N + block_y * PARALLEL_N + i * PARALLEL_N + j) / VALUE_NUM_PER_ROW_URAM;
                    // int idx_y = (block_x * PARALLEL_K * PARALLEL_N + block_y * PARALLEL_N + i * PARALLEL_N + j) % VALUE_NUM_PER_ROW_URAM;
                    int idx_x = (i * PARALLEL_N + j) / VALUE_NUM_PER_ROW_URAM;  //remove 
                    int idx_y = (i * PARALLEL_N + j) % VALUE_NUM_PER_ROW_URAM;
                    ap_uint<VALUE_NUM_PER_ROW_URAM*VALUE_DATAWIDTH> load = value_buffer[base_offset+idx_x];
                    value_stream[i].write(load.range((idx_y + 1) * VALUE_DATAWIDTH - 1, (idx_y) * VALUE_DATAWIDTH));
                    //std::cout << base_offset+idx_x << "Address Offset!!!!!" << std::endl;
                    //std::cout << VALUEs_stream[i] << ", ";
                }
                //std::cout << std::endl;
            }
            //std::cout << base_offset << "Address Offset!!!!!" << std::endl;
        }
    }
    //std::cout << "VALUE written: " << count << std::endl;

    //std::cout << "ATTENs streaming starts" << std::endl;
    //stream the ATTEN values
    for (int i = 0; i < X*Y; ++i) {
        data_t tmp = 0;
        tmp.range(ATTEN_DATAWIDTH-1, 0) = atten_buffer[i].range(ATTEN_BRAM_WIDTH-1, 0);
        //std::cout << "ATTEN_BUFFER: " << ATTEN_buffer[i] << ", " << "TMP: " << tmp << std::endl;
        atten_stream.write(tmp); //same size for now
    }
}

void CreateBitMask(hls::stream<ap_uint<WEIGHTS_DATAWIDTH>> weights_stream[PARALLEL_K], ap_uint<PARALLEL_N * WEIGHTS_DATAWIDTH> processing_buffer[PARALLEL_K],
            ap_uint<PARALLEL_N> bit_buffer_weights[PARALLEL_K])
{
    //std::cout << "MASK CALLED" << std::endl;
    for (int i = 0; i < PARALLEL_K; i++) {
        while (true) {
            if (!weights_stream[i].empty()) {
                break;
            }
        }
    }
//  change loop order to smooth the pipeline
    for (int i = 0; i < PARALLEL_K; ++i)
    { 
        ap_uint<PARALLEL_N*WEIGHTS_DATAWIDTH> payload = 0; //payload -> PARALLEL_N*WEIGHTS_DATA_WIDTH*(WT_Y/PARALLEL_N)
        ap_uint<PARALLEL_N> bitmask = 0;
        for (int j = 0; j < PARALLEL_N;++j) {
#pragma HLS PIPELINE II = 1
            ap_uint<WEIGHTS_DATAWIDTH> data = weights_stream[i].read();
            payload.range((j+1)*WEIGHTS_DATAWIDTH-1, j*WEIGHTS_DATAWIDTH) = data;
            bitmask.range(j, j) = data == 0? 0 : 1;
        }
        processing_buffer[i] = payload;
        bit_buffer_weights[i] = bitmask;
    }
}

void DPEUnit(data_t atten_value, int atten_idx, ap_uint<PARALLEL_N * VALUE_DATAWIDTH> processing_buffer[PARALLEL_K],
                    ap_uint<PARALLEL_N> bit_buffer_weights[PARALLEL_K], data_t output_buf[PARALLEL_K][PARALLEL_N], int k_idx, int y) {
}

//batch_num here is for recording how many groups of PARALLEL_K we have processed
void DPEComputation(int batch_num, data_t ATTEN_TEMP_BUFFER[PARALLEL_K], int block_idx_x, int block_idx_y,  ap_uint<PARALLEL_N * VALUE_DATAWIDTH> processing_buffer[PARALLEL_K],
                    ap_uint<PARALLEL_N> bit_buffer_weights[PARALLEL_K],  data_t buffer_out[MAX_WT_Y/PARALLEL_N][PARALLEL_N], int Wt_X, int Wt_Y) {
    //CALL DPEUNIT to perform computation within each row
    //SKIP Entire row in the dense matrix if the pruned attention score is zero

    //Better stream the entire processing buffer with weights corresponding to non-zero attention score

}


void OutputBuffer(
    ap_uint<HP_IFC_BANDWIDTH> *oacts_ifc,
    hls::stream<data_t> &output_stream,
    int X,
    int Wt_Y,
    int address,
    data_t output_buf[MAX_WT_Y/PARALLEL_N][PARALLEL_N])
{
    int overall_addr = address;
    int loop_count = X * Wt_Y;
    //Do I need output stream?

    for (int i = 0; i < Wt_Y * X / (HP_IFC_BANDWIDTH / OACTS_DATAWIDTH); ++i) { //1000/4 = 250
        ap_uint<HP_IFC_BANDWIDTH> payload = 0;
        for (int j = 0; j < HP_IFC_BANDWIDTH / OACTS_DATAWIDTH; ++j) {
            //payload.range((j+1) * OACTS_DATAWIDTH-1, j*OACTS_DATAWIDTH) = output_stream.read().range(OACTS_DATAWIDTH-1, 0);
            payload.range((j+1) * OACTS_DATAWIDTH-1, j*OACTS_DATAWIDTH) = output_buf[i][j].range(OACTS_DATAWIDTH-1, 0);
        }
        oacts_ifc[overall_addr++] = payload;
    }
}

inline void ReadAttenBuff(hls::stream<data_t> &atten_stream, data_t ATTEN_TEMP_BUFFER[PARALLEL_K]) {
#ifdef LOGGIGNG_READATTENBUFF_FUNC
    std::cout << "ATTEN BUFFER" << std::endl;
#endif
    for (int i = 0; i < PARALLEL_K; ++i) {
        data_t attention_score = atten_stream.read();
//do the data pruning here!!!!!!!
        if (attention_score < THRESHOLD_VALUE) {
            attention_score = 0;
        }
        ATTEN_TEMP_BUFFER[i] = attention_score;
#ifdef LOGGIGNG_READATTENBUFF_FUNC
        std::cout << ATTEN_TEMP_BUFFER[i] << ", ";
#endif

    }
}

void RunDataFlow(int block_num_x, int block_num_y, hls::stream<data_t> &atten_stream, data_t ATTEN_TEMP_BUFFER[PARALLEL_K], hls::stream<ap_uint<VALUE_DATAWIDTH>> value_stream[PARALLEL_K],
            data_t output_buf[MAX_WT_Y/PARALLEL_N][PARALLEL_N], ap_uint<PARALLEL_N * VALUE_DATAWIDTH> first_processing_buffer[PARALLEL_K], ap_uint<PARALLEL_N * VALUE_DATAWIDTH> second_processing_buffer[PARALLEL_K],
            ap_uint<PARALLEL_N> first_bit_buffer_value[PARALLEL_K], ap_uint<PARALLEL_N> second_bit_buffer_value[PARALLEL_K], int Wt_X, int Wt_Y)
{
#pragma HLS DATAFLOW
    CreateBitMask(value_stream, first_processing_buffer, first_bit_buffer_value);
    for (int i = 0; i < block_num_x; ++i) {

#ifdef LOGGING_RUNDATAFLOW_FUNC
        std::cout << "BLOCK_NUM: " << i << std::endl;
        std::cout << "Read Buffer" << std::endl;
#endif

        ReadAttenBuff(atten_stream, ATTEN_TEMP_BUFFER);

#ifdef LOGGING_RUNDATAFLOW_FUNC
        std::cout << "Double buffering" << std::endl;
#endif

        for (int j = 0; j < block_num_y; ++j) {
            int batch = i * block_num_y + j;
            if (batch != block_num_x*block_num_y - 1) {
                if (batch & 1 != 0) { //even case
                    CreateBitMask(value_stream, second_processing_buffer, second_bit_buffer_value);
                    DPEComputation(batch, ATTEN_TEMP_BUFFER, i, j, first_processing_buffer, first_bit_buffer_value, output_buf, Wt_X, Wt_Y);
                } else {
                    CreateBitMask(value_stream, first_processing_buffer, first_bit_buffer_value);
                    DPEComputation(batch, ATTEN_TEMP_BUFFER, i, j, second_processing_buffer, second_bit_buffer_value, output_buf, Wt_X, Wt_Y);
                }
            }
            else {
                DPEComputation(batch, ATTEN_TEMP_BUFFER, i, j, second_processing_buffer, second_bit_buffer_value, output_buf, Wt_X, Wt_Y); //depends on batch_num
            }
#ifdef LOGGING_RUNDATAFLOW_FUNC
            if (j == 0) std::cout << "IJ: " << i <<" , " << j << std::endl;
#endif
        }
    }
}

void LINEAR(
    ap_uint<HP_IFC_BANDWIDTH> *ifc1, //ifc1[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH> *ifc2, //ifc2[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH> *ifc3, //ifc3[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH> *ifc4, //ifc4[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH> *ifc5, //ifc5[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH> *ifc6, //ifc6[MAX_IFC_ENTRY],
    ap_uint<HP_IFC_BANDWIDTH> *ifc7, //ifc7[MAX_X*MAX_WT_Y],
    int X, // Input Shape X
    int Y,                                               // Attention Shape Y
    int Wt_X,                                            // VALUE Shape X
    int Wt_Y,                                            // VALUE Shape Y
    int bias)
{
// #define FINAL_DIM0 X
// #define FINAL_DIM1 Y
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc1 bundle = ifc1
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc2 bundle = ifc2
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc3 bundle = ifc3
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc4 bundle = ifc4
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc5 bundle = ifc5
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc6 bundle = ifc6
#pragma HLS INTERFACE m_axi depth = MAX_IFC_ENTRY offset = slave port = ifc7 bundle = ifc1

#pragma HLS INTERFACE s_axilite port = X
#pragma HLS INTERFACE s_axilite port = Y
#pragma HLS INTERFACE s_axilite port = Wt_X
#pragma HLS INTERFACE s_axilite port = Wt_Y
#pragma HLS INTERFACE s_axilite port = bias

    //Assign URAM to VALUE buffer
    ap_uint<VALUE_URAM_WIDTH> value_buffer[MAX_VALUE_URAM_ROW]; // need 240 urams in total (layout 60 * 4)  Each row stores 9 VALUE numbers.
#pragma HLS BIND_STORAGE variable = value_buffer type = ram_t2p impl = uram latency = 1
#pragma HLS array_partition variable = value_buffer type = cyclic factor = READ_PARALLEL_VALUE dim = 0

    ap_uint<PARALLEL_N * VALUE_DATAWIDTH> first_processing_buffer[PARALLEL_K];
    ap_uint<PARALLEL_N * VALUE_DATAWIDTH> second_processing_buffer[PARALLEL_K];
#pragma HLS BIND_STORAGE variable = first_processing_buffer type = ram_t2p impl = bram latency = 1  // 4 brams -> 128 * 128 * 1 + 32 * 512 * 4
#pragma HLS BIND_STORAGE variable = second_processing_buffer type = ram_t2p impl = bram latency = 1 // 160 * 20 as one blocks

    ap_uint<PARALLEL_N> first_bit_buffer_VALUE[PARALLEL_K];
    ap_uint<PARALLEL_N> second_bit_buffer_VALUE[PARALLEL_K];


    ap_uint<ATTEN_BRAM_WIDTH> atten_buffer[MAX_ATTEN_BRAM_ROW]; // need 4 brams (8 * 2k) in total
#pragma HLS BIND_STORAGE variable = atten_buffer type = ram_t2p impl = bram latency = 1
//#pragma HLS array_partition variable = ATTEN_buffer type = cyclic factor = PARALLEL_K dim = 1 // read 32 elements at one time

    hls::stream<data_t> atten_stream;
#pragma HLS STREAM variable = atten_stream depth = PARALLEL_K type = fifo

    hls::stream<ap_uint<VALUE_DATAWIDTH>> value_stream[PARALLEL_K];
#pragma HLS STREAM variable = value_stream depth = 10*PARALLEL_N type = fifo

    hls::stream<data_t> output_stream;
#pragma HLS STREAM variable = output_stream depth = PARALLEL_N type = fifo

    data_t output_buf[MAX_WT_Y/PARALLEL_N][PARALLEL_N];
    // ap_uint<OACTS_DATAWIDTH> output_buf[MAX_WT_Y/PARALLEL_N][PARALLEL_N];
#pragma HLS BIND_STORAGE variable = output_buf type = ram_t2p impl = bram latency = 1
#pragma HLS array_partition variable = output_buf type = complete dim = 1 // read 32 elements at one time

    int block_num_x = Wt_X / PARALLEL_K;
    int block_num_y = Wt_Y / PARALLEL_N;
    data_t ATTEN_TEMP_BUFFER[PARALLEL_K];

    int address_ifc = block_num_x * block_num_y * BLOCK_VALUE_TOTAL_IFC;
#pragma HLS DATAFLOW

#ifdef LOGGING_LINEAR_FUNC
    std::cout << "Read from memory" << std::endl;
#endif

    ReadFromMem(ifc1, ifc2, ifc3, ifc4, ifc5, ifc6, value_buffer, value_stream, atten_buffer, atten_stream, X, Y, Wt_X, Wt_Y);
    //compute_systolic(ATTENs_stream, VALUEs_stream, bias, output_stream, X, Y, Wt_X, Wt_Y);
#ifdef LOGGING_LINEAR_FUNC
    std::cout << "Create Bitmask" << std::endl;
#endif
    RunDataFlow(block_num_x, block_num_y, atten_stream, ATTEN_TEMP_BUFFER, value_stream, output_buf, 
                 first_processing_buffer, second_processing_buffer, first_bit_buffer_value, second_bit_buffer_value, Wt_X, Wt_Y);
#ifdef LOGGING_LINEAR_FUNC
    std::cout << "Output Buffer" << std::endl;
    for (int i = 0; i < 50; ++i) {
        for (int j = 0; j < 20; ++j) {
            std::cout << output_buf[i][j] << ", ";
        }
        std::cout << std::endl;
    }
#endif
    //OutputBuffer(ifc7, output_stream, X, Wt_Y, address_ifc, output_buf);
    OutputBuffer(ifc7, output_stream, X, Wt_Y, 0, output_buf);
}