#include "flat.hpp"
#include "sstream"
#include <iostream>
#include <fstream>

#include <cmath>

#define LOG_READ_DATA 0
#define DEBUG_MAIN 0
#define LOOP_TILING_PARTS 4

constexpr int LOOP_TILE_X = LOOP_TILING_PARTS / 2;
constexpr int LOOP_TILE_Y = LOOP_TILING_PARTS / 2;

using namespace std;

float attention_input[MAX_X*MAX_Y];
float value_input[MAX_WT_X*MAX_WT_Y];
float output[MAX_X*MAX_WT_Y];

data_t attention[MAX_X*MAX_Y];
data_t value[MAX_WT_X*MAX_WT_Y];
data_t oacts[MAX_X*MAX_WT_Y];

void read_bin_files(int X, int Y, int Wt_X, int Wt_Y){
    std::ifstream attention_file("/nethome/zchen752/temp/attentions.csv", ios::in);
    std::ifstream value_file("/nethome/zchen752/temp/values.csv", ios::in);
    std::ifstream reference_oacts_file("/nethome/zchen752/temp/output.csv", ios::in);
    std::cout <<"Import data" <<std::endl;
    std::cout << "attention data: " << std::endl;
    for(int i = 0; i < X; ++i){
        for (int j = 0; j < Y; ++j){
            std::string temp_data;
            std::stringstream ss;
            std::getline(attention_file, temp_data, ',');
            attention_input[i*Y+j] = std::stof(temp_data);
#ifdef LOG_READ_DATA
            std::cout<<attention_input[i*Y+j]<<", ";
#endif
        }
    }
    std::cout << std::endl;
    
    std::cout << "value data: " << std::endl;
    for(int i = 0; i < Wt_X; i++){ //1000 or 2048?
        for (int j = 0; j < Wt_Y; j++) { //2048
            std::string temp_data;
            std::stringstream ss;
            std::getline(value_file, temp_data, ',');
#ifdef LOG_READ_DATA
            value_input[i*1000+j] = std::stof(temp_data);
#endif
        }
    }
    std::cout << std::endl;

    std::cout << "oacts data: " << std::endl;
    for(int i = 0; i < X; i++){
        for (int j = 0; j < Wt_Y; j++) {
            std::string temp_data;
            std::stringstream ss;
            std::getline(reference_oacts_file, temp_data, ',');
            output[i*Wt_Y+j] = std::stof(temp_data);
#ifdef LOG_READ_DATA
           std::cout<<output[i*Wt_Y+j]<<", ";
#endif
        }
    }

    std::cout << std::endl;
}

void convert_data_type(int X, int Y, int Wt_X, int Wt_Y){
    std::cout << "Attention_score" << std::endl;
    for (int i = 0; i < X; ++i) {
        //std::cout << std::endl;
        for (int j = 0; j < Y; ++j) {
            //std::cout << attention_input[i*Y+j] << ", ";
            attention[i*Y+j] = (data_t)attention_input[i*Y+j];
        }
    }

    std::cout << "Value Matrix" << std::endl;
    for (int i = 0; i < Wt_X; ++i) {
        //std::cout << std::endl;
        for (int j = 0; j < Wt_Y; ++j) {
            value[i*Wt_Y+j] = (data_t)value_input[i*Wt_Y+j];
            //std::cout << value[i*Wt_Y+j] << ", ";
        }
    }

    std::cout << "Output" << std::endl;
    for (int i = 0; i < X; ++i) {
        std::cout << std::endl;
        for (int j = 0; j < Wt_Y; ++j) {
            //std::cout << output[i*Wt_Y+j] << ", ";
            oacts[i*Wt_Y+j] = (data_t)output[i*Wt_Y+j];
        }
    }
}

int main()
{
//Set input dim to 1*2048 & 2048*1000 for now
    int X = 1;
    int Y = 2048;
    int Wt_X = 2048;
    int Wt_Y = 1000;
    std::cout <<"Program Starts!!!!!" <<std::endl;
    read_bin_files(X, Y, Wt_X, Wt_Y);  //readin input files
    convert_data_type(X, Y, Wt_X, Wt_Y);  //convert data type 
    
    int overall_addr  = 0; //set overall output address

    std::cout << "Start Data Memory Layout" << std::endl;

    uint64_t value_matrix_size_bits = VALUE_DATAWIDTH * Wt_X * Wt_Y;
    uint64_t attention_matrix_size_bits = ATTENTION_DATAWIDTH * X * Y;
    uint64_t total_num_ifc_entries = (weight_matrix_size_bits + iact_matrix_size_bits) / HP_IFC_BANDWIDTH;
    
    assert(HP_IFC_BANDWIDTH % value_DATAWIDTH == 0); //check if datawidth is divisible by the IFC_Bandwidth    assert(HP_IFC_BANDWIDTH % attention_DATAWIDTH == 0);
    assert(total_num_ifc_entries == MAX_IFC_ENTRY);

    ap_uint<HP_IFC_BANDWIDTH> ifc[MAX_IFC_ENTRY];
    

    //weight data should be read column by column plus tiling in column dimension (2048 elements per column needs to be read; depends on how many columns get read at a time)
    int value_parallel_transmission = HP_IFC_BANDWIDTH / VALUE_DATAWIDTH; //4
    int value_loop_count_X = Wt_X / PARALLEL_K; //2048 / 32 = 64
    int value_loop_count_Y = Wt_Y / PARALLEL_N; //1000 / 20 = 50
    int value_loop_residual = Wt_X * Wt_Y % value_parallel_transmission; 

    int value_complete_loop_count = value_loop_count_X * value_loop_count_Y;
    int value_tile_num_x = value_loop_count_X / LOOP_TILE_X;
    int value_tile_num_y = value_loop_count_Y / LOOP_TILE_Y;
    int atten_parallel_transmission = HP_IFC_BANDWIDTH / ATTENTION_DATAWIDTH; //128 / 32 = 4
    int atten_complete_loop_count = X * Y / atten_parallel_transmission; //2048 / 4
    int atten_loop_residual = X * Y % atten_parallel_transmission; 
    int atten_tile_num = iact_complete_loop_count / LOOP_TILING_PARTS; // 512 / 4 = 128 transmission per tile

    data_t output[X][Wt_Y];
    ap_uint<HP_IFC_BANDWIDTH> output_ifc[MAX_WT_Y];
    int output_address = 0;
    int output_parallel_transmission = HP_IFC_BANDWIDTH / OUTPUT_DATAWIDTH; //4
    int output_tile_num = (X*Wt_Y) / output_parallel_transmission / LOOP_TILING_PARTS; // 512 / 4 = 128 transmission per tile

    assert(iact_loop_residual == 0);
    assert(value_loop_residual == 0);
    //Read 32 * 20 / 4 = 160 block a time --> tile the loop to 8 parts
    int count = 0;
    std::cout << "start transferring value" << std::endl;
    //Perform Loop tiling
    for (int part_x = 0; part_x < LOOP_TILE_X; ++part_x) {
        for (int part_y = 0; part_y < LOOP_TILE_Y; ++part_y) {
            int value_addr_offset = 0;
            for (int i = part_x*value_tile_num_x; i < (part_x+1)*value_tile_num_x-1; i++) { // 64 /2 = 32
                for (int j = part_y*value_tile_num_y; j < (part_y+1)*value_tile_num_y; j++) { // 50 /2 = 25 
                    for (int row = 0; row < PARALLEL_K; ++row) { //32 
                        for (int col = 0; col < PARALLEL_N / value_parallel_transmission; ++col) { //20 /4 = 5
                            ap_uint<HP_IFC_BANDWIDTH> temp = 0;
                            int idx_x = i*PARALLEL_K+row;
                            int idx_y = j*PARALLEL_N+col*value_parallel_transmission+idx;
                            temp.range((idx+1)*VALUE_DATAWIDTH-1, idx*VALUE_DATAWIDTH) = value[idx_x*Wt_Y+idx_y].range(31, 0);

#ifdef DEBUG_MAIN
                            std::cout << "idx_x: " << idx_x << " idx_y" << idx_y << " temp:" << temp << std::endl;
                            std::cout << value[idx_x*Wt_Y+idx_y] << ", " << std::endl;
#endif

                        }
                        ifc[overall_addr + count] = temp.range(HP_IFC_BANDWIDTH-1, 0);
                        count++;
                        value_addr_offset++;

#ifdef DEBUG_MAIN
                        std::cout << "ifc: " << ifc[overall_addr + count] << std::endl;
                        std::cout << "Count" << count << std::endl;
#endif
                    }
                }
            }
            overall_addr += value_addr_offset;

            int atten_addr_offset = 0;
            int atten_tile_idx = part_x*LOOP_TILE_Y + part_y;
            for (int i = atten_tile_idx*atten_tile_num; i < (atten_tile_idx+1)*atten_tile_num-1; i++) {
                ap_uint<HP_IFC_BANDWIDTH> temp;
                for (int j = 0; j < atten_parallel_transmission; ++j) {
                    int idx = i*atten_parallel_transmission+j;
                    temp.range((j+1)*ATTENTION_DATAWIDTH-1, j*ATTENTION_DATAWIDTH) = attention[idx].range(31, 0); //QUESTION: BIG_ENDIAN OR SMALL_ENGIAN READ
                    //std::cout << "idx: " << idx << " temp:" << temp << std::endl;
                }
                ifc[overall_addr + i] = temp;
                atten_addr_offset++;
            }
            overall_addr += atten_addr_offset;

            std::cout << "Start processing linear layer" << std::endl;
            LINEAR(ifc, ifc, ifc, ifc, ifc, ifc, output_ifc, X, Y, Wt_X, Wt_Y, 0);

            int output_tile_idx = part_x*LOOP_TILE_Y + part_y;
            std::cout << "Finish running linear function on board" <<std::endl;
            for (int i = 0; i < X; ++i) {
                for (int j = output_tile_idx*output_tile_num; j < (output_tile_idx+1)*output_tile_num-1; j++) {
                output[i][j].range(OACTS_DATAWIDTH-1, 0) = output_ifc[output_address].range(1 * OACTS_DATAWIDTH -1, 0);
                output[i][j+1].range(OACTS_DATAWIDTH-1, 0) = output_ifc[output_address].range(2 * OACTS_DATAWIDTH -1,  1*OACTS_DATAWIDTH);
                output[i][j+2].range(OACTS_DATAWIDTH-1, 0) = output_ifc[output_address].range(3 * OACTS_DATAWIDTH -1, 2*OACTS_DATAWIDTH);
                output[i][j+3].range(OACTS_DATAWIDTH-1, 0) = output_ifc[output_address].range(4 * OACTS_DATAWIDTH -1, 3*OACTS_DATAWIDTH);
                //overall_addr++;
                output_address++;
                }
            }
        }
    }


    std::cout << "Print Out OUTPUT Data" << std::endl;
    for (int i = 0; i < X; ++i)
    {
        for (int j = 0; j < Wt_Y; ++j)
        {
            std::cout << output[i][j] << ", ";
            if (j%20 == 0) std::cout << std::endl;
        }
        // std::cout << std::endl;
    }
    long double mse = 0.0;
    //Compute MSE
    for (int i = 0; i < X; ++i)
    {
        for (int j = 0; j < Wt_Y; ++j)
        {
            mse += std::pow((output[i][j]
                    - (data_t)oacts[i+j]), 2.0);
        }
    }
    
    mse = mse / (X*Wt_Y);

    std::cout << "Output MSE:  " << mse << std::endl;
    return 0;  
}

//     for (int i = 0; i < value_loop_count_X; i++) { //64
//         for (int j = 0; j < value_loop_count_Y; j++) { //50
//             for (int row = 0; row < PARALLEL_K; ++row) { //32 
//                 for (int col = 0; col < PARALLEL_N / value_parallel_transmission; ++col) { //20/4 = 5  50 is not divisible
//                     ap_uint<HP_IFC_BANDWIDTH> temp = 0;
//                     for (int idx = 0; idx < value_parallel_transmission; ++idx) {
//                         int idx_x = i*PARALLEL_K+row;
//                         int idx_y = j*PARALLEL_N+col*value_parallel_transmission+idx;
//                         temp.range((idx+1)*VALUE_DATAWIDTH-1, idx*VALUE_DATAWIDTH) = value[idx_x*Wt_Y+idx_y].range(31, 0);

// #ifdef DEBUG_MAIN
//                         std::cout << "idx_x: " << idx_x << " idx_y" << idx_y << " temp:" << temp << std::endl;
//                         std::cout << value[idx_x*Wt_Y+idx_y] << ", " << std::endl;
// #endif

//                     }
//                     ifc[overall_addr + count] = temp.range(HP_IFC_BANDWIDTH-1, 0);
//                     count++;
//                     weight_addr_offset++;

// #ifdef DEBUG_MAIN
//                     std::cout << "ifc: " << ifc[overall_addr + count] << std::endl;
//                     std::cout << "Count" << count << std::endl;

// #endif
//                 }
//             }   
//         }
//     }

//     if (weight_loop_residual != 0) {
//     //address the residual elements (Should not trigger the following part)
//         std::cout << "Weight residual happens!! Check what happens" << std::endl;
//         ap_uint<HP_IFC_BANDWIDTH> temp;
//         for (int i = 0; i < weight_parallel_transmission; ++i) {
//             if (i < weight_loop_residual) {
//                 int idx_x = (weight_complete_loop_count*weight_parallel_transmission+i) / Wt_X;
//                 int idx_y = (weight_complete_loop_count*weight_parallel_transmission+i) % Wt_X;
//                 temp.range(i*value_DATAWIDTH, (i+1)*value_DATAWIDTH-1) = value[idx_x*Wt_Y+idx_y].range(31, 0);
//             }
//             else
//             {
//                 temp.range(i*value_DATAWIDTH, (i+1)*value_DATAWIDTH-1) = 0;
//             }
//         }
//         ifc[overall_addr] = temp;
//         weight_addr_offset++;
//     }

//     overall_addr += weight_addr_offset;

//     //Layout Iact data
//     std::cout << "Start layout input activation data" << std::endl;
//     int iact_parallel_transmission = HP_IFC_BANDWIDTH / attention_DATAWIDTH; //128 / 32 = 4
//     int iact_complete_loop_count = X * Y / iact_parallel_transmission; //2048 / 4
//     int iact_loop_residual = X * Y % iact_parallel_transmission; 
//     int iact_addr_offset = 0;
//     for (int i = 0; i < iact_complete_loop_count; ++i) {
//         ap_uint<HP_IFC_BANDWIDTH> temp;
//         for (int j = 0; j < iact_parallel_transmission; ++j) {
//             int idx = i*iact_parallel_transmission+j;
//             //int idx_y = (i*iact_parallel_transmission+j) % X;
//             temp.range((j+1)*attention_DATAWIDTH-1, j*attention_DATAWIDTH) = attention[idx].range(31, 0); //QUESTION: BIG_ENDIAN OR SMALL_ENGIAN READ
//             //std::cout << "idx: " << idx << " temp:" << temp << std::endl;
//         }
//         ifc[overall_addr + i] = temp;
//         iact_addr_offset++;
//     }
//     if (iact_loop_residual != 0) {
//         //address the residual elements  -> Not triggered this part for now
//         std::cout << "IACT residual happens!! Check what happens" << std::endl;
//         for (int i = 0; i < iact_parallel_transmission; ++i) {
//             ap_uint<HP_IFC_BANDWIDTH> temp;
//             if (i < iact_loop_residual) {
//                 int idx = iact_complete_loop_count * iact_parallel_transmission + i;
//                 // int idx_x = (iact_complete_loop_count*iact_parallel_transmission+i) / X;
//                 // int idx_y = (iact_complete_loop_count*iact_parallel_transmission+i) % X;
//                 temp.range(i*attention_DATAWIDTH, (i+1)*attention_DATAWIDTH-1) = attention[idx].range(31, 0);
//             }
//             else
//             {
//                 temp.range(i*attention_DATAWIDTH, (i+1)*attention_DATAWIDTH-1) = 0;
//             }
//             ifc[overall_addr + i] = temp;
//             iact_addr_offset++;
//         }

//     }

//     overall_addr += iact_addr_offset;
