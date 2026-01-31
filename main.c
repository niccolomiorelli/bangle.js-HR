#include <stdio.h>
//#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "types.h"
#include "dirent.h" //If dirent.h is not pre included (in Microsoft it could not be included)
#include "algorithms/dummy/dummyHeartRate.h"
#include "algorithms/espruino/espruinoHeartRate.h"
#include "algorithms/fft/fftHeartRate.h"
#include "algorithms/autocorrelation/autocorrelationHeartRate.h"
#include "algorithms/oxford/oxfordHeartRate.h"
#include "algorithms/pantompkins/pantompkinsHeartRate.h"
#include "algorithms/autocorrelation2/autocorrelation2HeartRate.h"
#include "algorithms/algo1/algo1HeartRate.h"
#include "algorithms/fft2/fft2HeartRate.h"
#include "algorithms/spectralTracking/spectralTrackingHeartRate.h"
#include "algorithms/final/finalHeartRate.h"
#include "algorithms/finalOpt/finalOptHeartRate.h"
#include "algorithms/finalShort/finalShortHeartRate.h"
#include "algorithms/debug1/debug1HeartRate.h"
#include "algorithms/debug2/debug2HeartRate.h"
#include "algorithms/debug3/debug3HeartRate.h"
#include "algorithms/debug4/debug4HeartRate.h"
#include "algorithms/debug5/debug5HeartRate.h"
#include "algorithms/debug6/debug6HeartRate.h"
#include "algorithms/debug7/debug7HeartRate.h"
#include "algorithms/debug8/debug8HeartRate.h"
#include "algorithms/debug9/debug9HeartRate.h"
#include "algorithms/debug10/debug10HeartRate.h"
#include "algorithms/debug11/debug11HeartRate.h"
#include "algorithms/debug88/debug88HeartRate.h"
#include "algorithms/debug89/debug89HeartRate.h"
#include "algorithms/debug82/debug82HeartRate.h"
#include "algorithms/debug882/debug882HeartRate.h"
#include "algorithms/debug892/debug892HeartRate.h"



//Aggiungo
#define algoN 9

typedef struct Algo
{
    char *name;
    void (*init)();
    int (*get_heartrate)(time_delta_ms_t delta_ms, ppg_t ppg, accel_t accx, accel_t accy, accel_t accz);
    clock_t total_time;
} Algo;

////////////////////////////////////
// START MODIFY HERE TO ADD NEW ALGO


// all algorithms:
//const int algoN = 1; // change this to algo number!
Algo algos[algoN];

void createAlgos()
{
    // algos[0] = (Algo){
    //     .name = "Dummy",
    //     .init = dummy_heartrate_init,
    //     .get_heartrate = dummy_heartrate,
    //     .total_time = 0,
    // };

    // algos[1] = (Algo){
    //     .name = "Espruino",
    //     .init = espruino_heartrate_init,
    //     .get_heartrate = espruino_heartrate,
    //     .total_time = 0,
    // };

    // algos[2] = (Algo){
    //     .name = "FFT",
    //     .init = fft_heartrate_init,
    //     .get_heartrate = fft_heartrate,
    //     .total_time = 0,
    // };

    // algos[2] = (Algo){
    //     .name = "Autocorrelation",
    //     .init = autocorrelation_heartrate_init,
    //     .get_heartrate = autocorrelation_heartrate,
    //     .total_time = 0,
    // };

    // algos[4] = (Algo){
    //     .name = "Oxford",
    //     .init = oxford_heartrate_init,
    //     .get_heartrate = oxford_heartrate,
    //     .total_time = 0,
    // };

    // algos[5] = (Algo){
    //     .name = "PanTompkins",
    //     .init = pantompkins_heartrate_init,
    //     .get_heartrate = pantompkins_heartrate,
    //     .total_time = 0,
    // };
    // algos[6] = (Algo){
    //     .name = "Autocorrelation2",
    //     .init = autocorrelation2_heartrate_init,
    //     .get_heartrate = autocorrelation2_heartrate,
    //     .total_time = 0,
    // };
    // algos[7] = (Algo){
    //     .name = "Algo1",
    //     .init = algo1_heartrate_init,
    //     .get_heartrate = algo1_heartrate,
    //     .total_time = 0,
    // };
    algos[0] = (Algo){
        .name = "FFT2",
        .init = fft2_heartrate_init,
        .get_heartrate = fft2_heartrate,
        .total_time = 0,
    };
    // algos[9] = (Algo){
    //     .name = "SpectralTracking",
    //     .init = spectralTracking_heartrate_init,
    //     .get_heartrate = spectralTracking_heartrate,
    //     .total_time = 0,
    // };
    // algos[10] = (Algo){
    //     .name = "Final",
    //     .init = final_heartrate_init,
    //     .get_heartrate = final_heartrate,
    //     .total_time = 0,
    // };
    algos[1] = (Algo){
        .name = "FinalOpt",
        .init = finalOpt_heartrate_init,
        .get_heartrate = finalOpt_heartrate,
        .total_time = 0,
    };
    // algos[12] = (Algo){
    //     .name = "FinalShort",
    //     .init = finalShort_heartrate_init,
    //     .get_heartrate = finalShort_heartrate,
    //     .total_time = 0,
    // };
    algos[2] = (Algo){
         .name = "Debug1",
         .init = debug1_heartrate_init,
         .get_heartrate = debug1_heartrate,
         .total_time = 0,
    };
    // algos[6] = (Algo){
    //      .name = "Debug2",
    //      .init = debug2_heartrate_init,
    //      .get_heartrate = debug2_heartrate,
    //      .total_time = 0,
    // };
    // algos[7] = (Algo){
    //      .name = "Debug3",
    //      .init = debug3_heartrate_init,
    //      .get_heartrate = debug3_heartrate,
    //      .total_time = 0,
    // };
    // algos[8] = (Algo){
    //      .name = "Debug4",
    //      .init = debug4_heartrate_init,
    //      .get_heartrate = debug4_heartrate,
    //      .total_time = 0,
    // };
    // algos[9] = (Algo){
    //      .name = "Debug5",
    //      .init = debug5_heartrate_init,
    //      .get_heartrate = debug5_heartrate,
    //      .total_time = 0,
    // };
    // algos[10] = (Algo){
    //      .name = "Debug6",
    //      .init = debug6_heartrate_init,
    //      .get_heartrate = debug6_heartrate,
    //      .total_time = 0,
    // };
    // algos[3] = (Algo){
    //      .name = "Debug7",
    //      .init = debug7_heartrate_init,
    //      .get_heartrate = debug7_heartrate,
    //      .total_time = 0,
    // };
    algos[3] = (Algo){
         .name = "Version3_2",
         .init = debug8_heartrate_init,
         .get_heartrate = debug8_heartrate,
         .total_time = 0,
    };
    // algos[5] = (Algo){
    //      .name = "Debug9",
    //      .init = debug9_heartrate_init,
    //      .get_heartrate = debug9_heartrate,
    //      .total_time = 0,
    // };
    // algos[5] = (Algo){
    //      .name = "Debug10",
    //      .init = debug10_heartrate_init,
    //      .get_heartrate = debug10_heartrate,
    //      .total_time = 0,
    // };
    // algos[6] = (Algo){
    //      .name = "Debug11",
    //      .init = debug11_heartrate_init,
    //      .get_heartrate = debug11_heartrate,
    //      .total_time = 0,
    // };
    algos[4] = (Algo){
         .name = "Version2_2",
         .init = debug88_heartrate_init,
         .get_heartrate = debug88_heartrate,
         .total_time = 0,
    };
    algos[5] = (Algo){
         .name = "Version1_2",
         .init = debug89_heartrate_init,
         .get_heartrate = debug89_heartrate,
         .total_time = 0,
    };
    algos[6] = (Algo){
         .name = "Version3_1",
         .init = debug82_heartrate_init,
         .get_heartrate = debug82_heartrate,
         .total_time = 0,
    };
    algos[7] = (Algo){
         .name = "Version2_1",
         .init = debug882_heartrate_init,
         .get_heartrate = debug882_heartrate,
         .total_time = 0,
    };
    algos[8] = (Algo){
         .name = "Version1_1",
         .init = debug892_heartrate_init,
         .get_heartrate = debug892_heartrate,
         .total_time = 0,
    };

}

int main(int argc, char *argv[])
{

    
    if (argc < 3)
    {
        printf("Usage: %s <input directory> <output directory>\n", argv[0]);
        return 1;
    }

    DIR *input_dir;
    struct dirent *entry;
    // open the input directory
    input_dir = opendir(argv[1]);
    if (input_dir == NULL)
    {
        perror("Cannot open input directory");
        return 1;
    }

    createAlgos();

    // input ppg file
    FILE *ppg_fp;
    char line[1024];

    // read each file
    while ((entry = readdir(input_dir)) != NULL)
    {
        if (strstr(entry->d_name, ".csv") == 0)
        {
            // exclude non csv files
            continue;
        }

        // add path to the filename
        char filename[256];
        snprintf(filename, sizeof(filename), "%s/%s", argv[1], entry->d_name);

        // open the file
        ppg_fp = fopen(filename, "r");
        if (ppg_fp == NULL)
        {
            perror("Cannot open input file");
            printf("%s", filename);
            continue;
        }

        // initialise the algorithms
        for (int i = 0; i < algoN; i++)
        {
            algos[i].total_time = 0;
            algos[i].init();
        }

        // open the output file
        // Construct the full path to the file
        char hr_filepath[256];
        snprintf(hr_filepath, sizeof(hr_filepath), "%s/HR_%s", argv[2], entry->d_name);


        // Open the file for writing
        FILE *out_fp = fopen(hr_filepath, "w");
        if (out_fp == NULL)
        {
            perror("Error opening output file");
            return 1;
        }
        // write header of output file
        fprintf(out_fp, "time,");
        fprintf(out_fp, "time (abs),");
        for (int i = 0; i < algoN; i++)
        {
            fprintf(out_fp, "%s", algos[i].name);
            fprintf(out_fp, ",");
        }
        fprintf(out_fp, "GT_Bangle,");
        fprintf(out_fp, "GT_polar,");
        // fprintf(out_fp,"GT_cosmed");
        fprintf(out_fp, "\n");

        // counter of the line number
        unsigned int lineN = 0;
        unsigned int previous_ms = 0;

        while (fgets(line, sizeof(line), ppg_fp) != NULL)
        {
            // Process the line here:
            // printf("%s", line);
            //printf("Test2\n");


            lineN++;
            // discard first line, used for header
            if (lineN > 1)
            {
                // Remove trailing newline
                line[strcspn(line, "\n")] = 0;

                int ms, ppg, accx, accy, accz;

                int GT_bangle, GT_polar, GT_cosmed;
                long long ms_abs;

                // Parse integer values using sscanf
                if (sscanf(line, "%d,%lld,%d,%d,%d,%d,%d,%d", &ms, &ms_abs, &ppg, &accx, &accy, &accz, &GT_bangle, &GT_polar) != 8) //c'era un 5
                {
                    printf("Error parsing line: %s\n", line);
                    continue;
                }

                // if (lineN < 10 && ms > 5000)
                // {
                //     continue;
                // }

                // Process the extracted integer values
                printf("Values: %d, %lld, %d, %d, %d, %d, %d, %d\n", ms, ms_abs, ppg, accx, accy, accz, GT_bangle, GT_polar);
                printf("Debug\n");
                

                int delta_ms = 0;
                if (lineN > 2)
                    delta_ms = ms - previous_ms;

                fprintf(out_fp, "%d,", ms);
                fprintf(out_fp, "%lld,",ms_abs);

                
                // call all algorithms here:
                for (int i = 0; i < algoN; i++)
                {
                    clock_t start_t;
                    start_t = clock();
                    int hr = algos[i].get_heartrate(delta_ms, ppg, accx, accy, accz);
                    algos[i].total_time += clock() - start_t;
                    // write output file
                    fprintf(out_fp, "%d", hr);
                    fprintf(out_fp, ",");
                    
                }
                //Printing also the Ground Truths 
                fprintf(out_fp, "%d", GT_bangle);
                fprintf(out_fp, ",");
                fprintf(out_fp, "%d", GT_polar);
                // fprintf(out_fp, ",");
                // fprintf(out_fp, "%d", GT_cosmed);
                fprintf(out_fp, "\n");

                previous_ms = ms;
            }
        }
        fclose(ppg_fp);
        fclose(out_fp);
    }

    closedir(input_dir);

    return 0;
}