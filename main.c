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
    algos[0] = (Algo){
        .name = "Dummy",
        .init = dummy_heartrate_init,
        .get_heartrate = dummy_heartrate,
        .total_time = 0,
    };

    algos[1] = (Algo){
        .name = "Espruino",
        .init = espruino_heartrate_init,
        .get_heartrate = espruino_heartrate,
        .total_time = 0,
    };

    algos[2] = (Algo){
        .name = "FFT",
        .init = fft_heartrate_init,
        .get_heartrate = fft_heartrate,
        .total_time = 0,
    };

    algos[3] = (Algo){
        .name = "Autocorrelation",
        .init = autocorrelation_heartrate_init,
        .get_heartrate = autocorrelation_heartrate,
        .total_time = 0,
    };

    algos[4] = (Algo){
        .name = "Oxford",
        .init = oxford_heartrate_init,
        .get_heartrate = oxford_heartrate,
        .total_time = 0,
    };

    algos[5] = (Algo){
        .name = "PanTompkins",
        .init = pantompkins_heartrate_init,
        .get_heartrate = pantompkins_heartrate,
        .total_time = 0,
    };
    algos[6] = (Algo){
        .name = "Autocorrelation2",
        .init = autocorrelation2_heartrate_init,
        .get_heartrate = autocorrelation2_heartrate,
        .total_time = 0,
    };
    algos[7] = (Algo){
        .name = "Algo1",
        .init = algo1_heartrate_init,
        .get_heartrate = algo1_heartrate,
        .total_time = 0,
    };
    algos[8] = (Algo){
        .name = "FFT2",
        .init = fft2_heartrate_init,
        .get_heartrate = fft2_heartrate,
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
        for (int i = 0; i < algoN; i++)
        {
            fprintf(out_fp, "%s", algos[i].name);
            if (i < algoN - 1)
                fprintf(out_fp, ",");
        }
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

                // Parse integer values using sscanf
                if (sscanf(line, "%d,%d,%d,%d,%d", &ms, &ppg, &accx, &accy, &accz) != 5) //c'era un 4
                {
                    printf("Error parsing line: %s\n", line);
                    continue;
                }

                // if (lineN < 10 && ms > 5000)
                // {
                //     continue;
                // }

                // Process the extracted integer values
                printf("Values: %d, %d, %d, %d, %d\n", ms, ppg, accx, accy, accz);
                printf("Debug688\n");
                

                int delta_ms = 0;
                if (lineN > 2)
                    delta_ms = ms - previous_ms;

                fprintf(out_fp, "%d,", ms);

                
                // call all algorithms here:
                for (int i = 0; i < algoN; i++)
                {
                    clock_t start_t;
                    start_t = clock();
                    int hr = algos[i].get_heartrate(delta_ms, ppg, accx, accy, accz);
                    algos[i].total_time += clock() - start_t;
                    // write output file
                    fprintf(out_fp, "%d", hr);
                    if (i < algoN - 1)
                    {
                        fprintf(out_fp, ",");
                    }
                }
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