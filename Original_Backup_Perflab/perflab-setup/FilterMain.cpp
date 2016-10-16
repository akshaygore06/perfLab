/*
* Program  : Performance Lab
* Programmer : Akshay Gore
* CSUC Id : 006951205 
* References : https://en.wikipedia.org/wiki/Loop_unrolling
               http://www.openmp.org/mp-documents/OpenMP3.0-SummarySpec.pdf
               https://msdn.microsoft.com/en-us/library/0ca2w8dk.aspx
               https://www.ibm.com/support/knowledgecenter/SSGH2K_13.1.0/com.ibm.xlc131.aix.doc/compiler_ref/prag_omp_parallel.html
*/


#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <stdint.h>
#include <fstream>
#include "Filter.h"
#include <omp.h>
#include "rtdsc.h"
#include <string>

using namespace std;

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string FilteredOutput = filtername;
  size_t loc = FilteredOutput.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    FilteredOutput = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + FilteredOutput + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
	int value;
	input >> value;
	filter -> set(i,j,value);
      }
    }
    return filter;
  }
}

#if defined(__arm__)
static inline unsigned int get_cyclecount (void)
{
 unsigned int value;
 // Read CCNT Register
 asm volatile ("MRC p15, 0, %0, c9, c13, 0\t\n": "=r"(value));
 return value;
}

static inline void init_perfcounters (int32_t do_reset, int32_t enable_divider)
{
 // in general enable all counters (including cycle counter)
 int32_t value = 1;

 // peform reset:
 if (do_reset)
 {
   value |= 2;     // reset all counters to zero.
   value |= 4;     // reset cycle counter to zero.
 }

 if (enable_divider)
   value |= 8;     // enable "by 64" divider for CCNT.

 value |= 16;

 // program the performance-counter control-register:
 asm volatile ("MCR p15, 0, %0, c9, c12, 0\t\n" :: "r"(value));

 // enable all counters:
 asm volatile ("MCR p15, 0, %0, c9, c12, 1\t\n" :: "r"(0x8000000f));

 // clear overflows:
 asm volatile ("MCR p15, 0, %0, c9, c12, 3\t\n" :: "r"(0x8000000f));
}



#endif


double
applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output)
{
  #if defined(__arm__)
  init_perfcounters (1, 1);
  #endif

  long long cycStart, cycStop;
  #if defined(__arm__)
  cycStart = get_cyclecount();
  #else
  cycStart = rdtscll();
  #endif
  output->width  = input->width;
  output->height = input->height;

  int col, row, plane, finalOutput;
  int filterDivisor = filter->getDivisor();
  int inputWidth    = input->width - 1;
  int inputHeight   = input->height - 1;

int filterType = 0;
int filterValue = filter->get(2,1); // taking value of pixel(2,1)

if( filterValue == 4) // applied filter is  GAUSS FILTER
{
  filterType = 1;
}
else if( filterValue == -1) // applied filter is  EMBOSS FILTER
{
  filterType = 2;
}
else if( filterValue == 1) // applied filter is  AVERAGE FILTER
{
  filterType = 3;
}
else // applied filter is  HLINE FILTER
{
  filterType = 4;
}   

  switch (filterType) {
    case 1:
#pragma omp parallel for private(plane, row, col)  // creating private threads for parallel processing
    for (plane = 0; plane < 3; ++plane) {
      for (row = 1; row < inputHeight; ++row) {
        for (col = 1; col < inputWidth; ++col) {

          output->color[plane][row][col] = 0;
          finalOutput = output->color[plane][row][col];
          // Applying gauss filter
          finalOutput = (input->color[plane][row-1][col] << 2) +
                                   (input->color[plane][row][col-1] << 2) +
                                   (input->color[plane][row][col] << 3) +
                                   (input->color[plane][row][col+1] << 2) +
                                   (input->color[plane][row+1][col] << 2);

          finalOutput = finalOutput / filterDivisor;

          if (finalOutput < 0) {
            finalOutput = 0;
          } else if (finalOutput > 255) {
            finalOutput = 255;
          }
          output->color[plane][row][col] = finalOutput;
        }
      }
    }
    break;
    case 2:
  #pragma omp parallel for private(plane, row, col) // creating private threads for parallel processing
    for (plane = 0; plane < 3; ++plane) {
      for (row = 1; row < inputHeight; ++row) {
        for (col = 1; col < inputWidth; ++col) {

          output->color[plane][row][col] = 0;
          finalOutput = output->color[plane][row][col];
          // Applying emboss filter
          finalOutput = (input->color[plane][row-1][col-1] << 0) +
                                   (input->color[plane][row-1][col] << 0) -
                                   (input->color[plane][row-1][col+1] <<0) +
                                   (input->color[plane][row][col-1] <<0)+
                                   (input->color[plane][row][col]<<0) -
                                   (input->color[plane][row][col+1] <<0)+
                                   (input->color[plane][row+1][col-1] <<0)-
                                   (input->color[plane][row+1][col] <<0) -
                                   (input->color[plane][row+1][col+1] <<0);

          if (finalOutput < 0) {
            finalOutput = 0;
          } else if (finalOutput > 255) {
            finalOutput = 255;
          }
          output->color[plane][row][col] = finalOutput;
        }
      }
    }
    break;
    case 3:
#pragma omp parallel for private(plane, row, col) // creating private threads for parallel processing
    for (plane = 0; plane < 3; ++plane) {
      for (row = 1; row < inputHeight; ++row) {
        for (col = 1; col < inputWidth; ++col) {

          output->color[plane][row][col] = 0;
          finalOutput = output->color[plane][row][col];
          // Applying avg filter
          finalOutput = input->color[plane][row-1][col-1] +
                                   input->color[plane][row-1][col] +
                                   input->color[plane][row-1][col+1] +
                                   input->color[plane][row][col-1] +
                                   input->color[plane][row][col] +
                                   input->color[plane][row][col+1] +
                                   input->color[plane][row+1][col-1] +
                                   input->color[plane][row+1][col] +
                                   input->color[plane][row+1][col+1];

          finalOutput /= filterDivisor;

          if (finalOutput < 0) {
            finalOutput = 0;
          } else if (finalOutput > 255) {
            finalOutput = 255;
          }  


         output->color[plane][row][col] = finalOutput;

        }
      }
    }
    break;
    case 4:
#pragma omp parallel for private(plane, row, col)  // creating private threads for parallel processing
    for (plane = 0; plane < 3; ++plane) {
      for (row = 1; row < inputHeight; ++row) {
        for (col = 1; col < inputWidth; ++col) {

          output->color[plane][row][col] = 0;
          finalOutput = output->color[plane][row][col];
          // Applying HLINE filter
          finalOutput = -input->color[plane][row-1][col-1] -
                                   (input->color[plane][row-1][col] << 1) -
                                   input->color[plane][row-1][col+1] +
                                   input->color[plane][row+1][col-1] +
                                   (input->color[plane][row+1][col] << 1) +
                                   input->color[plane][row+1][col+1];

          if (finalOutput < 0) {
            finalOutput = 0;
          } else if (finalOutput > 255) {
            finalOutput = 255;
          }
          output->color[plane][row][col] = finalOutput;
        }
      }
    }
    break;
  }

  #if defined(__arm__)
  cycStop = get_cyclecount();
  #else
  cycStop = rdtscll();
  #endif

  double diff = cycStop-cycStart;
  #if defined(__arm__)
  diff = diff * 64;
  #endif
  double diffPerPixel = diff / (output -> width * output -> height);
  fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
        diff, diff / (output -> width * output -> height));
  return diffPerPixel ;
}
